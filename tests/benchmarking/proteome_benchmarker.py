import os
import subprocess
import traceback
import csv
import time
from genedesign.seq_utils.Translate import Translate
from genedesign.transcript_designer import TranscriptDesigner
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker

def parse_fasta(fasta_file):
    """
    Parses the FASTA file to extract gene names and protein sequences.
    """
    sequences = {}
    current_gene = None
    current_sequence = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_gene:
                    sequences[current_gene] = ''.join(current_sequence)
                gene_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith("GN="):
                        gene_name = part.split("=")[1]
                        break
                if not gene_name:
                    gene_name = line.split('|')[2].split(' ')[0]
                current_gene = gene_name
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_gene:
            sequences[current_gene] = ''.join(current_sequence)
    
    return sequences

def benchmark_proteome(fasta_file):
    """
    Benchmarks the proteome using TranscriptDesigner.
    """
    designer = TranscriptDesigner()
    designer.initiate()

    proteome = parse_fasta(fasta_file)
    successful_results = []
    error_results = []

    for gene, protein in proteome.items():
        try:
            print(f"Processing gene: {gene} with protein sequence: {protein[:30]}...")
            ignores = set()
            transcript = designer.run(protein, ignores)
            successful_results.append({
                'gene': gene,
                'protein': protein,
                'transcript': transcript
            })
        except Exception as e:
            error_results.append({
                'gene': gene,
                'protein': protein,
                'error': f"Error: {str(e)}\nTraceback: {traceback.format_exc()}"
            })
    
    return successful_results, error_results

def analyze_errors(error_results):
    """
    Write the error analysis to a text file.
    """
    error_summary = {}
    with open('error_summary.txt', 'w') as f:
        for error in error_results:
            error_message = error['error'].split("\n")[0]
            error_summary[error_message] = error_summary.get(error_message, 0) + 1
            f.write(f"Gene: {error['gene']}\n{error['error']}\n\n")
    
    return error_summary

def validate_transcripts(successful_results):
    """
    Validate the successful transcripts using various checkers.
    Returns (validation_failures, checker_passes) where checker_passes maps
    checker name to count of sequences that passed. Sequences that crash the
    designer never enter here, so they are implicit failures for all checks.
    """
    forbidden_checker = ForbiddenSequenceChecker()
    forbidden_checker.initiate()
    promoter_checker = PromoterChecker()
    promoter_checker.initiate()
    translator = Translate()
    translator.initiate()
    codon_checker = CodonChecker()
    codon_checker.initiate()

    validation_failures = []
    checker_passes = {
        'translation': 0,
        'forbidden': 0,
        'hairpin': 0,
        'promoter': 0,
        'codon': 0,
    }

    for result in successful_results:
        cds = ''.join(result['transcript'].codons)
        try:
            if len(cds) % 3 != 0:
                raise ValueError("CDS length is not a multiple of 3.")
            original_protein = result['protein']
            translated_protein = translator.run(cds)
            if original_protein != translated_protein:
                raise ValueError(f"Translation mismatch: Original {original_protein}, Translated {translated_protein}")
            if not (cds.startswith(("ATG", "GTG", "TTG")) and cds.endswith(("TAA", "TGA", "TAG"))):
                raise ValueError("CDS does not start with a valid start codon or end with a valid stop codon.")
        except ValueError as e:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': cds,
                'site': f"Translation or completeness error: {str(e)}"
            })
            continue

        checker_passes['translation'] += 1

        transcript_dna = result['transcript'].rbs.utr.upper() + cds

        passed_hairpin, hairpin_string = hairpin_checker(transcript_dna)
        if passed_hairpin:
            checker_passes['hairpin'] += 1
        else:
            formatted_hairpin = hairpin_string.replace('\n', ' ').replace('"', "'")
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Hairpin detected: {formatted_hairpin}"
            })

        passed_forbidden, forbidden_site = forbidden_checker.run(transcript_dna)
        if passed_forbidden:
            checker_passes['forbidden'] += 1
        else:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Forbidden sequence: {forbidden_site}"
            })

        passed_promoter, found_promoter = promoter_checker.run(transcript_dna)
        if passed_promoter:
            checker_passes['promoter'] += 1
        else:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Constitutive promoter detected: {found_promoter}" if found_promoter else "Constitutive promoter detected"
            })

        codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(result['transcript'].codons)
        if codons_above_board:
            checker_passes['codon'] += 1
        else:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': cds,
                'site': f"Codon usage check failed: Diversity={codon_diversity}, Rare Codons={rare_codon_count}, CAI={cai_value}"
            })

    return validation_failures, checker_passes

def write_validation_report(validation_failures):
    """
    Writes validation results to a TSV file.
    """
    with open('validation_failures.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['gene', 'protein', 'cds', 'site'])
        for failure in validation_failures:
            writer.writerow([failure['gene'], failure['protein'], failure['cds'], failure['site']])

def find_student_files(repo_root):
    """
    Detects the GitHub repo URL, current branch, and any files students added
    beyond the original scaffold — new checkers and new tests.
    Returns a dict with keys: base_url, branch, transcript_designer_url,
    new_checkers (list of (filename, url)), new_tests (list of (filename, url)).
    On any failure, returns None.
    """
    ORIGINAL_CHECKERS = {
        'codon_checker.py',
        'forbidden_sequence_checker.py',
        'hairpin_checker.py',
        'internal_promoter_checker.py',
    }
    ORIGINAL_TESTS = {
        'test_codon_checker.py',
        'test_forbidden_sequence_checker.py',
        'test_internal_promoter_checker.py',
        'test_operon_designer.py',
        'test_transcript_designer.py',
        'test_hairpin_counter.py',
    }

    try:
        raw = subprocess.check_output(
            ['git', 'remote', 'get-url', 'origin'],
            cwd=repo_root, stderr=subprocess.DEVNULL
        ).decode().strip()
        # Normalise SSH → HTTPS and strip .git
        if raw.startswith('git@github.com:'):
            raw = 'https://github.com/' + raw[len('git@github.com:'):]
        base_url = raw.removesuffix('.git')

        branch = subprocess.check_output(
            ['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
            cwd=repo_root, stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return None

    def gh(rel_path):
        return f"{base_url}/blob/{branch}/{rel_path}"

    # TranscriptDesigner is always at the same place
    td_path = 'genedesign/transcript_designer.py'
    result = {
        'base_url': base_url,
        'branch': branch,
        'transcript_designer_url': gh(td_path),
        'new_checkers': [],
        'new_tests': [],
    }

    # Scan genedesign/checkers/ for anything not in the original set
    checker_dir = os.path.join(repo_root, 'genedesign', 'checkers')
    if os.path.isdir(checker_dir):
        for fname in sorted(os.listdir(checker_dir)):
            if fname.endswith('.py') and fname not in ORIGINAL_CHECKERS and fname != '__init__.py':
                rel = f'genedesign/checkers/{fname}'
                result['new_checkers'].append((fname, gh(rel)))

    # Scan the whole tests tree for anything not in the original set
    tests_dir = os.path.join(repo_root, 'tests')
    if os.path.isdir(tests_dir):
        for dirpath, _, filenames in os.walk(tests_dir):
            for fname in sorted(filenames):
                if fname.endswith('.py') and fname not in ORIGINAL_TESTS and fname != '__init__.py':
                    # Skip the benchmarker itself
                    if fname == 'proteome_benchmarker.py':
                        continue
                    rel = os.path.relpath(os.path.join(dirpath, fname), repo_root).replace(os.sep, '/')
                    result['new_tests'].append((fname, gh(rel)))

    return result


def generate_grade_report(total_genes, checker_passes, parsing_time, execution_time, proteome_name, student_files=None):
    """
    Generates a grade report with per-checker pass rates, thresholds, and a point score.
    Prints the report and writes it to grade_report.txt.

    Scoring table (55 pts total):
      20 pts  Valid ORFs produced          >= 80%
      20 pts  Forbidden sequence checker   >= 90%
       5 pts  Codon checker (CAI/bias)     >= 90%
       5 pts  Promoter checker             >= 90%
       2 pts  Hairpin checker              >= 90%
       2 pts  [6th check — TBD]            >= 90%
       1 pt   All of the above             >= 99%
    """
    THRESHOLDS = {
        'translation': 0.80,
        'forbidden':   0.90,
        'codon':       0.90,
        'promoter':    0.90,
        'hairpin':     0.90,
    }
    MAX_POINTS = {
        'translation': 20,
        'forbidden':   20,
        'codon':        5,
        'promoter':     5,
        'hairpin':      2,
        'sixth':        2,
        'all99':        1,
    }

    rates = {key: checker_passes[key] / total_genes for key in THRESHOLDS}

    score = sum(MAX_POINTS[k] for k in THRESHOLDS if rates[k] >= THRESHOLDS[k])

    # Bonus point if every defined check is >= 99%
    all_above_99 = all(rates[k] >= 0.99 for k in THRESHOLDS)
    if all_above_99:
        score += MAX_POINTS['all99']
    all99_pts = MAX_POINTS['all99'] if all_above_99 else 0

    def _row(label, n, rate, threshold_pct, max_pts):
        pts = max_pts if rate >= threshold_pct / 100 else 0
        if pts > 0:
            return f"[PASS] {label}: {n}/{total_genes} ({rate*100:.1f}%) — {pts}/{max_pts} pts"
        else:
            return (f"[FAIL] {label}: {n}/{total_genes} ({rate*100:.1f}%) — 0/{max_pts} pts"
                    f" — needs >={threshold_pct}%, {total_genes - n} failures")

    total_runtime = parsing_time + execution_time
    lines = [
        "TranscriptDesigner Benchmarker Grade Report",
        f"Proteome: {proteome_name} ({total_genes} proteins)",
        f"Runtime: {total_runtime:.1f} seconds",
        "",
        f"SCORE: {score} / 55",
        "",
        _row("Valid ORFs",          checker_passes['translation'], rates['translation'], 80,  MAX_POINTS['translation']),
        _row("Forbidden sequences", checker_passes['forbidden'],   rates['forbidden'],   90,  MAX_POINTS['forbidden']),
        _row("Codon bias",          checker_passes['codon'],       rates['codon'],       90,  MAX_POINTS['codon']),
        _row("Promoters",           checker_passes['promoter'],    rates['promoter'],    90,  MAX_POINTS['promoter']),
        _row("Hairpins",            checker_passes['hairpin'],     rates['hairpin'],     90,  MAX_POINTS['hairpin']),
        "[ N/A] 6th check (TBD) — 0/2 pts",
        f"[{'PASS' if all_above_99 else 'FAIL'}] All checks >=99% bonus — {all99_pts}/1 pts",
        "",
        "--- Code Links ---",
    ]

    if student_files is None:
        lines.append("[WARNING] Could not detect GitHub repo — add your repo URL here manually.")
    else:
        lines.append(f"TranscriptDesigner: {student_files['transcript_designer_url']}")
        if student_files['new_checkers']:
            for fname, url in student_files['new_checkers']:
                lines.append(f"New checker — {fname}: {url}")
        else:
            lines.append("[WARNING] No new checker found in genedesign/checkers/")
        if student_files['new_tests']:
            for fname, url in student_files['new_tests']:
                lines.append(f"New test — {fname}: {url}")
        else:
            lines.append("[WARNING] No new test file found in tests/")

    report = "\n".join(lines)
    print(report)
    with open('grade_report.txt', 'w') as f:
        f.write(report + "\n")

def run_benchmark(fasta_file):
    """
    Runs the complete benchmark process: parsing, running TranscriptDesigner, validating, and generating reports.
    """
    proteome_name = os.path.splitext(os.path.basename(fasta_file))[0]
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(fasta_file))))

    # Detect GitHub repo URL and any student-added files
    student_files = find_student_files(repo_root)

    # Benchmark the proteome
    parsing_start = time.time()
    successful_results, error_results = benchmark_proteome(fasta_file)
    parsing_time = time.time() - parsing_start

    # Analyze and log exceptions to error_summary.txt
    analyze_errors(error_results)

    # Validate the successful transcripts
    validation_start = time.time()
    validation_failures, checker_passes = validate_transcripts(successful_results)
    execution_time = time.time() - validation_start

    # Sequences that crashed the designer are implicit failures for every check.
    # Fold their count into the translation (ORF) denominator by using total_genes
    # as the denominator throughout; checker_passes already excludes crash victims.
    total_genes = len(successful_results) + len(error_results)

    # Write validation failures for instructor spot-checking
    write_validation_report(validation_failures)

    # Print and save the grade report (includes code links)
    generate_grade_report(total_genes, checker_passes, parsing_time, execution_time, proteome_name, student_files)

if __name__ == "__main__":
    fasta_file = "tests/benchmarking/uniprotkb_proteome_UP000054015_2024_09_24.fasta"
    run_benchmark(fasta_file)
