
# UCB_BioE134_GeneDesign

This repository provides a suite of tools for designing genetic constructs, including operons, transcripts, and ribosome binding sites (RBS) for synthetic biology applications.

## Project Structure

The repository contains the following main directories and files:

```
UCB_BioE134_GeneDesign/
│
├── genedesign/
│   ├── operon_to_seq.py
│   ├── operon_designer.py
│   ├── rbs_chooser.py
│   ├── transcript_designer.py
│   ├── transcript_to_seq.py
│   ├── checkers/
│   │   ├── codon_checker.py
│   │   ├── forbidden_sequence_checker.py
│   │   ├── hairpin_checker.py
│   │   └── internal_promoter_checker.py
│   ├── data/
│   │   └── codon_usage.txt
│   ├── models/
│   │   ├── composition.py
│   │   ├── host.py
│   │   ├── operon.py
│   │   ├── rbs_option.py
│   │   └── transcript.py
│   └── seq_utils/
│       ├── translate.py
│       ├── calc_edit_distance.py
│       ├── hairpin_counter.py
│       └── reverse_complement.py
│
├── tests/
│   ├── benchmarking/
│   │   ├── proteome_benchmarker.py
│   │   └── uniprotkb_proteome_UP000054015_2024_09_24.fasta
│   └── unit/
│       ├── checkers/
│       │   ├── test_codon_checker.py
│       │   ├── test_forbidden_sequence_checker.py
│       │   ├── test_internal_promoter_checker.py
│       ├── designer/
│       │   ├── test_operon_designer.py
│       │   └── test_transcript_designer.py
│       └── seq_utils/
│           └── test_hairpin_counter.py
│
├── README.md
├── requirements.txt
├── summary_report.txt
├── validation_failures.tsv
└── error_summary.txt
```

### Key Components

- **genedesign/**: This directory contains the core functionality for designing genetic constructs, including operons, transcripts, and RBS sequences.
  - `operon_designer.py`: Constructs a multi-gene operon sequence by arranging genes, promoters, and terminators based on a given composition. It allows for the design of complex genetic constructs.
  - `transcript_designer.py`: Designs individual transcripts by integrating a ribosome binding site (RBS), coding sequence (CDS), and other elements to ensure proper translation of the gene.
  - `rbs_chooser.py`: Selects optimal ribosome binding site (RBS) sequences to control translation initiation, optimizing gene expression based on the design.
  - `operon_to_seq.py`: Converts operon models into DNA sequences by combining genetic elements into a single continuous sequence ready for synthesis.
  - `transcript_to_seq.py`: Converts designed transcript objects into DNA sequences, generating the final nucleotide sequence of the transcript.

- **checkers/**: Contains sequence validation modules that ensure the designed constructs are free from errors and potential regulatory issues.
  - `codon_checker.py`: Validates the codon usage in a sequence, checking codon diversity, rare codon count, and calculating the Codon Adaptation Index (CAI) to ensure the sequence is optimized for the host organism.
  - `forbidden_sequence_checker.py`: Detects forbidden sequences that may interfere with proper gene function, including restriction sites or undesired motifs.
  - `hairpin_checker.py`: Detects secondary structures like hairpins in the sequence, which can cause issues in gene expression.
  - `internal_promoter_checker.py`: Detects internal promoter sequences that could lead to unintended gene expression within the construct.

- **models/**: Contains data models used across the project to represent genetic components and structures.
  - `composition.py`: Represents a genetic composition, including its parts (e.g., promoter, genes).
  - `host.py`: Defines different host organisms (e.g., _E. coli_, _S. cerevisiae_).
  - `operon.py`: Represents a genetic operon, which consists of multiple transcripts, a promoter, and a terminator.
  - `rbs_option.py`: Describes RBS sequences as modular components to control translation initiation.
  - `transcript.py`: Represents a transcript, including its RBS and coding sequence (CDS).

- **seq_utils/**: Utility scripts for handling DNA and protein sequence operations.
  - `translate.py`: Handles the translation of DNA sequences into corresponding protein sequences.
  - `calc_edit_distance.py`: Computes the edit distance between two sequences, useful for comparing genetic variants.
  - `hairpin_counter.py`: Detects potential hairpin structures in nucleotide sequences that could disrupt transcription or translation.
  - `reverse_complement.py`: Computes the reverse complement of a DNA sequence, often needed in cloning or analysis workflows.


## How to Install and Run

This project has minimal dependencies, which can be installed using `pip` and a virtual environment (venv) to avoid conflicts with other Python projects.

### Requirements

- **Python 3.x** is required. Make sure it's installed on your system.
- **Dependencies**: Listed in the `requirements.txt` file.

### Steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/UCB-BioE-Anderson-Lab/UCB_BioE134_GeneDesign.git
   cd UCB_BioE134_GeneDesign
   ```

2. Create a virtual environment:

   #### On Windows:
   ```bash
   python -m venv venv
   venv\Scripts\activate
   ```

   #### On macOS/Linux:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```

3. Install dependencies using `pip`:
   ```bash
   pip install -r requirements.txt
   ```

4. Run the scripts directly. For example, to design an operon:
   ```bash
   python genedesign/operon_designer.py
   ```

5. **Running Tests**: You can run the tests using `pytest` to ensure everything works correctly:
   ```bash
   pytest
   ```

If you run into errors finding paths, such as "ModuleNotFoundError: No module named 'genedesign" try putting this into the command line:
   ```bash
   export PYTHONPATH=$(pwd)
   ```

6. To deactivate the virtual environment when finished:
   ```bash
   deactivate
   ```

### Usage

To design your genetic constructs:
- Instantiate a `Composition` object to reflect your desired genetic elements (promoters, genes, terminators, etc.).
- Run OperatorDesigner (or TranscriptDesigner, or other lower functions) to compute the resulting Operon
- These function scripts have main methods with runnable examples
   ```python
    from genedesign.operon_designer import OperonDesigner
    from genedesign.construct_to_seq import operon_to_seq
    from genedesign.models.composition import Composition
    from genedesign.models.host import Host

    # Example promoter and terminator sequences
    Pbad_Promoter = "TTATGACAACTTGACGGCTACATCATTCACTTTTTCTTCACAACCGGCACG..."
    TrrnB_Terminator = "TGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACT..."

    # Protein sequences (translated into codon-optimized DNA sequences in the final output)
    PaIPDS = "MALLSSSLSSQIPTGSHPLTHTQCIPHFSTTINAGISAGKPRSF..."
    crtE = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERD..."

    # Create the Composition object with promoter, proteins, and terminator
    proteins = [PaIPDS, crtE]
    comp = Composition(host=Host.Ecoli, promoter=Pbad_Promoter, proteins=proteins, terminator=TrrnB_Terminator)

    # Initialize and run the algorithm to design the operon
    designer = OperonDesigner()
    designer.initiate()
    construct_result = designer.run(comp)

    # Convert the operon design into a complete DNA sequence
    output_seq = operon_to_seq(construct_result)

    # Print the resulting DNA sequence
    print(output_seq)
   ```
### Expected Output
The output of running the scripts will be a complete DNA sequence, representing either an entire operon or individual mRNA transcripts. These outputs consist of sequences for the promoter, ribosome binding sites (RBS), coding sequences for proteins, and terminators.

Multiple operons may be incorporated into larger genetic structures. For example, a typical plasmid might contain distinct operons for expressing a protein of interest, antibiotic resistance genes, and an origin of replication, all within a single circular DNA sequence.