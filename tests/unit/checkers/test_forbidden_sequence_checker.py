import pytest
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker

@pytest.fixture
def checker():
    checker = ForbiddenSequenceChecker()
    checker.initiate()
    return checker

def test_forbidden_sites(checker):
    # Sequences with forbidden sites
    forbidden_seqs = [
        "TTGACAATTgaattcCGAACTAGTATAAT",
        "TTTTTTTTTTGTCGAGAAATTTATAAT",
        "TTGACATTACCGTCTCGAGCGCCTATAAT",
        "TTGACAGCGAACGCTTCAGACTGCAGAT",
        "TTGACTGCAGTTGTAACTTATATAAT",
        "TTGACAAGAAGCGGCCGCTCAATTATAAT",
        "TTGACATTATGACACCTGCTTATTATAAT",
        "TTGACACTCGAGCACAGGCTCTATAAT",
        "TTGACATCGGGGGGGGTTTTACCATGGTCGTTATAAT",
        "TTGACAAAGTCGATTTTTTTTTTCCTTCGATTATAAT",
        "TTGACACGGTCTCATTCACTAGGTTATAAT",
        "TTGACAGTCTAGAGTCTGAACAAGGAGATCTTAAT",
    ]

    print(">> Testing problematic sequences (expected False)")
    for seq in forbidden_seqs:
        result, site = checker.run(seq.upper())
        print(f"result: {result} on {seq}")
        assert result == False

def test_allowed_sites(checker):
    # Sequences without forbidden sites
    allowed_seqs = [
        "AAACTGTAATCCACCACAAGTCAAGCCAT",
        "GCCTCTCTGAGGACGCCGTATGAATTAATA",
        "GTAAACTTTGCGCGGGTTCACTGCGATCC",
        "TTCAGTCTCGTCCAAGGGCACAATCGAAT",
        "ATCCCCCGAAGTTTAGCAGGTCGTGAGGT",
        "TCATGGAGGCTCTCGTTCATCCCGTGGGA",
        "ATCAAGGCTTCGCCTTGATAAAGCACCCCG",
        "TCGGGTGTAGCAGAGAAGACGCCTACTGA",
        "TTGTGCGATCCCTCCACCTCAGCTAAGGT",
        "GCTACCAATATTTAGTTTTTTAGCCTTGC",
        "ACAGACCTCCTACTTAGATTGCCACGCAT",
        "GTTCCGCTGGCGATCCATCGTTGGCGGCCG",
    ]

    print("\n>> Testing random sequences (expected True)")
    for seq in allowed_seqs:
        result, site = checker.run(seq)
        print(f"result: {result} on {seq}")
        assert result == True
