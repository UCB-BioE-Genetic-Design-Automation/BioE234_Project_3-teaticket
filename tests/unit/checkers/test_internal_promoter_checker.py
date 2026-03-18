import pytest
from genedesign.checkers.internal_promoter_checker import PromoterChecker

@pytest.fixture
def promoter_checker():
    p = PromoterChecker()
    p.initiate()
    return p

def test_constitutive_promoters(promoter_checker):
    seq1 = [
        "TTGACAATTAATCATCGAACTAGTATAAT",
        "TTGACATCTACTGTCGAGAAATTTATAAT",
        "TTGACATTACTGACTTGAGCGCCTATAAT",
        "TTGACAGCGAACGCTTCAGATGTTATAAT",
        "TTGACATCCCAAATTGTAACTTATATAAT",
        "TTGACAAGAAGGAGTGAATCAATTATAAT",
        "TTGACATTATGTTTTATTATTATTATAAT",
        "TTGACAGTATTCTGCACAGGCTCTATAAT",
        "TTGACATCTTTTACCATGGTCGTTATAAT",
        "TTGACAAAGTCGATTCCTTCGATTATAAT",
        "TTGACACCGCTCATTCACTAGGTTATAAT",
        "TTGACAGGGTGGCTGAACAAGGGTATAAT"
    ]
    
    print("\n>> Testing constitutive promoters (expected False)")
    for seq in seq1:
        result, promoter = promoter_checker.run(seq)
        print(f"Result: {result} on {seq}, Promoter: {promoter}")
        assert result == False  # Expecting constitutive promoters to return False

def test_random_sequences(promoter_checker):
    seq2 = [
        "AAACTGTAATCCACCACAAGTCAAGCCAT",
        "GCCTCTCTGAGACGCCGTATGAATTAATA",
        "GTAAACTTTGCGCGGGTTCACTGCGATCC",
        "TTCAGTCTCGTCCAAGGGCACAATCGAAT",
        "ATCCCCCGAAGTTTAGCAGGTCGTGAGGT",
        "TCATGGAGGCTCTCGTTCATCCCGTGGGA",
        "ATCAAGCTTCGCCTTGATAAAGCACCCCG",
        "TCGGGTGTAGCAGAGAAGACGCCTACTGA",
        "TTGTGCGATCCCTCCACCTCAGCTAAGGT",
        "GCTACCAATATTTAGTTTTTTAGCCTTGC",
        "GTTCCGCTGGGATCCATCGTTGGCGGCCG"
    ]

    print("\n>> Testing random sequences (expected True)")
    for seq in seq2:
        result, promoter = promoter_checker.run(seq)
        print(f"Result: {result} on {seq}, Promoter: {promoter}")
        assert result == True  # Expecting random sequences to return True

def test_j23119_promoters(promoter_checker):
    # List of sequences with expected results
    # Selected from https://parts.igem.org/Promoters/Catalog/Anderson highest (False) and lowest (True) strength promoters
    j23119_sequences = [
        ('TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC',	 False),
        ('TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC',	 False),
        ('TTGACAGCTAGCTCAGTCCTAGGTACTGTGCTAGC',	 False),
        ('TTGACAGCTAGCTCAGTCCTAGGTATTGTGCTAGC',	 False),
        ('TTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGC',	 False),
        ('TTGACGGCTAGCTCAGTCCTAGGTATAGTGCTAGC',	 False),
        ('TTGACGGCTAGCTCAGTCCTAGGTATTGTGCTAGC',	 False),
        ('CTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC',	 False),
        ('TTTACGGCTAGCTCAGTCCTAGGTATAGTGCTAGC',	 False),
        ('TTTACGGCTAGCTCAGCCCTAGGTATTATGCTAGC',	 False),
        ('TTTACAGCTAGCTCAGTCCTAGGGACTGTGCTAGC',	 True),
        ('CTGATAGCTAGCTCAGTCCTAGGGATTATGCTAGC',	 True),
        ('CTGATGGCTAGCTCAGTCCTAGGGATTATGCTAGC',	 True),
        ('CTGATAGCTAGCTCAGTCCTAGGGATTATGCTAGC',	 True),
    ]

    print("\n>> Testing J23119 promoters with expected results")
    for seq, expected in j23119_sequences:
        result, promoter = promoter_checker.run(seq)
        print(f"Sequence: {seq}, Expected: {expected}, Got: {result}, Promoter: {promoter}")
        assert result == expected, f"Test failed for sequence: {seq}. Expected {expected} but got {result}."
