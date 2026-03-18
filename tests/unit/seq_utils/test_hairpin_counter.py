import pytest
from genedesign.seq_utils.hairpin_counter import hairpin_counter

def test_no_hairpin():
    sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAA"
    count, hairpins = hairpin_counter(sequence)
    assert count == 0, "Expected no hairpins, but got a non-zero count."
    assert hairpins is None, "Expected no hairpin string, but got one."

def test_hairpin_detected():
    sequence = "CCCCCTTTCCCCCCAAACCCCCC"  # Example of a weak hairpin
    count, hairpins = hairpin_counter(sequence)
    assert count > 0, "Expected a hairpin, but got zero."
    assert hairpins is not None, "Expected a hairpin string, but got None."

def test_strong_hairpin():
    sequence = "AAAAACCCCCAAAAAAAAGGGGGAAA"  # Example of a strong hairpin
    weak_sequence = "AAAAACCCAAAAAAAAAGGGAAAAAA"  # Example of a weaker hairpin
    
    strong_count, strong_hairpins = hairpin_counter(sequence)
    weak_count, weak_hairpins = hairpin_counter(weak_sequence)
    
    assert strong_count > weak_count, "Expected stronger hairpin to have a higher count than weaker one."
    assert strong_count > 0, "Expected a hairpin in the strong sequence, but got zero."
    assert weak_count > 0, "Expected a hairpin in the weak sequence, but got zero."
    assert strong_hairpins is not None, "Expected a hairpin string in the strong sequence, but got None."
    assert weak_hairpins is not None, "Expected a hairpin string in the weak sequence, but got None."

@pytest.mark.parametrize("sequence", [
    "GGCTAATTTAGCCATTAAGGCTAATAGGCTAA",  # Example: may or may not form a hairpin
    "AAAAACACGAAAAAAAACGTGAAAAAA",       # Potential hairpin
    "AAAACCCCCAAAAAAAAGGGGGAAA",         # Strong hairpin
])
def test_varied_hairpins(sequence):
    count, hairpins = hairpin_counter(sequence)
    assert count >= 0, "Hairpin count should never be negative."
    if count > 0:
        assert hairpins is not None, "Expected a hairpin string, but got None."
    else:
        assert hairpins is None, "Expected no hairpin string, but got one."
