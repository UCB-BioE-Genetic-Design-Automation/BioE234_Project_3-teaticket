def reverse_complement(dna_sequence: str) -> str:
    """
    Returns the reverse complement of a DNA sequence.

    Parameters:
        dna_sequence (str): The DNA sequence to reverse complement.

    Returns:
        str: The reverse complement of the DNA sequence.
    """
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'
    }
    return ''.join(complement[base] for base in reversed(dna_sequence))

def main():
    # Example usage of reverse_complement
    dna_sequences = [
        "ATGCGACGTTAA",  # Example 1
        "ATGTTTCCC",     # Example 2
        "ATGTTTTGA"      # Example 3
    ]
    
    for seq in dna_sequences:
        rev_comp = reverse_complement(seq)
        print(f"DNA sequence: {seq} -> Reverse complement: {rev_comp}")

if __name__ == "__main__":
    main()
