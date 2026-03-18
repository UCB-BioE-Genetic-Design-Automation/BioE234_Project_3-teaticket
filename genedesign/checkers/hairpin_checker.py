from genedesign.seq_utils.hairpin_counter import hairpin_counter

def hairpin_checker(dna):
    """
    Checks for bad hairpin structures in the DNA sequence by splitting it into 50 bp chunks with
    an overlap of 25 bp, and passing each chunk to hairpin_counter. If any chunk has more than
    1 hairpin, it returns False and the problematic hairpin string. Otherwise, it returns True and None.

    Parameters:
        dna (str): The DNA sequence to analyze.

    Returns:
        tuple: (bool, str or None)
            - True and None if no problematic hairpins are found.
            - False and the problematic hairpin string if more than one hairpin is found in any chunk.
    """
    chunk_size = 50  # 50 bp window
    overlap = 25     # Overlap by 25 bp
    min_stem = 3     # Minimum number of bases in the stem
    min_loop = 4     # Minimum number of bases in the loop
    max_loop = 9     # Maximum number of bases in the loop
    
    # Iterate over the sequence in 50 bp chunks with 25 bp overlap
    for i in range(0, len(dna) - chunk_size + 1, overlap):
        chunk = dna[i:i + chunk_size]
        
        # Get the count of hairpins and the hairpin string from hairpin_counter
        hairpin_count, hairpin_string = hairpin_counter(chunk, min_stem, min_loop, max_loop)
        
        # If more than 1 hairpin is found, return False and the problematic hairpin string
        if hairpin_count > 1:
            return False, hairpin_string
    
    # If no problematic hairpin chunk is found, return True and None
    return True, None

# Example usage
if __name__ == "__main__":
    result, hairpin = hairpin_checker("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCAAAAAAAGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    print(result, hairpin)
