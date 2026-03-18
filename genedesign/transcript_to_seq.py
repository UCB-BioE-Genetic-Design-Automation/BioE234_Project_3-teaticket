from models.transcript import Transcript

def transcript_to_seq(transcript: Transcript) -> str:
    """
    Converts a Transcript object into its mRNA sequence by concatenating the RBS and the codons.

    Parameters:
        transcript (Transcript): The transcript object containing RBS, peptide, and codons.

    Returns:
        str: The mRNA sequence with the RBS in lowercase and the CDS in uppercase.
    """
    # Build the mRNA sequence from the codons
    codon_sequence = ''.join(transcript.codons)
    
    # Return the RBS (lowercase) concatenated with the codon sequence (uppercase)
    return transcript.rbs.utr.lower() + codon_sequence.upper()
