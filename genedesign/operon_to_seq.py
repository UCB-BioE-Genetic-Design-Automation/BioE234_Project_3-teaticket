from models.operon import Operon
from transcript_to_seq import transcript_to_seq

def operon_to_seq(operon: Operon) -> str:
    """
    Converts a Construct object into its full DNA sequence by concatenating the promoter,
    the sequences of the mRNAs, and the terminator.
    
    Parameters:
        operon (Operon): The construct object containing mRNAs, promoter, and terminator.
    
    Returns:
        str: The full DNA sequence of the construct.
    """
    # Start with the promoter
    out = [operon.promoter]
    
    # Append the sequence of each mRNA by calling transcript_to_seq on each one
    out.extend(transcript_to_seq(mrna) for mrna in operon.transcripts)
    
    # Finally, append the terminator
    out.append(operon.terminator)
    
    # Join the list into a single string and return
    return ''.join(out)
