from dataclasses import dataclass
from typing import List
from .rbs_option import RBSOption  # Assuming RBSOption is defined in rbs_option.py

@dataclass(frozen=True)
class Transcript:
    """
    Encodes a monocistronic mRNA from an RBS and a coding sequence.
    """
    rbs: RBSOption
    peptide: str
    codons: List[str]
