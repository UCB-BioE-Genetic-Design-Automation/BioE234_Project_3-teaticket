from dataclasses import dataclass
from typing import List
from .transcript import Transcript  # Assuming Transcript is defined in transcript.py

@dataclass(frozen=True)
class Operon:
    """
    Encodes a genetic construct described in terms of a single operon
    encoding multiple cocistronic transcripts.
    """
    transcripts: List[Transcript]
    promoter: str
    terminator: str
