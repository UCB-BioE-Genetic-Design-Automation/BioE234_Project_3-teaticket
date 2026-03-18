from dataclasses import dataclass
from typing import List
from .host import Host

@dataclass(frozen=True)
class Composition:
    """
    Describes the specification for a genetically engineered organism derived
    from a cell line specified by the Host augmented with a list of proteins.
    """
    host: Host
    promoter: str
    proteins: List[str]
    terminator: str
