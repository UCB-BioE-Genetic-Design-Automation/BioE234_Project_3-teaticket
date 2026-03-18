from enum import Enum

class Host(Enum):
    """
    The Host is the organism into which DNAs will be inserted or removed.
    """
    Ecoli = "Escherichia coli"      # A strain of Escherichia coli
    Scerevisiae = "Saccharomyces cerevisiae"  # A strain of Saccharomyces cerevisiae
