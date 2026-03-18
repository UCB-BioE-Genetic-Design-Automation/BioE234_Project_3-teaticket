from dataclasses import dataclass

@dataclass(frozen=True)
class RBSOption:
    """
    Encapsulates Ribosome Binding Site (RBS) encoding DNAs as modular components (parts) for synthetic biology,
    representing essential sequence elements to facilitate selection algorithms.

    Attributes:
        utr (str): The 5' untranslated region (5' UTR) sequence
        cds (str): The coding sequence of the source gene
        gene_name (str): The name of the source gene
        first_six_aas (str): The precalculated first six amino acids of the source gene's protein sequence
    """
    utr: str
    cds: str
    gene_name: str
    first_six_aas: str
