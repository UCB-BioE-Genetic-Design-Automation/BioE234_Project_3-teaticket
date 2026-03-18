from genedesign.seq_utils.reverse_complement import reverse_complement

class ForbiddenSequenceChecker:
    def __init__(self):
        self.forbidden = []

    def initiate(self):
        # Populate forbidden sequences
        self.forbidden = [
            "AAAAAAAA",  # poly(A)
            "TTTTTTTT",  # poly(T)
            "CCCCCCCC",  # poly(C)
            "GGGGGGGG",  # poly(G)
            "ATATATAT",  # poly(AT)
            "CAATTG",    # MfeI
            "GAATTC",    # EcoRI
            "GGATCC",    # BamHI
            "AGATCT",    # BglII
            "ACTAGT",    # SpeI
            "TCTAGA",    # XbaI
            "GGTCTC",    # BsaI
            "CGTCTC",    # BsmBI
            "CACCTGC",   # AarI
            "CTGCAG",    # PstI
            "CTCGAG",    # XhoI
            "GCGGCCGC",  # NotI
            "AAGCTT",    # HindIII
        ]

    def run(self, dnaseq):
        # Use the reverse_complement function from seq_utils
        rc = reverse_complement(dnaseq)
        combined = (dnaseq + "x" + rc).upper()

        for site in self.forbidden:
            if site in combined:
                return False, site

        return True, None

def main():
    checker = ForbiddenSequenceChecker()
    checker.initiate()
    result = checker.run("GGGGGGGGG")  # returns False due to poly(G)
    print(result)

if __name__ == "__main__":
    main()
