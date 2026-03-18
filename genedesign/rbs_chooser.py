from genedesign.models.rbs_option import RBSOption

class RBSChooser:
    """
    A simple RBS selection algorithm that chooses an RBS from a list of options, excluding any RBS in the ignore set.
    """

    def __init__(self):
        self.rbsOptions = []

    def initiate(self) -> None:
        """
        Populates the RBS options list with predefined options.
        """
        # Add RBS options based on their sequences and properties
        opt1 = RBSOption(
            utr="aaagaggagaaatactag",
            cds="atggcttcctccgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtgaaggtcgtccgtacgaaggtacccagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagtacggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaagacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtaccaacttcccgtccgacggtccggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaactgaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggacatcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgcttaa",
            gene_name="BBa_b0034",
            first_six_aas="MASSED"
        )
        opt2 = RBSOption(
            utr="tcacacaggaaagtactag",
            cds="atgactcaacgtatcgcatatgtaactggtggtatgggtggtatcggtactgcaatttgccagcgtctggcgaaagacggtttccgtgttgttgcgggctgcggtccgaactccccgcgtcgtgaaaagtggctggaacaacagaaagccctgggcttcgacttcattgcctccgagggtaatgtagctgactgggattccaccaagactgccttcgataaagttaaatctgaagtgggcgaagtagatgtactgatcaacaacgccggtattactcgtgatgtcgtattccgcaaaatgacccgtgcagactgggatgcagttatcgacaccaacctgacgtctctgttcaacgttaccaaacaggttattgatggtatggctgaccgtggctggggccgcatcgtgaacatctctagcgttaacggccaaaaaggccaatttggtcagacgaattacagcacggctaaagcaggcctgcacggtttcaccatggcactggcgcaggaagtggcgaccaaaggtgttaccgttaataccgtttctccaggttacatcgccaccgatatggttaaggctatccgccaagatgttctggacaagatcgtggctaccattccggttaaacgcctgggcctgccggaagaaattgcgtccatctgtgcgtggctgagctccgaagagtctggtttttccaccggtgcggatttctctctgaacggtggtctgcacatgggttga",
            gene_name="BBa_b0032",
            first_six_aas="MTQRIA"
        )
        opt3 = RBSOption(
            utr="CCATACCCGTTTTTTTGGGCTAACAGGAGGAATTAAcc",
            cds="atgGacacAattaacatcgctaagaacgacttctctgacatcgaactggctgctatcccgttcaacactctggctgaccattacggtgagcgtttagctcgcgaacagttggcccttgagcatgagtcttacgagatgggtgaagcacgcttccgcaagatgtttgagcgtcaacttaaagctggtgaggttgcggataacgctgccgccaagcctctcatcactaccctactccctaagatgattgcacgcatcaacgactggtttgaggaagtgaaagctaagcgcggcaagcgcccgacagccttccagttcctgcaagaaatcaagccggaagccgtagcgtacatcaccattaagaccactctggcttgcctaaccagtgctgacaatacaaccgttcaggctgtagcaagcgcaatcggtcgggccattgaggacgaggctcgcttcggtcgtatccgtgaccttgaagctaagcacttcaagaaaaacgttgaggaacaactcaacaagcgcgtagggcacgtctacaagaaagcatttatgcaagttgtcgaggctgacatgctctctaagggtctactcggtggcgaggcgtggtcttcgtggcataaggaagactctattcatgtaggagtacgctgcatcgagatgctcattgagtcaaccggaatggttagcttacaccgccaaaatgctggcgtagtaggtcaagactctgagactatcgaactcgcacctgaatacgctgaggctatcgcaacccgtgcaggtgcgctggctggcatctctccgatgttccaaccttgcgtagttcctcctaagccgtggactggcattactggtggtggctattgggctaacggtcgtcgtcctctggcgctggtgcgtactcacagtaagaaagcactgatgcgctacgaagacgtttacatgcctgaggtgtacaaagcgattaacattgcgcaaaacaccgcatggaaaatcaacaagaaagtcctagcggtcgccaacgtaatcaccaagtggaagcattgtccggtcgaggacatccctgcgattgagcgtgaagaactcccgatgaaaccggaagacatcgacatgaatcctgaggctctcaccgcgtggaaacgtgctgccgctgctgtgtaccgcaaggacaaggctcgcaagtctcgccgtatcagccttgagttcatgcttgagcaagccaataagtttgctaaccataaggccatctggttcccttacaacatggactggcgcggtcgtgtttacgctgtgtcaatgttcaacccgcaaggtaacgatatgaccaaaggactgcttacgctggcgaaaggtaaaccaatcggtaaggaaggttactactggctgaaaatccacggtgcaaactgtgcgggtgtcgataaggttccgttccctgagcgcatcaagttcattgaggaaaaccacgagaacatcatggcttgcgctaagtctccactggagaacacttggtgggctgagcaagattctccgttctgcttccttgcgttctgctttgagtacgctggggtacagcaccacggcctgagctataactgctcccttccgctggcgtttgacgggtcttgctctggcatccagcacttctccgcgatgctccgagatgaggtaggtggtcgcgcggttaacttgcttcctagtgaaaccgttcaggacatctacgggattgttgctaagaaagtcaacgagattctacaagcagacgcaatcaatgggaccgataacgaagtagttaccgtgaccgatgagaacactggtgaaatctctgagaaagtcaagctgggcactaaggcactggctggtcaatggctggcttacggtgttactcgcagtgtgactaagcgttcagtcatgacgctggcttacgggtccaaagagttcggcttccgtcaacaagtgctggaagataccattcagccagctattgattccggcaagggtctgatgttcactcagccgaatcaggctgctggatacatggctaagctgatttgggaatctgtgagcgtgacggtggtagctgcggttgaagcaatgaactggcttaagtctgctgctaagctgctggctgctgaggtcaaagataagaagactggagagattcttcgcaagcgttgcgctgtgcattgggtaactcctgatggtttccctgtgtggcaggaatacaagaagcctattcagacgcgcttgaacctgatgttcctcggtcagttccgcttacagcctaccattaacaccaacaaagatagcgagattgatgcacacaaacaggagtctggtatcgctcctaactttgtacacagccaagacggtagccaccttcgtaagactgtagtgtgggcacacgagaagtacggaatcgaatcttttgcactgattcacgactccttcggtaccattccggctgacgctgcgaacctgttcaaagcagtgcgcgaaactatggttgacacatatgagtcttgtgatgtactggctgatttctacgaccagttcgctgaccagttgcacgagtctcaattggacaaaatgccagcacttccggctaaaggtaacttgaacctccgtgacatcttagagtcggacttcgcgttcgcAtaa",
            gene_name="Pbad_rbs",
            first_six_aas="MDTINI"
        )
        self.rbsOptions.extend([opt1, opt2, opt3])

    def run(self, cds: str, ignores: set) -> RBSOption:
        """
        Selects an RBS that is not in the ignore set.
        
        Parameters:
            cds (str): The coding sequence.
            ignores (set): A set of RBS options to ignore.
        
        Returns:
            RBSOption: The selected RBS option.
        """
        for rbsopt in self.rbsOptions:
            if rbsopt not in ignores:
                return rbsopt
        raise Exception("No valid RBS option available.")

if __name__ == "__main__":
    # Example usage of RBSChooser
    cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATT..."

    # Initialize the chooser
    chooser = RBSChooser()
    chooser.initiate()

    # Choose RBS with no ignores
    ignores = set()
    selected1 = chooser.run(cds, ignores)
    
    # Add the first selection to the ignore list
    ignores.add(selected1)
    
    # Choose another RBS option after ignoring the first
    selected2 = chooser.run(cds, ignores)

    # Print the selected RBS options
    print("Selected1:", selected1)
    print("Selected2:", selected2)
