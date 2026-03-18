import math
from genedesign.seq_utils.reverse_complement import reverse_complement

class PromoterChecker:
    """
    A class to check for the presence of constitutive sigma70 promoters in a DNA sequence.

    The class uses a Position Weight Matrix (PWM) derived from a hardcoded Position Frequency Matrix (PFM)
    to scan a sequence of DNA and evaluate whether a constitutive promoter is present. The `run` method
    evaluates both the input sequence and its reverse complement.

    Attributes:
        pwm: A 2D list representing the Position Weight Matrix (PWM) used to score sequences.
    """

    def __init__(self):
        """
        Initializes the PromoterChecker by setting pwm to None.
        The PWM will be computed later in the initiate method.
        """
        self.pwm = None

    def initiate(self):
        """
        Initializes the Position Weight Matrix (PWM) based on a hardcoded Position Frequency Matrix (PFM).
        The PFM represents the nucleotide frequencies in 12 known sequences with constitutive promoters.
        This matrix is converted to a PWM, which is used to score other DNA sequences.
        """
        # Hard code a PFM matrix representing the nucleotide frequencies at 29 positions in known promoter sequences.
        pfm = [
            [0, 0, 0, 12, 0, 12, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 12, 0, 12, 12, 0],  # A
            [0, 0, 0, 0, 12, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0],  # C
            [0, 0, 12, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0],  # G
            [12, 12, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 12, 0, 12, 0, 0, 12]   # T
        ]

        # Create an empty PWM matrix with the same dimensions as the PFM.
        ncols = len(pfm[0])
        self.pwm = [[0] * ncols for _ in range(4)]

        # Convert the PFM to a PWM by calculating the weight of each base at each position.
        for x in range(ncols):
            total = sum(pfm[y][x] for y in range(4))  # Total count for each position
            for y in range(4):
                prob_base = 0.25  # Assuming equal probability for each base in a random sequence
                freq = pfm[y][x]
                # Compute the log-odds weight for each nucleotide at each position
                w = (math.log((freq + math.sqrt(total) * prob_base) / (total + math.sqrt(total)) / prob_base)) / math.log(2)
                self.pwm[y][x] = w

    def run(self, seq):
        """
        Checks if the given DNA sequence contains a constitutive sigma70 promoter.

        The method scores the input sequence and its reverse complement using a sliding window approach.
        If a windowed sequence has a score above a certain threshold, it is considered to contain a promoter.

        Parameters:
            seq (str): A DNA sequence to check.

        Returns:
            tuple: (bool, str or None)
                - bool: True if no promoter is found, False if a promoter is found.
                - str: The promoter sequence if found, None otherwise.
        """
        # Convert the sequence to uppercase and compute its reverse complement.
        seq = seq.upper()
        rc = reverse_complement(seq)
        combined = seq + "x" + rc  # Concatenate the original sequence and its reverse complement.

        sliding_frame = 29  # The sliding window size is 29 nucleotides.
        threshold = 9.134   # A threshold score for detecting promoter activity.

        # Slide over the sequence and calculate the score for each window.
        for i in range(len(combined) - sliding_frame + 1):
            score = 0.0
            partseq = combined[i:i + sliding_frame]
            for x, base in enumerate(partseq):
                # Map the base to its corresponding row in the PWM
                y = {'A': 0, 'C': 1, 'G': 2, 'T': 3}.get(base, -1)
                if y != -1:
                    score += self.pwm[y][x]
            # If the score exceeds the threshold, the sequence likely contains a constitutive promoter.
            if score >= threshold:
                return False, partseq  # Promoter found, return the sequence
        return True, None  # No promoter detected in the sequence


if __name__ == "__main__":
    checker = PromoterChecker()
    checker.initiate()

    constitutive = "TTGACAATTAATCATCGAACTAGTATAAT"
    constitutiveBroken = "TTctgAATTAATCATCGAACTAGgcgAAT"

    # Example checks
    result, promoter = checker.run(constitutive)
    print(f"constitutive: {result}, Promoter: {promoter}")  # Expected output: False (contains promoter)

    result, promoter = checker.run(constitutiveBroken)
    print(f"constitutiveBroken: {result}, Promoter: {promoter}")  # Expected output: True (no promoter)

    # Example of testing a list of sequences
    seq3 = [
        "ttgacagctagctcagtcctaggtataatgctagc",
        "ttgacggctagctcagtcctaggtacagtgctagc",
        "ttgacagctagctcagtcctaggtactgtgctagc",
        "ttgacagctagctcagtcctaggtattgtgctagc",
        "tttacagctagctcagtcctaggtattatgctagc",
        "ttgacggctagctcagtcctaggtatagtgctagc",
        "ttgacggctagctcagtcctaggtattgtgctagc",
        "ctgacagctagctcagtcctaggtataatgctagc",
        "tttacggctagctcagtcctaggtatagtgctagc",
        "tttacggctagctcagccctaggtattatgctagc",
        "tttacggctagctcagtcctaggtacaatgctagc",
        "tttacggctagctcagtcctaggtactatgctagc",
        "ttgacagctagctcagtcctagggactatgctagc",
        "tttatagctagctcagcccttggtacaatgctagc",
        "tttatggctagctcagtcctaggtacaatgctagc",
        "ttgacagctagctcagtcctagggattgtgctagc",
        "tttacagctagctcagtcctagggactgtgctagc",
        "ctgatagctagctcagtcctagggattatgctagc",
        "ctgatggctagctcagtcctagggattatgctagc",
        "ctgatagctagctcagtcctagggattatgctagc"
    ]

    print("\n>> Testing J23119 promoters (mixed results)")
    for seq in seq3:
        result, promoter = checker.run(seq)
        print(f"Result: {result}, Promoter: {promoter}")  # Outputs either True (no promoter) or False (contains promoter)
