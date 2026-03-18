from dataclasses import dataclass

@dataclass
class Translate:
    """
    Translates a DNA sequence into a protein sequence using the standard genetic code, halting at the first stop codon encountered and throwing an error for invalid codons.

    Attributes:
        codon_table (dict): Maps each DNA codon to its corresponding single-letter amino acid code.
    """
    codon_table: dict = None

    def initiate(self) -> None:
        """
        Initializes the codon table with the genetic code for translating nucleotide triplets into amino acids.
        """
        self.codon_table = {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "TAT": "Y", "TAC": "Y", "TAA": "Stop", "TAG": "Stop",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "TGT": "C", "TGC": "C", "TGA": "Stop", "TGG": "W",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        }

    def run(self, dna_sequence: str) -> str:
        """
        Translates a DNA sequence into a protein sequence using the codon table.

        Parameters:
            dna_sequence (str): The DNA sequence to translate.

        Returns:
            str: The corresponding amino acid sequence.

        Raises:
            ValueError: If the DNA sequence length is not a multiple of three, contains untranslated sequence after a stop codon, or contains invalid codons.
        """
        if len(dna_sequence) % 3 != 0:
            raise ValueError("The DNA sequence length must be a multiple of 3.")

        protein = []
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]
            if codon not in self.codon_table:
                raise ValueError(f"Invalid codon '{codon}' encountered in DNA sequence.")
            amino_acid = self.codon_table[codon]
            if amino_acid == "Stop":
                if i + 3 != len(dna_sequence):
                    raise ValueError("Untranslated sequence after stop codon.")
                break
            protein.append(amino_acid)

        return ''.join(protein)

def main():
    # Example usage
    translator = Translate()
    translator.initiate()
    
    dna_sequences = [
        "ATGCGACGTTAA",  # Example with a stop codon at the end
        "ATGTTTCCC",     # Another example without a stop codon
        "ATGTTTTGA"      # Example with a stop codon in the middle
    ]
    
    for seq in dna_sequences:
        try:
            protein_sequence = translator.run(seq)
            print(f"DNA sequence: {seq} -> Protein sequence: {protein_sequence}")
        except ValueError as e:
            print(f"Error for sequence {seq}: {e}")

if __name__ == "__main__":
    main()
