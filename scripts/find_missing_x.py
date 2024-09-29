from Bio import SeqIO
from Bio.Seq import Seq
import argparse

def find_missing_aa(protein_file, nucleotide_file):
    # Load the protein and nucleotide sequences
    protein_record = SeqIO.read(protein_file, "fasta")
    nucleotide_record = SeqIO.read(nucleotide_file, "fasta")

    protein_sequence = str(protein_record.seq)
    nucleotide_sequence = str(nucleotide_record.seq)

    # Find the position of 'X'
    position_x = protein_sequence.find("X")

    if position_x != -1:
        # Calculate the corresponding nucleotide position (3 nucleotides per AA)
        start_codon_index = position_x * 3
        codon = nucleotide_sequence[start_codon_index:start_codon_index + 3]
        
        print(f"Codon for 'X' at position {position_x}: {codon}")

        # Translate the codon
        if len(codon) == 3:  # Ensure the codon is complete
            amino_acid = Seq(codon).translate()
            print(f"Corresponding Amino Acid: {amino_acid}")
        else:
            print("Incomplete codon. Unable to translate.")
    else:
        print("No 'X' found in the protein sequence.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find codon for 'X' in protein sequence")
    parser.add_argument("protein_file", help="Path to the input protein FASTA file")
    parser.add_argument("nucleotide_file", help="Path to the input nucleotide FASTA file")
    
    args = parser.parse_args()
    
    find_missing_aa(args.protein_file, args.nucleotide_file)
