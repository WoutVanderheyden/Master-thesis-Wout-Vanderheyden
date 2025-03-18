import os
import shutil
from Bio import PDB

# Define paths
source_dir = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/2_pdb_files"
dest_dir = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/2_pdb_files_full"

# Ensure destination directory exists
os.makedirs(dest_dir, exist_ok=True)

# Set up PDB parser
parser = PDB.PDBParser(QUIET=True)

# Process each PDB file in the source directory
for filename in os.listdir(source_dir):
    if filename.endswith('.pdb'):
        file_path = os.path.join(source_dir, filename)
        
        try:
            # Parse the PDB file
            structure = parser.get_structure('structure', file_path)
            
            # Check if chain A exists
            if 'A' in structure[0]:
                # Count residues in chain A
                chain_a = structure[0]['A']
                residue_count = len(list(chain_a.get_residues()))
                
                print(f"File: {filename}, Chain A residue count: {residue_count}")
                
                # Copy files with more than 450 residues
                if residue_count > 450:
                    dest_path = os.path.join(dest_dir, filename)
                    shutil.copy2(file_path, dest_path)
                    print(f"  → Copied to {dest_dir}")
            else:
                print(f"File: {filename}, Chain A not found")
                
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")

print("\nProcessing complete. Files with >450 residues in chain A were copied.")