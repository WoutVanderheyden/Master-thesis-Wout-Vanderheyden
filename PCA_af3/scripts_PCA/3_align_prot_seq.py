import os
from Bio.Align.Applications import MafftCommandline
import subprocess

def run_mafft_alignment(input_fasta, output_aln):
    """
    Run MAFFT alignment on sequences from a FASTA file
    """
    try:
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"Input file {input_fasta} not found")
            
        output_dir = os.path.dirname(output_aln)
        os.makedirs(output_dir, exist_ok=True)
        
        mafft_cline = MafftCommandline(input=input_fasta, auto=True)
        
        cmd = str(mafft_cline)
        print(f"Running MAFFT command: {cmd}")
        
        # Run MAFFT 
        result = subprocess.run(
            cmd.split(),
            capture_output=True,
            text=True,
            check=True
        )
        
        # Save alignment
        with open(output_aln, "w") as aln_file:
            aln_file.write(result.stdout)
            
        print(f"Alignment completed and saved to: {output_aln}")
        
        # Print warnings (if there are any)
        if result.stderr:
            print("MAFFT Warnings/Messages:")
            print(result.stderr)
            
    except FileNotFoundError as e:
        print(f"Error: {e}")
        raise
    except subprocess.CalledProcessError as e:
        print(f"MAFFT execution failed with error: {e}")
        if e.stderr:
            print(f"MAFFT error output: {e.stderr}")
        raise
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise

input_fasta = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/4_YAM_PCA_FULL/4_align_coordinates/YAM.fasta"
output_aln = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/4_YAM_PCA_FULL/4_align_coordinates/YAM_aligned.fasta"

if __name__ == "__main__":
    run_mafft_alignment(input_fasta, output_aln)
