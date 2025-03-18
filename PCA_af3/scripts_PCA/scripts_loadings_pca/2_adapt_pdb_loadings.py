import pandas as pd

def update_b_factors(pdb_file, csv_file, output_pdb):
    """
    Update the B-factor column in a PDB file based on residue contribution values from a CSV file.
    
    Args:
    - pdb_file (str): Path to the input PDB file.
    - csv_file (str): Path to the CSV file with residue contributions.
    - output_pdb (str): Path to save the modified PDB file.
    """

    # Load residue contribution data
    contributions = pd.read_csv(csv_file)
    contribution_dict = dict(zip(contributions["residue_nr"], contributions["total_contribution"]))

    updated_lines = []

    with open(pdb_file, "r") as pdb:
        for line in pdb:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract residue number (positions 22-26 in PDB format)
                residue_number = int(line[22:26].strip())

                ##############################################################
                #################### IMPORTANT ###############################
                ##############################################################

                # start is coordinate 17 but in the file its different!!!!!!!!!!
                residue_number = residue_number - 15

                # If the residue is in our contribution data, update the B-factor (positions 61-66)
                if residue_number in contribution_dict:
                    contribution_value = f"{contribution_dict[residue_number]:6.2f}"  # Format to fit PDB width
                    new_line = line[:60] + contribution_value.rjust(6) + line[66:]
                    updated_lines.append(new_line)
                else:
                    updated_lines.append(line)
            else:
                updated_lines.append(line)  # Keep non-ATOM lines unchanged

    # Save the modified PDB
    with open(output_pdb, "w") as out_pdb:
        out_pdb.writelines(updated_lines)

    print(f"Updated PDB file saved as: {output_pdb}")


# Example Usage
pdb_file = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/4_YAM_PCA_FULL/12_changed_pdb_files/fold_b_finland_33_2010_epi_isl_86993_model_0_filtered.pdb"  
csv_file = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/4_YAM_PCA_FULL/10_PCA_results_cut/residue_contributions_PC1.csv"
output_pdb = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/4_YAM_PCA_FULL/12_changed_pdb_files/finland_colour.pdb"
update_b_factors(pdb_file, csv_file, output_pdb)
