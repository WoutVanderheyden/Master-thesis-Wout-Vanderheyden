from Bio import PDB

def filter_pdb(input_pdb, output_pdb, start_res, end_res, chain_id='A'):
    """
    Extracts residues from start_res to end_res for the given chain_id.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    
    io = PDB.PDBIO()
    class SelectResidues(PDB.Select):
        def accept_residue(self, residue):
            return start_res <= residue.id[1] <= end_res and residue.parent.id == chain_id
    
    io.set_structure(structure)
    io.save(output_pdb, select=SelectResidues())
    print(f"Filtered PDB saved as {output_pdb} (Residues {start_res} to {end_res}, Chain {chain_id})")


#### EXAMPLE ON HOW TO USE THIS SCRIPT ####

input_pdb = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/4_YAM_PCA_FULL/2_pdb_files_full/fold_b_finland_33_2010_epi_isl_86993_model_0.pdb"
output_pdb = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/4_YAM_PCA_FULL/12_changed_pdb_files/fold_b_finland_33_2010_epi_isl_86993_model_0_filtered.pdb"
start_res = 16
end_res = 524
chain_id = "A"

filter_pdb(input_pdb, output_pdb, start_res, end_res, chain_id)
