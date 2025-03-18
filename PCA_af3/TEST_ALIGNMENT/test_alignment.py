#!/usr/bin/env python3
"""
Script to align chain A of a PDB structure to chain A of a reference PDB structure
by calculating an optimal rotation matrix and translation vector.
Only uses chain A CA atoms for alignment and RMSD calculation.
Only writes chain A to the output PDB.
Also creates a copy of the reference PDB with only chain A.

Usage:
    python align_pdb.py input.pdb reference.pdb output.pdb
"""

import sys
import numpy as np
import os


def parse_pdb_coordinates(pdb_file, chain_id='A'):
    """
    Extract CA atom coordinates from a PDB file for a specific chain.
    Returns a numpy array where each row contains the x, y, z coordinates of a CA atom.
    
    :param pdb_file: Path to the PDB file
    :param chain_id: Chain ID to extract (default 'A')
    :return: numpy array of coordinates
    """
    coordinates = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and " CA " in line and line[21] == chain_id:
                # Extract x, y, z coordinates (columns 31-38, 39-46, 47-54)
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coordinates.append([x, y, z])
    
    if not coordinates:
        print(f"Warning: No CA atoms found for chain {chain_id} in {pdb_file}")
        sys.exit(1)
        
    return np.array(coordinates)


def extract_chain_to_file(pdb_file, output_file, chain_id='A'):
    """
    Extract a specific chain from a PDB file and write it to a new file.
    
    :param pdb_file: Path to the input PDB file
    :param output_file: Path to the output PDB file
    :param chain_id: Chain ID to extract (default 'A')
    """
    with open(pdb_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Keep HEADER, TITLE, etc. lines and only chain A ATOM/HETATM lines
            if not (line.startswith("ATOM") or line.startswith("HETATM")) or line[21] == chain_id:
                outfile.write(line)
    
    print(f"Chain {chain_id} from {pdb_file} written to {output_file}")


def write_aligned_pdb(input_pdb, output_pdb, rot_matrix, input_com, ref_com, chain_id='A'):
    """
    Apply rotation matrix and translation to chain A atoms in the input PDB
    and write the transformed coordinates to the output PDB.
    
    :param input_pdb: Path to the input PDB file
    :param output_pdb: Path to the output PDB file
    :param rot_matrix: 3x3 rotation matrix
    :param input_com: Center of mass of input structure
    :param ref_com: Center of mass of reference structure
    :param chain_id: Chain ID to process (default 'A')
    """
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[21] == chain_id:
                # Extract coordinates
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                
                # Apply rotation and translation
                coords = np.array([x, y, z])
                
                # First center the coordinates
                centered_coords = coords - input_com
                
                # Then rotate
                rotated_coords = rot_matrix @ centered_coords
                
                # Finally translate to reference center of mass
                new_coords = rotated_coords + ref_com
                
                # Write transformed coordinates back to PDB format
                new_line = f"{line[:30]}{new_coords[0]:8.3f}{new_coords[1]:8.3f}{new_coords[2]:8.3f}{line[54:]}"
                outfile.write(new_line)
            elif not (line.startswith("ATOM") or line.startswith("HETATM")):
                # Keep header lines, etc.
                outfile.write(line)

################################################################################
#### STARTING FROM HERE IS JUST THE SAME SCRIPT I USED IN THE OTHER ANALYSIS ###
################################################################################

def normalize(array):
    """Normalizes an 1d array."""
    norm = np.linalg.norm(array)
    if norm == 0:
        return array
    return array / norm


def calc_center_of_mass(coordinate_matrix):
    """Calculates the center of mass of a system with coordinates of each
    residue stored in rows of a 2D numpy array.

    :param coordinate_matrix: matrix whose rows are coordinates of each residue
    :return: vector of the center of mass
    """
    return coordinate_matrix.sum(axis=0) / coordinate_matrix.shape[0]


def shift_center_of_mass(coordinate_matrix):
    """Shifts the center of mass of a system to the origin.

    :param coordinate_matrix: matrix whose rows are coordinates of each residue
    :return: a coordinate_matrix with the center of mass shifted to the origin
    """
    center_of_mass = calc_center_of_mass(coordinate_matrix)
    
    return coordinate_matrix - center_of_mass


def calc_corr_matrix(coord_matrix_0, coord_matrix_1):
    """Calculates the correlation matrix from two coordinate matrices. They
    should have the center of mass aligned.

    :param coord_matrix_0, coord_matrix_1: coordinate matrices whose rows are
        coordinates of each residue
    :return: the 3x3 correlation matrix
    """
    return coord_matrix_0.T @ coord_matrix_1


def calc_f_matrix(R):
    """Calculates a 4x4 F matrix from the correlation matrix R whose
    eigenvectors will correspod to the optimal rotational quaternion.

    :param R: 3x3 correlation matrix calculated by calc_corr_matrix
    :return f_matrix: explained above
    """
    f_matrix = np.array([
        [
            R[0, 0] + R[1, 1] + R[2, 2],
            R[1, 2] - R[2, 1],
            R[2, 0] - R[0, 2],
            R[0, 1] - R[1, 0],
        ],
        [
            R[1, 2] - R[2, 1],
            R[0, 0] - R[1, 1] - R[2, 2],
            R[0, 1] + R[1, 0],
            R[0, 2] + R[2, 0],
        ],
        [
            R[2, 0] - R[0, 2],
            R[0, 1] + R[1, 0],
            R[1, 1] - R[0, 0] - R[2, 2],
            R[1, 2] + R[2, 1],
        ],
        [
            R[0, 1] - R[1, 0],
            R[0, 2] + R[2, 0],
            R[1, 2] + R[2, 1],
            R[2, 2] - R[1, 1] - R[0, 0],
        ],
    ])
    return f_matrix


def calc_rot_matrix(f_matrix):
    """Calculates the optimal rotational matrix from the eigenvector of the F
    matrix corresponding to the biggest eigenvalue.

    :param f_matrix: matrix generated by calc_f_matrix
    :return rot_matrix: the optimal rotational matrix
    """
    e_val, e_vec = np.linalg.eigh(f_matrix)
    q = e_vec[:, -1]
    rot_matrix = np.array([
        [
            (q[0] ** 2 + q[1] ** 2 - q[2] ** 2 - q[3] ** 2),
            2 * (q[1] * q[2] - q[0] * q[3]),
            2 * (q[1] * q[3] + q[0] * q[2]),
        ],
        [
            2 * (q[1] * q[2] + q[0] * q[3]),
            (q[0] ** 2 - q[1] ** 2 + q[2] ** 2 - q[3] ** 2),
            2 * (q[2] * q[3] - q[0] * q[1]),
        ],
        [
            2 * (q[1] * q[3] - q[0] * q[2]),
            2 * (q[2] * q[3] + q[0] * q[1]),
            (q[0] ** 2 - q[1] ** 2 - q[2] ** 2 + q[3] ** 2),
        ],
    ])
    return rot_matrix


def rotate(coordinate_matrix, rot_matrix):
    """Rotates the system described by the coordinate_matrix.

    :param coordinate_matrix: matrix whose rows are coordinates of each residue
    :param rot_matrix: a 3x3 rotation matrix
    :return: rotated coordinate_matrix
    """
    rotated_coordinate_matrix = rot_matrix @ coordinate_matrix.T
    return rotated_coordinate_matrix.T


def calc_coordinate_difference(coord_matrix_0, coord_matrix_1):
    """Calculates the difference between coord_matrix_0 and coord_matrix_1.

    :param coord_matrix_0, coord_matrix_1: coordinate matrices whose rows are
        coordinates of each residue
    :return: flattend 1D coordinates difference vector
    """
    return (coord_matrix_0 - coord_matrix_1).flatten()


def calc_rmsd(coordinate_difference):
    """Calculates the rmsd from the coordinate difference vector of the aligned
    structures.

    :param coordinate_difference: 1D coordinates difference vector
    """
    rmsd = np.linalg.norm(coordinate_difference) / np.sqrt(
        (np.size(coordinate_difference) / 3)
    )
    return rmsd


def main():
    if len(sys.argv) != 4:
        print("Usage: python align_pdb.py input.pdb reference.pdb output.pdb")
        sys.exit(1)
    
    input_pdb = sys.argv[1]
    reference_pdb = sys.argv[2]
    output_pdb = sys.argv[3]
    chain_id = 'A' 

    # Generate filenames for chain A extracts
    input_name = os.path.splitext(os.path.basename(input_pdb))[0]
    ref_name = os.path.splitext(os.path.basename(reference_pdb))[0]
    ref_chain_a_pdb = f"{ref_name}_chainA.pdb"
    
    # Check if files exist
    if not os.path.exists(input_pdb):
        print(f"Error: Input PDB file '{input_pdb}' not found.")
        sys.exit(1)
    if not os.path.exists(reference_pdb):
        print(f"Error: Reference PDB file '{reference_pdb}' not found.")
        sys.exit(1)
    
    # Extract chain A from reference PDB
    extract_chain_to_file(reference_pdb, ref_chain_a_pdb, chain_id)
    
    print(f"Aligning chain {chain_id} of {input_pdb} to chain {chain_id} of {reference_pdb}...")
    
    # Extract chain A CA coordinates
    input_coords = parse_pdb_coordinates(input_pdb, chain_id)
    ref_coords = parse_pdb_coordinates(ref_chain_a_pdb, chain_id)
    
    # Check if the number of CA atoms match
    if input_coords.shape[0] != ref_coords.shape[0]:
        print(f"Warning: Number of CA atoms differs between input ({input_coords.shape[0]}) "
              f"and reference ({ref_coords.shape[0]})")
        print("Continuing with alignment, but RMSD might not be meaningful.")
    
    # Calculate center of mass for both structures
    input_center_of_mass = calc_center_of_mass(input_coords)
    ref_center_of_mass = calc_center_of_mass(ref_coords)
    
    # Center both structures at origin
    centered_input_coords = shift_center_of_mass(input_coords)
    centered_ref_coords = shift_center_of_mass(ref_coords)
    
    # Calculate correlation matrix
    R = calc_corr_matrix(centered_input_coords, centered_ref_coords)
    
    # Calculate F matrix
    F = calc_f_matrix(R)
    
    # Calculate optimal rotation matrix
    rot_matrix = calc_rot_matrix(F)
    print(rot_matrix)
    print(ref_center_of_mass)
    # Align the CA coordinates with the optimal rotation and translation
    aligned_ca_coords = rotate(centered_input_coords, rot_matrix) + ref_center_of_mass
    
    # Calculate RMSD between aligned and reference CA coordinates
    coord_diff = calc_coordinate_difference(aligned_ca_coords, ref_coords)
    rmsd = calc_rmsd(coord_diff)
    
    print(f"RMSD between aligned and reference structures (chain {chain_id} CA atoms only): {rmsd:.4f} Å")
    print(f"Number of CA atoms used for alignment: {input_coords.shape[0]}")
    
    # Write output PDB with rotated and translated coordinates (chain A only)
    write_aligned_pdb(input_pdb, output_pdb, rot_matrix, input_center_of_mass, ref_center_of_mass, chain_id)
    
    print(f"Successfully aligned chain {chain_id} PDB saved to {output_pdb}")
    print(f"Reference chain {chain_id} PDB saved to {ref_chain_a_pdb}")


if __name__ == "__main__":
    main()