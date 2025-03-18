import numpy as np
import os
from motion_overlap import align_coordinates

# FOR NOW JUST CHOOSE THIS 
reference_matrix = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/6_numpy_matrices/A_Brazil_11_1978_4.npy'
coord_matrix_1 = np.load(reference_matrix)

coordinate_directory = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/6_numpy_matrices'
output_directory = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/7_rotated_coordinates'
os.makedirs(output_directory, exist_ok=True)

# Save the reference matrix aswell in output directory
np.save(os.path.join(output_directory, f"aligned_{os.path.basename(reference_matrix)}"), coord_matrix_1)


# Iterate over all .npy files (excluding the reference matrix)
for filename in os.listdir(coordinate_directory):
    if filename.endswith('.npy') and filename != os.path.basename(reference_matrix):
        file_2 = os.path.join(coordinate_directory, filename)
        print(filename)
        # Load the coordinate matrix to be aligned
        coord_matrix_2 = np.load(file_2)

        # Align the matrices using the align_coordinates function
        aligned_matrix_1 = align_coordinates(coord_matrix_2, coord_matrix_1)

        # Define the output file path for the aligned matrix
        output_file = os.path.join(output_directory, f"{filename}")

        # Save aligned matrix 
        np.save(output_file, aligned_matrix_1)
        print(f"Aligned coordinates for {filename} have been saved to {output_file}")
