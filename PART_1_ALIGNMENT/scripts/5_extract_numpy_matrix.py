import numpy as np
import os

# Dir containing CSV files
coordinate_directory = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/5_aligned_coordinate_files_cut'

# Output dir for .npy files
output_directory = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/6_numpy_matrices'

# Create the output directory (if it doesn't exist)
os.makedirs(output_directory, exist_ok=True)


for filename in os.listdir(coordinate_directory):
    if filename.endswith('.csv'): 
        file_path = os.path.join(coordinate_directory, filename)
        
        # Load the CSV data into a NumPy array, skip header + exclude first column
        coordinates = np.genfromtxt(file_path, delimiter=',', skip_header=1, usecols=(1, 2, 3))
        
        print(f"Loaded matrix for {filename} with shape {coordinates.shape}")
        
        # Generate output file path 
        output_file = os.path.join(output_directory, os.path.basename(filename).replace('.csv', '.npy'))
        print(output_file)
        
        # Save the matrix as .npy file
        np.save(output_file, coordinates)

        print(f"Saved NumPy array for {filename} to {output_file}")

