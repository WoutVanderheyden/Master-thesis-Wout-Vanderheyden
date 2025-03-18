import numpy as np
import os

aligned_directory = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/7_rotated_coordinates'

csv_output_directory = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/8_rotated_csv'
os.makedirs(csv_output_directory, exist_ok=True)


for filename in os.listdir(aligned_directory):
    if filename.endswith('.npy'):  
        file_path = os.path.join(aligned_directory, filename)

        # Load the aligned coordinate matrix from the .npy file
        aligned_matrix = np.load(file_path)

        # Define the output file path for the CSV file
        csv_output_file = os.path.join(csv_output_directory, f"{os.path.splitext(filename)[0]}.csv")
        print(csv_output_file)

        # Save the matrix as a CSV file
        np.savetxt(csv_output_file, aligned_matrix, delimiter=',', header='x,y,z', comments='', fmt='%8.3f')
        print(f"CSV file for {filename} has been saved to {csv_output_file}")
