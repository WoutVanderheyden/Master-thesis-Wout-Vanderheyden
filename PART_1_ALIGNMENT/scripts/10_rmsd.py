import numpy as np
import os
import csv
from motion_overlap import calc_coordinate_difference, calc_rmsd

# Define directories
matrix_dir = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/7_rotated_coordinates"
output_file = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/11_RMSD_calculated/rmsd_results.csv"

# Reference matrix
ref_matrix_name = "A_Brisbane_59_2007_0.npy"
ref_matrix_path = os.path.join(matrix_dir, ref_matrix_name)
ref_matrix = np.load(ref_matrix_path)

# Get all .npy files in the directory
matrix_files = [f for f in os.listdir(matrix_dir) if f.endswith(".npy") and f != ref_matrix_name]

# Store results in a list
rmsd_values = []

for matrix_file in matrix_files:
    matrix_path = os.path.join(matrix_dir, matrix_file)
    matrix = np.load(matrix_path)
    
    # Compute coordinate difference and RMSD
    coord_diff = calc_coordinate_difference(ref_matrix, matrix)
    print(coord_diff)
    rmsd_value = calc_rmsd(coord_diff)
    print(rmsd_value)
    
    # Append results to the list
    rmsd_values.append([matrix_file, rmsd_value])

# Normalize RMSD values using min-max normalization
if rmsd_values:
    min_rmsd = min(r[1] for r in rmsd_values)
    max_rmsd = max(r[1] for r in rmsd_values)
    
    if max_rmsd > min_rmsd:  # Avoid division by zero
        for r in rmsd_values:
            r.append((r[1] - min_rmsd) / (max_rmsd - min_rmsd))
    else:
        for r in rmsd_values:
            r.append(0.0)  # If all values are the same, set normalized RMSD to 0

# Sort results by RMSD (2nd column)
rmsd_values.sort(key=lambda x: x[1])

# Write sorted results to CSV
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Compared Matrix", "RMSD", "Normalized RMSD"])
    writer.writerows(rmsd_values)

print(f"RMSD calculations completed. Results saved to {output_file}")
