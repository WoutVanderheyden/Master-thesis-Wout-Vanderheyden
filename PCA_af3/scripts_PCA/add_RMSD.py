import csv

# File paths
rmsd_file = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/1_H1N1_PCA_FULL/11_RMSD_calculated_cut/rmsd_results.csv"
final_file = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/1_H1N1_PCA_FULL/9_combined_file_cut/updated.csv"
output_file = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/1_H1N1_PCA_FULL/9_combined_file_cut/final.csv"

# Load RMSD values into a dictionary, stripping .npy from filenames
rmsd_dict = {}
with open(rmsd_file, "r") as rmsd_csv:
    reader = csv.reader(rmsd_csv)
    header = next(reader)  # Skip header
    
    for row in reader:
        if len(row) >= 3:  # Ensure normalized RMSD column exists
            filename = row[0].replace(".npy", "")  # Remove .npy
            rmsd_dict[filename] = row[2]  # Store normalized RMSD

# Read final.csv and insert normalized RMSD as column 4
updated_rows = []
with open(final_file, "r") as final_csv:
    reader = csv.reader(final_csv)
    header = next(reader)

    # Insert "Normalized RMSD" at column index 3 (zero-based, so it's the 4th column)
    header.insert(3, "Normalized RMSD")
    updated_rows.append(header)

    for row in reader:
        filename = row[0]
        norm_rmsd = rmsd_dict.get(filename, "N/A")  # Default to "N/A" if not found

        # Insert the RMSD value at column index 2 (before the old column 3)
        row.insert(2, norm_rmsd)
        updated_rows.append(row)

# Write updated data to a new file
with open(output_file, "w", newline="") as output_csv:
    writer = csv.writer(output_csv)
    writer.writerows(updated_rows)

print(f"Updated file with Normalized RMSD as column 4 saved to {output_file}")
