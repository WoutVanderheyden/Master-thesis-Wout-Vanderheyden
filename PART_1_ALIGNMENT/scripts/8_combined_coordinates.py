import numpy as np
import os
import pandas as pd

# Directory containing the rotated CSV files
rotated_csv_directory = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/8_rotated_csv'

# Output directory for combined CSV file
combined_output_directory = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/9_combined_file'
os.makedirs(combined_output_directory, exist_ok=True)

# Lists to store structure data and filenames
combined_data_by_structure = []
filenames_list = []

for filename in os.listdir(rotated_csv_directory):
    if filename.endswith('.csv'):  
        file_path = os.path.join(rotated_csv_directory, filename)

        # Load the CSV data into pandas DF
        df = pd.read_csv(file_path)

        # Extract x,y,z columns
        coordinates = df[['x', 'y', 'z']].values

        # Reshape -> x1, y1, z1, x2, y2, z2,...
        reshaped_coordinates = coordinates.flatten()

        # Process the filename
        processed_filename = filename.replace('.csv', '')

        # Store the processed filename
        filenames_list.append(processed_filename)

        # Append the reshaped coordinates as a new row
        combined_data_by_structure.append(reshaped_coordinates)

# Convert to DataFrame (to handle filenames as first column)
combined_df = pd.DataFrame(combined_data_by_structure)

# Insert filenames as the first column
combined_df.insert(0, 'filename', filenames_list)

# Output
combined_csv_file_by_structure = os.path.join(combined_output_directory, 'combined_coordinates_by_structure.csv')

# Save using DataFrame (to handle mixed string/float format)
combined_df.to_csv(combined_csv_file_by_structure, index=False, float_format='%.3f')
print(combined_df)

print(f"Combined CSV file by structure has been saved to {combined_csv_file_by_structure}")

