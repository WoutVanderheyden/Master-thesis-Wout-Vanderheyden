import pandas as pd

# Load the second CSV
df2 = pd.read_csv("/data/leuven/358/vsc35887/master_thesis/PCA_af3/1_H1N1_PCA_FULL/9_combined_file_cut/combined_coordinates_by_structure.csv")

# Extract everything starting from "b_" onwards
df2["base_filename"] = df2["filename"].str.extract(r"(b_.*)").fillna("")  # Avoid NaN issues
df2["base_filename"] = df2["base_filename"].str.replace('_', '')

# Extract the final part of the filename (after the last '_')
df2["year"] = df2["filename"].str.split('_').str[-1]

# Insert the new 'year' column as the 2rd column
df2.insert(1, "year", df2.pop("year"))

df2.drop(columns=["base_filename"], inplace=True)

print(df2)

# Save the modified CSV
df2.to_csv("/data/leuven/358/vsc35887/master_thesis/PCA_af3/1_H1N1_PCA_FULL/9_combined_file_cut/updated.csv", index=False)