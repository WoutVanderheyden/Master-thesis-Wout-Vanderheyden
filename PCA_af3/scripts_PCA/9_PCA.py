import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import os
from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import MinMaxScaler

# Input combined CSV file 
combined_file_path = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/9_combined_file/combined_coordinates_by_structure.csv'

# Load data
data = pd.read_csv(combined_file_path)

# Extract filenames (first column) and numerical data
filenames = data.iloc[:, 0]  # First column contains names
# year = data.iloc[:, 1]
# rmsd = data.iloc[:, 2].astype(float)  # Ensure RMSD is numeric
numeric_data = data.iloc[:, 1:]  # Exclude non-numeric columns for PCA

# Perform PCA
pca = PCA()
pca_result = pca.fit_transform(numeric_data)

# Create DF for PCA results
pca_df = pd.DataFrame(pca_result, columns=[f"PC{i+1}" for i in range(pca_result.shape[1])])
pca_df.insert(0, 'filename', filenames)
# pca_df.insert(1, 'year', year)
# pca_df.insert(2, 'rmsd', rmsd)  # Include RMSD in PCA dataframe

# Output directory
output_directory = '/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/10_PCA_results'
os.makedirs(output_directory, exist_ok=True)

# Save PCA results
pca_output_file = f"{output_directory}/output_PCA.csv"
pca_df.to_csv(pca_output_file, index=False)
print(f"PCA results saved to {pca_output_file}")

####  Scree plot (only first 15 PCs) #######
num_pcs_to_show = 15
plt.figure(figsize=(10, 6))
plt.plot(range(1, num_pcs_to_show + 1), pca.explained_variance_ratio_[:num_pcs_to_show], marker='o', linestyle='--')
plt.title(f'Scree Plot (First {num_pcs_to_show} PCs)')
plt.xlabel('Principal Component')
plt.ylabel('Explained Variance Ratio')
plt.xticks(range(1, num_pcs_to_show + 1))
scree_plot_path = f"{output_directory}/scree_plot_first_{num_pcs_to_show}_PCs.png"
plt.savefig(scree_plot_path)
plt.close()
print(f"Scree plot saved to {scree_plot_path}")

#### 2D Scatter plot of the first two PCs #######
plt.figure(figsize=(10, 6))
plt.scatter(pca_df['PC1'], pca_df['PC2'], alpha=0.7)
for i, txt in enumerate(pca_df['filename']):
    plt.annotate(txt, (pca_df['PC1'][i], pca_df['PC2'][i]), fontsize=8, alpha=0.7)
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA Scatter Plot (PC1 vs PC2)')
scatter_plot_path = f"{output_directory}/pca_2D_scatter.png"
plt.savefig(scatter_plot_path)
plt.close()
print(f"2D PCA scatter plot saved to {scatter_plot_path}")


###############################################################################################
################### TAKE A LOOK AT THE LOADINGS OF PCA ########################################
###############################################################################################

# Load dataset (excluding filename column)
df = pd.read_csv(combined_file_path)
data = df.iloc[:, 3:].values 

# Perform PCA
pca = PCA(n_components=6)  # Adjust as needed
pca.fit(data)

# Get the loadings (contributions of each coordinate)
loadings = pca.components_

# Convert to a DataFrame for easier interpretation
loadings_df = pd.DataFrame(loadings.T, columns=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6"])
print(loadings_df)

# Number of residues
num_residues = loadings_df.shape[0] // 3  # Each residue has 3 coordinates (X, Y, Z)

# Initialize arrays for residue contributions
residue_loadings_PC1 = np.zeros(num_residues)
residue_loadings_PC2 = np.zeros(num_residues)
residue_loadings_PC3 = np.zeros(num_residues)

for i in range(num_residues):
    x_index = i * 3
    y_index = i * 3 + 1
    z_index = i * 3 + 2
    
    # Compute Euclidean norms   
    norm_x_PC1 = np.abs(loadings_df.loc[x_index, "PC1"])
    norm_y_PC1 = np.abs(loadings_df.loc[y_index, "PC1"])
    norm_z_PC1 = np.abs(loadings_df.loc[z_index, "PC1"])
    
    norm_x_PC2 = np.abs(loadings_df.loc[x_index, "PC2"])
    norm_y_PC2 = np.abs(loadings_df.loc[y_index, "PC2"])
    norm_z_PC2 = np.abs(loadings_df.loc[z_index, "PC2"])
    
    norm_x_PC3 = np.abs(loadings_df.loc[x_index, "PC3"])
    norm_y_PC3 = np.abs(loadings_df.loc[y_index, "PC3"])
    norm_z_PC3 = np.abs(loadings_df.loc[z_index, "PC3"])
    
    # Individual PC contributions
    residue_loadings_PC1[i] = np.sqrt(norm_x_PC1**2 + norm_y_PC1**2 + norm_z_PC1**2)
    residue_loadings_PC2[i] = np.sqrt(norm_x_PC2**2 + norm_y_PC2**2 + norm_z_PC2**2)
    residue_loadings_PC3[i] = np.sqrt(norm_x_PC3**2 + norm_y_PC3**2 + norm_z_PC3**2)

# Function to save residue contributions with normalized values
def save_residue_contributions(filename, sorted_indices, contributions):
    # Normalize contributions to range [0-1]
    normalized_contributions = contributions / np.max(contributions)
    
    output_df = pd.DataFrame({
        "residue_nr": sorted_indices + 1,  # Residue numbers (1-based index)
        "total_contribution": contributions[sorted_indices],
        "normalized_contribution": normalized_contributions[sorted_indices]
    })
    output_df.to_csv(filename, index=False)

# Sort and save results
base_path = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/10_PCA_results/"

sorted_indices_PC1 = np.argsort(residue_loadings_PC1)[::-1]
sorted_indices_PC2 = np.argsort(residue_loadings_PC2)[::-1]
sorted_indices_PC3 = np.argsort(residue_loadings_PC3)[::-1]

save_residue_contributions(f"{base_path}residue_contributions_PC1.csv", sorted_indices_PC1, residue_loadings_PC1)
save_residue_contributions(f"{base_path}residue_contributions_PC2.csv", sorted_indices_PC2, residue_loadings_PC2)
save_residue_contributions(f"{base_path}residue_contributions_PC3.csv", sorted_indices_PC3, residue_loadings_PC3)

# Optional: Save a copy of the raw loadings for reference
loadings_df.to_csv(f"{base_path}raw_pc_loadings.csv", index=True)