# Step-1_Load_and_QC.py

import scanpy as sc
import pandas as pd
import os

# Set the path to the raw count matrix file
data_path = "./Dental Project/GSE227731_DPSC_PDLSC_scRNAseq_raw_count_matrix.csv.gz"

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(data_path, index_col=0)

# Transpose the matrix so that cells are rows and genes are columns
adata = sc.AnnData(df.T)

# Ensure gene names are unique
adata.var_names_make_unique()

# Show basic information about the dataset
print("Dataset shape (cells x genes):", adata.shape)
print("Number of genes:", adata.n_vars)
print("Number of cells:", adata.n_obs)

# Calculate quality control metrics
adata.obs['n_counts'] = adata.X.sum(axis=1)           # Total counts per cell
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)      # Number of genes expressed per cell

# Plot violin plots for quality metrics
sc.pl.violin(adata, ['n_counts', 'n_genes'], jitter=0.4, multi_panel=True, show=True)

# Filter out low-quality cells and genes
sc.pp.filter_cells(adata, min_genes=200)      # Keep cells with at least 200 genes expressed
sc.pp.filter_genes(adata, min_cells=3)        # Keep genes expressed in at least 3 cells

# Create output folder if it doesn't exist
os.makedirs("G:/Dental Project/Processed", exist_ok=True)

# Save the AnnData object for later use
adata.write("G:/Dental Project/Processed/adata_dental_raw_filtered.h5ad")

print("AnnData file successfully saved.")
