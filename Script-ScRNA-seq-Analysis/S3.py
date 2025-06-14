# Step_3_DentalProject_Normalize_PCA_UMAP.py

import scanpy as sc

# Load the labeled AnnData object
adata = sc.read("./Dental Project/Processed/adata_dental_labeled.h5ad")

# Step 1: Normalize total counts per cell to 10,000 and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Step 2: Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Keep only highly variable genes
adata = adata[:, adata.var.highly_variable]

# Step 3: Scale the data
sc.pp.scale(adata, max_value=10)

# Step 4: PCA
sc.tl.pca(adata, svd_solver='arpack')

# Step 5: Compute the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Step 6: UMAP projection
sc.tl.umap(adata)

# Step 7: Leiden clustering
sc.tl.leiden(adata, resolution=0.5)

# Step 8: Plot UMAP colored by clusters and CellSource
sc.pl.umap(adata, color=['leiden', 'CellSource'], save="_leiden_CellSource.png")

# Save the processed AnnData object
adata.write("G:/Dental Project/Processed/adata_dental_normalized_clustered.h5ad")

print("Normalization and clustering complete. UMAP plot saved.")
