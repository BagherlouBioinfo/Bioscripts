# Step_2:Label_DPSC_PDLSC.py

import scanpy as sc

# Load the filtered AnnData object
adata = sc.read("./Dental Project/Processed/adata_dental_raw_filtered.h5ad")

# Assign cell type based on whether the cell name includes "DPSC"
adata.obs['CellSource'] = adata.obs_names.to_series().apply(lambda x: 'DPSC' if 'DPSC' in x else 'PDLSC')

# Show summary of cell counts per group
print("Cell counts per group:")
print(adata.obs['CellSource'].value_counts())

# Save updated AnnData object
adata.write("G:/Dental Project/Processed/adata_dental_labeled.h5ad")

print("Labeled AnnData saved with 'CellSource' column.")
