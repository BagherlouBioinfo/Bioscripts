#Step_4
import scanpy as sc
import pandas as pd
import numpy as np

# Load processed data
adata = sc.read("Processed/adata_dental_normalized_clustered.h5ad")

# Ensure group label exists
if "CellSource" not in adata.obs.columns:
    raise ValueError("Missing 'CellSource' column in metadata.")

# Initialize result dictionary
results = []

# For each group: DPSC and PDLSC
for group in adata.obs["CellSource"].unique():
    print(f"Calculating for group: {group}")
    # Subset data
    group_adata = adata[adata.obs["CellSource"] == group]
    X = group_adata.X.toarray() if hasattr(group_adata.X, "toarray") else group_adata.X

    # Mean and std for each gene across cells
    gene_means = np.mean(X, axis=0)
    gene_stds = np.std(X, axis=0)
    gene_cvs = gene_stds / (gene_means + 1e-6)  # Add small value to avoid division by zero

    # Save results
    for gene, mean, std, cv in zip(adata.var_names, gene_means, gene_stds, gene_cvs):
        results.append({
            "Gene": gene,
            "Group": group,
            "Mean": mean,
            "Std": std,
            "CV": cv
        })

# Convert to DataFrame
df_noise = pd.DataFrame(results)

# Save to Excel
output_path = "Processed/gene_expression_noise_by_group.xlsx"
df_noise.to_excel(output_path, index=False)
print(f"Noise metrics saved to: {output_path}")

