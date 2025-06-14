# S6.py

import pandas as pd
import os

# Define file paths
input_file = "Processed/gene_noise_comparison.xlsx"
output_file = "Processed/genes_lower_noise_in_PDLSC.xlsx"

# Create output directory if it doesn't exist
os.makedirs("Processed", exist_ok=True)

# Load comparison data
df = pd.read_excel(input_file)

# Filter genes with significantly lower noise in PDLSC (i.e., CV_Diff > 0)
filtered_df = df[(df["CV_Diff"] > 0) & (df["p_adj"] < 0.05)].copy()
filtered_df.sort_values("CV_Diff", ascending=False, inplace=True)

# Save to Excel
filtered_df.to_excel(output_file, index=False)

print(f"✅ {filtered_df.shape[0]} genes with lower noise in PDLSC found.")
print(f"📁 Results saved to: {output_file}")
