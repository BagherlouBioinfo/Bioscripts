import pandas as pd
from scipy.stats import mannwhitneyu

# Load noise data
df = pd.read_excel("Processed/gene_expression_noise_by_group.xlsx")

# Pivot the data to compare DPSC vs PDLSC for each gene
pivot_df = df.pivot(index="Gene", columns="Group", values="CV").dropna()

# Perform statistical test
results = []
for gene, row in pivot_df.iterrows():
    try:
        stat, pval = mannwhitneyu([row["DPSC"]], [row["PDLSC"]], alternative="two-sided")
    except Exception:
        pval = 1.0
    results.append({
        "Gene": gene,
        "CV_DPSC": row["DPSC"],
        "CV_PDLSC": row["PDLSC"],
        "CV_Diff": row["DPSC"] - row["PDLSC"],
        "p_value": pval
    })

# Convert to DataFrame
df_compare = pd.DataFrame(results)

# Adjust p-values using Benjamini-Hochberg
df_compare["p_adj"] = pd.Series(df_compare["p_value"]).rank(method='first') / len(df_compare) * 0.05
df_compare = df_compare.sort_values("p_adj")

# Save results
output_path = "Processed/gene_noise_comparison.xlsx"
df_compare.to_excel(output_path, index=False)
print(f"📊 Noise comparison results saved to: {output_path}")
