# Transcriptomic Noise Analysis in Dental Stem Cells

This repository is part of a larger research project aimed at demonstrating bioinformatics and coding proficiency through the analysis of single-cell RNA sequencing data.

The current project uses Python to analyze transcriptional noise between human dental pulp stem cells (DPSC) and periodontal ligament stem cells (PDLSC) using the dataset GSE227731.

## Summary

The goal is to identify genes with significantly lower expression variability in PDLSC compared to DPSC, based on coefficient of variation (CV), and to functionally annotate them using external tools such as Enrichr.

## Contents

The project is organized into seven Python scripts:

| File | Purpose |
|------|---------|
| S1.py | Load raw scRNA-seq data and perform quality control |
| S2.py | Label cells as DPSC or PDLSC |
| S3.py | Normalize data, run PCA and UMAP, and cluster cells |
| S4.py | Calculate gene-wise expression noise (mean, SD, CV) per group |
| S5.py | Compare noise metrics between groups using statistical tests |
| S6.py | Filter genes with lower noise in PDLSC |
| S7.py | Submit filtered gene list to Enrichr for pathway analysis |

## Data

The raw dataset used in this project is available from GEO under accession GSE227731. The raw `.csv.gz` matrix is not included in this repository due to size.

## Dependencies

Python 3.8+  
Required packages are listed in `requirements.txt`. To install:

```bash
pip install -r requirements.txt
```

## Disclaimer

This repository is not a complete publication but a component of a larger pipeline under development. For inquiries, please contact the author.

