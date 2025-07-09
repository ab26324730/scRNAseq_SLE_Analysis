# scRNA-seq Data Analysis for SLE
Dataset : GSE174188

This repository contains R scripts for single-cell RNA-seq data processing and analysis, including:

# Data processing
- Quality control (filtering low-quality cells and doublets)
- Normalization and PCA
- Clustering and annotation (using SingleR with Blueprint/ENCODE)
- Visualization (UMAP, heatmaps, density plots, volcano plots)
- GSEA (Hallmark and C5 gene sets)

## Environment

- R version 4.5.0
- Seurat, SingleR, Nebulosa, msigdbr, and other relevant packages.

## How to run

1. Clone this repo
2. Install required packages
3. Please decompress the feature, barcode, and matrix files first.
4. Repalce the data files direction which you just decompressed.
5. Execute `analysis.R` in your R session.
