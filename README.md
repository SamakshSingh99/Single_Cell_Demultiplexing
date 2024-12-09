# 🧬 Single Cell Hashing Data Demultiplexing With HTO Tags
 
This repository contains a comprehensive R script for analyzing single-cell hashing data, enabling the separation of cells based on Hashtag Oligonucleotide (HTO) tags. It leverages the power of Seurat, Harmony, and other essential R packages to demultiplex single-cell data, detect doublets, and perform downstream analyses like dimensionality reduction and visualization.

## Data
The data for taken from a publically available literature "Putative regulators for the continuum of erythroid differentiation revealed by single-cell transcriptome of human BM and UCB cells" (https://doi.org/10.1073/pnas.1915085117)

## Steps Performed

* __Quality Control:__ Filter out cells with zero HTO tag counts.
* __Demultiplexing:__ Classify cells as singlets, doublets, or negatives using HTO tags.
* __Normalization:__ Apply CLR normalization for HTO data and log normalization for gene expression data.
* __Dimensionality Reduction:__ Perform PCA and UMAP for visualizing clusters.
* __Batch Correction:__ Use Harmony for batch effect correction in merged datasets.
* __Sample-Specific Analysis:__ Extract and analyze cells specific to individual HTO tags.
* __Visualization:__ Generate informative plots for exploratory data analysis.

