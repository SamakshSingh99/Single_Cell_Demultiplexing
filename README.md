# 🧬 Single Cell Hashing Data Demultiplexing With HTO Tags
 
This repository contains a comprehensive R script for analyzing single-cell hashing data, enabling the separation of cells based on Hashtag Oligonucleotide (HTO) tags. It leverages the power of Seurat, Harmony, and other essential R packages to demultiplex single-cell data, detect doublets, and perform downstream analyses like dimensionality reduction and visualization.

## About

__Single Cell Hashing__
Single-cell hashing is a technique used to label individual cells with unique molecular tags (e.g., DNA-barcoded antibodies or chemical tags) before pooling them for single-cell sequencing. This method allows researchers to multiplex multiple samples or conditions within the same sequencing run. By applying unique barcodes to each sample, single-cell hashing enables the identification of the sample of origin for each individual cell after sequencing. This technique helps reduce experimental costs and minimizes batch effects, making it a valuable tool in high-throughput single-cell analysis.

__Single Cell Demultiplexing__
Single-cell demultiplexing is the computational process used to assign individual cells back to their respective sample or donor of origin after sequencing. By analyzing the barcode sequences (from single-cell hashing) or genotype information, demultiplexing algorithms differentiate and group cells based on their origin, enabling the separation of pooled samples. This process is essential for accurately analyzing and comparing data from different samples or conditions within a single experiment, allowing for more precise biological insights from multiplexed single-cell sequencing.

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

## Outputs

### __Raw Data Contents__
<img width="930" alt="1_UBC_Seurat_HTO" src="https://github.com/user-attachments/assets/6b3e10d9-5e7b-4f30-a51d-786c118845a4">

### __UMAP Dimensional Reduction Plot After Demultiplexing With HTO Tags__
![UMAP_HTOs](https://github.com/user-attachments/assets/571c1437-763d-409b-977f-9cd750efc6b8)

### __UMAP Dimensional Reduction Plot After Harmony Dimentional Reduction__
![UBC_Harmony](https://github.com/user-attachments/assets/36608aed-85eb-4b7c-be0b-bbdcb8628c7e)



