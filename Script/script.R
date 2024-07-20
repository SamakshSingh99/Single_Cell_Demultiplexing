###########################
## LOADING THE LIBRARIES ##
###########################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringr)
library(harmony)
library(gridExtra)

#############################
# Setting Primary Directory #
#############################
prim_dir <- list.dirs(path = "/Users/samaksh/Repositories/scRNA_HTO_Demultiplexing/", recursive = T, full.names = F)

UCB23_File <- Read10X("/Users/samaksh/Repositories/scRNA_HTO_Demultiplexing/zipped_file_UCB2UCB3")

########################################
## OBJECTIVE 1. UNDERSTANDING AND     ##
## VISUALIZATION OF RAW UCB2UCB3 DATA ##
########################################
Seurat_UCB23 <- CreateSeuratObject(UCB23_File, project = "Raw_UCB")
View(Seurat_UCB23@meta.data)

Seurat_UCB23 <- CreateSeuratObject(UCB23_File)

# Normalization
Seurat_UCB23 <- NormalizeData(Seurat_UCB23, 
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)

Seurat_UCB23 <- ScaleData(Seurat_UCB23)

Seurat_UCB23 <- FindVariableFeatures(Seurat_UCB23,
                                     mean.function = ExpMean,
                                     dispersion.function = LogVMR)

Seurat_UCB23 <- RunPCA(Seurat_UCB23)

######################################
## OBJECTIVE 2. SEPARATING THE DATA ##
## FROM THE CELL HASHING FILE TO    ##
## PERFORM SCRNASeq ANALYSIS.       ##
######################################

Seurat_UCB23 <- CreateSeuratObject(UCB23_File, project = "Raw_UCB")
View(Seurat_UCB23@meta.data)

###############################
# Extracting Info on HTO Tags #
###############################

HTO_labs <- as.data.frame(UCB23_File$`Antibody Capture`)
rownames(HTO_labs)


 # "B0251 anti-human Hashtag1" "B0252 anti-human Hashtag2"

#########################################
# Performing QC For Downstream Analysis #
#########################################

pucb.htos <- HTO_labs[,colSums(HTO_labs)>0] # Barcode ID with 0 in both samples removed

# Gene Expression Matrix

pucb.umis <- UCB23_File$`Gene Expression`

# QC Perfromed on The Surface Protein, Performing Intersection

joint.bcs <- intersect(colnames(pucb.umis), colnames(pucb.htos))

# Subset the matrices to include only the joint barcodes
pucb.umis <- pucb.umis[, joint.bcs]
pucb.htos <- pucb.htos[, joint.bcs]

UCB_seurat <-  CreateSeuratObject(counts = pucb.umis, project = "UCB")

#Normalizing The Data 

UCB_seurat <- NormalizeData(UCB_seurat)

# Finding High Variables

UCB_seurat <- FindVariableFeatures(UCB_seurat, selection.method = "mean.var.plot")

# Standardization

UCB_seurat <- ScaleData(UCB_seurat, features = VariableFeatures(UCB_seurat))

View(UCB_seurat@meta.data)

UCB_seurat

# An object of class Seurat 
# 32738 features across 20885 samples within 1 assay 
# Active assay: RNA (32738 features, 0 variable features)
# 3 layers present: counts, data, scale.data

###########################
# Adding HTO Label Matrix #
###########################

UCB_seurat[["HTO"]] <- CreateAssayObject(counts = pucb.htos)

UCB_seurat

head(UCB_seurat@meta.data)

# An object of class Seurat 
# 32740 features across 20885 samples within 2 assays 
# Active assay: RNA (32738 features, 0 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: HTO

##############
# Core Steps #
##############

UCB_seurat <- HTODemux(UCB_seurat, assay = "HTO", positive.quantile = 0.99)


