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

######################################
## OBJECTIVE 1. SEPARATING THE DATA ##
## FROM THE CELL HASHING FILE TO    ##
## PERFORM SCRNASeq ANALYSIS.       ##
######################################

UCB23_File <- Read10X("/Users/samaksh/Repositories/scRNA_HTO_Demultiplexing/zipped_file_UCB2UCB3")

###############################
# Extracting Info on HTO Tags #
###############################

HTO_labs <- as.data.frame(UCB23_File$`Antibody Capture`)
rownames(HTO_labs)

 # "B0251 anti-human Hashtag1" "B0252 anti-human Hashtag2"

#########################################
# Performing QC For Downstream Analysis #
#########################################

pucb.htos <- HTO_labs[colSums(HTO_labs)>0] # Barcode ID with 0 in both samples removed

# Gene Expression Matrix

pucb.umis <- UCB23_File$`Gene Expression`

# QC Perfromed on The Surface Protein, Performing Intersection

joint.bcs <- intersect(colnames(pucb.umis), colnames(pucb.htos))

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



