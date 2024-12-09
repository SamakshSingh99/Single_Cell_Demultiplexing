#####################################################
# This script performes separation of a single-cell #
# hashing data on the basis of HTO-tags.            #
#####################################################

#######################
# 1. Loading Packages #
#######################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringr)
library(harmony)

##############################################
# 2. Loading the single-cell hashing dataset #
##############################################
setwd("/Users/samaksh/Repositories/Single_Cell_Demultiplexing/Data")

UBC <- Read10X("./")
UBC

#########################
# 3. Observing HTO tags #
#########################

# HTO labels enable sample multiplexing, doublet detection, and demultiplexing in single-cell hashing experiments.

HTO <- as.data.frame(UBC$`Antibody Capture`) 
rownames(HTO) # Checking the presence of HTO labels

# [1] "B0251 anti-human Hashtag1" "B0252 anti-human Hashtag2"

######################
# 4. Quality Control #
######################

# Checking the expression matrix to exclude barcode IDs which are 0 in each sample
ubc.hto <- HTO[colSums(HTO)>0]

UBC.ge <- UBC$`Gene Expression` # extracting gene expression from UBC seurat object

# Intersecting gene expression and antibody capture with barcode ID > 0 

joint.filt_barcodes <- intersect(colnames(UBC.ge), colnames(ubc.hto)) # Obtaining intersection information

UBC.ge <- UBC.ge[, joint.filt_barcodes]

ubc.hto <- as.matrix(ubc.hto[, joint.filt_barcodes])

##############################################
# 5. Exploratory data analysis and importing #
# Seurat object for downstream analysis.     #
##############################################

# Creating new seurat object
UBC.tags <- CreateSeuratObject(counts = UBC.ge)

# Normalization
UBC.tags <- NormalizeData(UBC.tags)

# Identifying features with high variability
UBC.tags <- FindVariableFeatures(UBC.tags, selection.method = "mean.var.plot")

# Scaling seurat object
UBC.tags <- ScaleData(UBC.tags, features = VariableFeatures(UBC.tags))
UBC.tags

# An object of class Seurat 
# 32738 features across 18777 samples within 1 assay 
# Active assay: RNA (32738 features, 1749 variable features)
# 3 layers present: counts, data, scale.data

##########################################
# 6. Adding HTO tag matrix for different #
# sample source identification.          #
##########################################

UBC.tags [["HTO"]] <- CreateAssayObject(counts = ubc.hto)
UBC.tags

# An object of class Seurat 
# 32740 features across 18777 samples within 2 assays 
# Active assay: RNA (32738 features, 1749 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: HTO

# Central log ratio (CLR) normalization for HTO labels
UBC.tags <- NormalizeData(UBC.tags,
                          assay = "HTO",
                          normalization.method = "CLR")


head(UBC.tags@meta.data)

###############################
# 7. Demultiplexing with HTOs #
###############################

UBC.tags <- HTODemux(UBC.tags,
                     assay = "HTO",
                     positive.quantile = 0.99)

table(UBC.tags$HTO_classification.global)

# Doublet Negative  Singlet 
# 152    14650     3975

UBC.tags@assays

# $RNA
# Assay (v5) data with 32738 features for 18777 cells
# Top 10 variable features:
#  HBE1, S100A9, CLC, HBD, PF4, PPBP, S100A8, HSPA1A, LYZ, NKG7 
# Layers:
#  counts, data, scale.data 

# $HTO
# Assay data with 2 features for 18777 cells
# First 2 features:
#  B0251 anti-human Hashtag1, B0252 anti-human Hashtag2 

# Confirmation for division
table(UBC.tags$hash.ID)

# Doublet                  Negative                  B0251 anti-human Hashtag1                 B0252 anti-human Hashtag2 
# 152                     14650                      2144                                      1831 

################################
# 8. Standardization Operation #
################################

# Filtering negatives and doublets

UBC.HTO.subset <- subset(UBC.tags,
                         idents = c("Negative", "Doublet"),
                         invert = T)

table(Idents(UBC.HTO.subset))

# B0251 anti-human Hashtag1 B0252 anti-human Hashtag2 
# 2144                      1831 

# Dimensionality reduction using PCA and UMAP

UBC.HTO.subset <- RunPCA(UBC.HTO.subset)

UBC.HTO.subset <- RunUMAP(UBC.HTO.subset,dims=1:10)

DimPlot(UBC.HTO.subset) 

#######################################
# 9. Sample specific cells extraction #
#######################################

UBC2 <- subset(UBC.tags,
               idents = c("B0251 anti-human Hashtag1"))

UBC3 <- subset(UBC.tags,
               idents = c("B0252 anti-human Hashtag2"))

# Merging UBC2 and UBC3 seurat objects

UBC2UBC3 <- merge(UBC2, y= UBC3, add.cell.ids = c('UBC2', 'UBC3'))

# Assigning the corresponding cell ids in the metadata

phe = as.data.frame(str_split(colnames(UBC2UBC3), '_', simplify = TRUE))

UBC2UBC3@meta.data$orig.ident = phe$V1

table(UBC2UBC3@meta.data$orig.ident )

# UBC2 UBC3 
# 2144 1831 

View(UBC2UBC3@meta.data)


# Merged Seurat standardization

###############
# 10. Harmony #
###############

# Running processes upto PCA as per harmony requirements

# Data Normalization
UBC2UBC3 <- NormalizeData(UBC2UBC3,
                          normalization.method = "LogNormalize",
                          scale.factor = 1e4) 

# Variable Feature Identification
UBC2UBC3 <- FindVariableFeatures(UBC2UBC3)

# Data Scaling
UBC2UBC3 <- ScaleData(UBC2UBC3)

# Dimensionality reduction 

UBC2UBC3 <- RunPCA(UBC2UBC3, features = VariableFeatures(object=UBC2UBC3)) # PCA

UBC.seurat <- RunHarmony(UBC2UBC3, "orig.ident") #Harmony Reduction

names(UBC.seurat@reductions)

UBC.seurat <- RunUMAP(UBC.seurat,  dims = 1:15, 
                     reduction = "harmony") 

UMAP <- UMAPPlot(UBC.seurat, group.by = "orig.ident")
DimPlot(UBC.seurat,reduction = "umap",label=F )
