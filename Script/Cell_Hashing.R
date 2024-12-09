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

