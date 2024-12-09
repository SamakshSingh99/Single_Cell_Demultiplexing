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

####################
# Working with HTO #
####################

# HTO labels enable sample multiplexing, doublet detection, and demultiplexing in single-cell hashing experiments.

HTO <- as.data.frame(UBC$`Antibody Capture`) 
rownames(HTO) # Checking the presence of HTO labels

# [1] "B0251 anti-human Hashtag1" "B0252 anti-human Hashtag2"



