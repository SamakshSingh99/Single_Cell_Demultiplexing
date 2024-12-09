####################################
# Filtering Raw Seurat UCB1 Object #
####################################

#######################
# 1. Loading Packages #
#######################

library(Seurat)
library(ggplot2)
library(tidyverse)

setwd("./parent_directory")

# Reading Raw UCB1 Seurat Data
load("./Raw_Format/UCB1_Seurat.Rdata")

#####################
# 2. Filtering Data #
#####################
view(UCB1_Seurat@meta.data)

# Calculating Mitochondrial Percentage
UCB1_Seurat[["mito.perc"]] <- PercentageFeatureSet(UCB1_Seurat, pattern = '^MT-') 


# Setting Library size of cells within 2 SDs around mean value

rna_assay <- UCB1_Seurat@assays$RNA
total_mrna <- Matrix::colSums(UCB1_Seurat@assays$RNA$counts)

# Adding Total_mRNA column in the meta.data of UCB1_Seurat
UCB1_Seurat <- AddMetaData(UCB1_Seurat, metadata = total_mrna, col.name = "Total_mRNAs")

# Setting 2SD threshold for filtering
upper_bound <- 10^(mean(log10(UCB1_Seurat@meta.data$Total_mRNAs)) +2*sd(log10(UCB1_Seurat@meta.data$Total_mRNAs)))
lower_bound <- 10^(mean(log10(UCB1_Seurat@meta.data$Total_mRNAs)) -2*sd(log10(UCB1_Seurat@meta.data$Total_mRNAs)))

## Filtering based on the threshold mentioned in the reference 

UCB1_Seurat <- subset(UCB1_Seurat, subset = mito_percent <= 5 &
                Total_mRNAs >= lower_bound &
                Total_mRNAs <= upper_bound)
UCB1_Seurat

#############################
# Normalization and Scaling #
#############################

# Log Normalization
UCB1_Seurat <- NormalizeData(UCB1_Seurat, 
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)

UCB1_Seurat <- ScaleData(UCB1_Seurat)

UCB1_Seurat <- FindVariableFeatures(UCB1_Seurat,
                                    mean.function = ExpMean,
                                    dispersion.function = LogVMR)

save(UCB1_Seurat, file = "./Raw_Format/Filtered_UCB1_Seurat.Rdata")
