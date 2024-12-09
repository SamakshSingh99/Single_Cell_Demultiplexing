####################################
# Filtering Raw Seurat UCB1 Object #
####################################

#######################
# 1. Loading Packages #
#######################

library(Seurat)
library(ggplot2)

setwd("./parent_directory")

# Reading Raw UCB1 Seurat Data
load("./Raw_Format/UCB2UCB3_Seurat.Rdata")

#####################
# 2. Filtering Data #
#####################

# Calculating Mitochondrial Percentage
UCB2UCB3_Seurat$mito_percent <- PercentageFeatureSet(UCB2UCB3_Seurat, pattern = '^MT-') 

# Setting Library size of cells within 2 SDs around mean value

rna_assay <- UCB2UCB3_Seurat@assays$RNA
total_mrna <- Matrix::colSums(UCB2UCB3_Seurat@assays$RNA$counts)

# Adding Total_mRNA column in the meta.data of UCB1_Seurat
UCB2UCB3_Seurat <- AddMetaData(UCB2UCB3_Seurat, metadata = total_mrna, col.name = "Total_mRNAs")

# Setting 2SD threshold for filtering
upper_bound <- 10^(mean(log10(UCB2UCB3_Seurat@meta.data$Total_mRNAs)) +2*sd(log10(UCB2UCB3_Seurat@meta.data$Total_mRNAs)))
lower_bound <- 10^(mean(log10(UCB2UCB3_Seurat@meta.data$Total_mRNAs)) -2*sd(log10(UCB2UCB3_Seurat@meta.data$Total_mRNAs)))

## Filtering based on the threshold mentioned in the reference 

UCB2UCB3_Seurat <- subset(UCB2UCB3_Seurat, subset = mito_percent <= 5 &
                Total_mRNAs >= lower_bound &
                Total_mRNAs <= upper_bound)
UCB2UCB3_Seurat

#############################
# Normalization and Scaling #
#############################

# Log Normalization
UCB2UCB3_Seurat <- NormalizeData(UCB2UCB3_Seurat, 
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)

UCB2UCB3_Seurat <- ScaleData(UCB2UCB3_Seurat)

UCB2UCB3_Seurat <- FindVariableFeatures(UCB2UCB3_Seurat,
                                    mean.function = ExpMean,
                                    dispersion.function = LogVMR)

save(UCB2UCB3_Seurat, file = "./Raw_Format/Filtered_UCB2UCB3_Seurat.Rdata")
