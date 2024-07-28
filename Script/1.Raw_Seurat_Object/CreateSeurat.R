###################################################
# Generating Seurat Objects for UCB1 and UCB2UCB3 #
###################################################

#######################
# 1. Loading Packages #
#######################

library(Seurat)
library(tidyverse)

setwd("./parent_directory")

##################################
# 2. Read UCB1 and UCB2UCB3 data #
##################################

UCB1 <- Read10X(data.dir = './zipped_file_UCB1')

UCB2UCB3 <- Read10X(data.dir = './zipped_file_UCB2UCB3')

############################
# 3. Create Seurat Objects #
############################

UCB1_Seurat <- CreateSeuratObject(counts = UCB1, project = 'UCB1')

UCB2UCB3_Seurat <- CreateSeuratObject(counts = UCB2UCB3, project = "UCB2UCB3")



# Saving Seurat Files

save(UCB1_Seurat, file = "./Raw_Format/UCB1_Seurat.Rdata")
save(UCB2UCB3_Seurat, file = "./Raw_Format/UCB2UCB3_Seurat.Rdata")
