##R script to calculate divergence data for different cytosince contexts
# Setting the paths
.libPaths(c("/projappl/project_2000350/rpackages", .libPaths()))

library(tidyverse)
library(data.table)
library(gtools)
library(igraph)
library(BiocParallel)
#library(AlphaBeta)
source("/scratch/project_2000350/genomics/methylation/buildPedigree_Regions.R") #Modified version of AlphaBeta (multicore setup needs to be different)

multicoreParam <- MulticoreParam(workers=4) #Needed to set up multicore parameters


nodeFile.matA <- "/scratch/project_2000350/genomics/methylation/methimpute/control/nodelist_matA_DMR_CHG.fn"
edgeFile.matA <- "/scratch/project_2000350/genomics/methylation/methimpute/control/edgelist_DMR_matA.fn"

nodeFile.mata <- "/scratch/project_2000350/genomics/methylation/methimpute/control/nodelist_mata_DMR_CHG.fn"
edgeFile.mata <- "/scratch/project_2000350/genomics/methylation/methimpute/control/edgelist_DMR_mata.fn"

#Calculate pedigree divergence data for both pedigrees separately
output.matA.DMR.all.CHG <- buildPedigreeRegions(nodelist = nodeFile.matA, edgelist = edgeFile.matA, cytosine = "CHG", posteriorMaxFilter = 0.99)

output.mata.DMR.all.CHG <- buildPedigreeRegions(nodelist = nodeFile.mata, edgelist = edgeFile.mata, cytosine = "CHG", posteriorMaxFilter = 0.99)

save(output.matA.DMR.all.CHG, output.mata.DMR.all.CHG, file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/DMR_all/MA_DMR_all_mat_CHG.RData")
