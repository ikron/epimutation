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

nodeFile <- "/scratch/project_2000350/genomics/methylation/methimpute/control/nodelist_DMR_H3K9_ex_cent_CHH.fn"
edgeFile <- "/scratch/project_2000350/genomics/methylation/methimpute/control/edgelist_DMR.fn"

#nodeFile <- "/scratch/project_2000350/genomics/methylation/nodelist_matA.fn"
#edgeFile <- "/scratch/project_2000350/genomics/methylation/edgelist_matA.fn"

#nodeFile <- "/scratch/project_2000350/genomics/methylation/nodelist_mata.fn"
#edgeFile <- "/scratch/project_2000350/genomics/methylation/edgelist_mata.fn"

output.DMR.H3K9_ex_cent.CHH <- buildPedigreeRegions(nodelist = nodeFile, edgelist = edgeFile, cytosine = "CHH", posteriorMaxFilter = 0.99)

save(output.DMR.H3K9_ex_cent.CHH, file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/DMR_H3K9_ex_cent/MA_DMR_H3K9_ex_cent_CHH.RData")
