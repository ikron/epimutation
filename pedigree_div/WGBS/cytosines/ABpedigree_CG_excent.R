##R script to calculate divergence data for different cytosince contexts
# Setting the paths
.libPaths(c("/projappl/project_2000350/rpackages", .libPaths()))

library(tidyverse)
library(data.table)
library(gtools)
library(igraph)
library(BiocParallel)
#library(AlphaBeta)
source("/scratch/project_2000350/genomics/methylation/ABmod.R") #Modified version of AlphaBeta (multicore setup needs to be different)

multicoreParam <- MulticoreParam(workers=4) #Needed to set up multicore parameters

nodeFile <- "/scratch/project_2000350/genomics/methylation/methimpute/control/nodelist_CG_excent.fn"
edgeFile <- "/scratch/project_2000350/genomics/methylation/methimpute/control/edgelist.fn"

#nodeFile <- "/scratch/project_2000350/genomics/methylation/nodelist_matA.fn"
#edgeFile <- "/scratch/project_2000350/genomics/methylation/edgelist_matA.fn"

#nodeFile <- "/scratch/project_2000350/genomics/methylation/nodelist_mata.fn"
#edgeFile <- "/scratch/project_2000350/genomics/methylation/edgelist_mata.fn"

output.CG <- buildPedigree(nodelist = nodeFile, edgelist = edgeFile, cytosine = "CG", posteriorMaxFilter = 0.99)

save(output.CG, file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/MA_CG_excent.RData")
