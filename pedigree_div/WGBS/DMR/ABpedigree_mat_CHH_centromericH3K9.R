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

nodeFile.matA <- "/scratch/project_2000350/genomics/methylation/methimpute/control/nodelist_matA_centromericH3K9.fn"
edgeFile.matA <- "/scratch/project_2000350/genomics/methylation/methimpute/control/edgelist_matA.fn"

nodeFile.mata <- "/scratch/project_2000350/genomics/methylation/methimpute/control/nodelist_mata_centromericH3K9.fn"
edgeFile.mata <- "/scratch/project_2000350/genomics/methylation/methimpute/control/edgelist_mata.fn"

#nodeFile <- "/scratch/project_2000350/genomics/methylation/nodelist_matA.fn"
#edgeFile <- "/scratch/project_2000350/genomics/methylation/edgelist_matA.fn"

#nodeFile <- "/scratch/project_2000350/genomics/methylation/nodelist_mata.fn"
#edgeFile <- "/scratch/project_2000350/genomics/methylation/edgelist_mata.fn"

output.matA.centromericH3K9.CHH <- buildPedigree(nodelist = nodeFile.matA, edgelist = edgeFile.matA, cytosine = "CHH", posteriorMaxFilter = 0.99)

output.mata.centromericH3K9.CHH <- buildPedigree(nodelist = nodeFile.mata, edgelist = edgeFile.mata, cytosine = "CHH", posteriorMaxFilter = 0.99)

save(output.matA.centromericH3K9.CHH, output.mata.centromericH3K9.CHH, file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/MA_mat_CHH_centromericH3K9.RData")
