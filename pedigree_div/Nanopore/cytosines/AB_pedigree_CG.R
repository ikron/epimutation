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

nodeFile <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/nodelist.fn"
edgeFile <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/edgelist.fn"

#run and save CG context
AB.output.CG <- buildPedigree(nodelist = nodeFile, edgelist = edgeFile, cytosine = "CG", posteriorMaxFilter = 0.99)
save(AB.output.CG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_CG.RData")

#run and save CHH context
#AB.output.CHH <- buildPedigree(nodelist = nodeFile, edgelist = edgeFile, cytosine = "CHH", posteriorMaxFilter = 0.99)
#save(AB.output.CHH, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_CHH.RData")

#run and save CHG context
#AB.output.CHG <- buildPedigree(nodelist = nodeFile, edgelist = edgeFile, cytosine = "CHG", posteriorMaxFilter = 0.99)
#save(AB.output.CHG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_CHG.RData")
