##R script to calculate divergence data for different cytosince contexts
# Setting the paths
.libPaths(c("/projappl/project_2000350/rpackages", .libPaths()))

library(tidyverse)
library(data.table)
library(gtools)
library(igraph)
library(BiocParallel)
#library(AlphaBeta)
source("/scratch/project_2000350/genomics/methylation/buildPedigree_Regions.R")

multicoreParam <- MulticoreParam(workers=4) #Needed to set up multicore parameters


edgeFile <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/edgelist.fn"


nodeFile.CG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_CG.fn"
#run and save CG context
AB.output.jDMR.CG <- buildPedigreeRegions(nodelist = nodeFile.CG, edgelist = edgeFile, cytosine = "CG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.CG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_CG.RData")

nodeFile.CHH <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_CHH.fn"
#run and save CHH context
AB.output.jDMR.CHH <- buildPedigreeRegions(nodelist = nodeFile.CHH, edgelist = edgeFile, cytosine = "CHH", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.CHH, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_CHH.RData")

nodeFile.CHG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_CHG.fn"
#run and save CHG context
AB.output.jDMR.CHG <- buildPedigreeRegions(nodelist = nodeFile.CHG, edgelist = edgeFile, cytosine = "CHG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.CHG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_CHG.RData")
