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

#H3K9cent
nodeFile.CG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_H3K9cent_CG.fn"
#run and save CG context
AB.output.jDMR.H3K9cent.CG <- buildPedigreeRegions(nodelist = nodeFile.CG, edgelist = edgeFile, cytosine = "CG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.H3K9cent.CG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_H3K9cent_CG.RData")

nodeFile.CHH <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_H3K9cent_CHH.fn"
#run and save CHH context
AB.output.jDMR.H3K9cent.CHH <- buildPedigreeRegions(nodelist = nodeFile.CHH, edgelist = edgeFile, cytosine = "CHH", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.H3K9cent.CHH, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_H3K9cent_CHH.RData")

nodeFile.CHG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_H3K9cent_CHG.fn"
#run and save CHG context
AB.output.jDMR.H3K9cent.CHG <- buildPedigreeRegions(nodelist = nodeFile.CHG, edgelist = edgeFile, cytosine = "CHG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.H3K9cent.CHG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_H3K9cent_CHG.RData")

#H3K9excent
nodeFile.exCG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_H3K9excent_CG.fn"
#run and save CG context
AB.output.jDMR.H3K9excent.CG <- buildPedigreeRegions(nodelist = nodeFile.exCG, edgelist = edgeFile, cytosine = "CG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.H3K9excent.CG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_H3K9excent_CG.RData")

nodeFile.exCHH <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_H3K9excent_CHH.fn"
#run and save CHH context
AB.output.jDMR.H3K9excent.CHH <- buildPedigreeRegions(nodelist = nodeFile.exCHH, edgelist = edgeFile, cytosine = "CHH", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.H3K9excent.CHH, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_H3K9excent_CHH.RData")

nodeFile.exCHG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_H3K9excent_CHG.fn"
#run and save CHG context
AB.output.jDMR.H3K9excent.CHG <- buildPedigreeRegions(nodelist = nodeFile.exCHG, edgelist = edgeFile, cytosine = "CHG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.H3K9excent.CHG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_H3K9excent_CHG.RData")

