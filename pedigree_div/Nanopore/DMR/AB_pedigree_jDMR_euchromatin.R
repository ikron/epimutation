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

#H3K9euchromatin
nodeFile.CG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_euchromatin_CG.fn"
#run and save CG context
AB.output.jDMR.euchromatin.CG <- buildPedigreeRegions(nodelist = nodeFile.CG, edgelist = edgeFile, cytosine = "CG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.euchromatin.CG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_euchromatin_CG.RData")

nodeFile.CHH <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_euchromatin_CHH.fn"
#run and save CHH context
AB.output.jDMR.euchromatin.CHH <- buildPedigreeRegions(nodelist = nodeFile.CHH, edgelist = edgeFile, cytosine = "CHH", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.euchromatin.CHH, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_euchromatin_CHH.RData")

nodeFile.CHG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_euchromatin_CHG.fn"
#run and save CHG context
AB.output.jDMR.euchromatin.CHG <- buildPedigreeRegions(nodelist = nodeFile.CHG, edgelist = edgeFile, cytosine = "CHG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.euchromatin.CHG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_euchromatin_CHG.RData")

#excent
nodeFile.exCG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_excent_CG.fn"
#run and save CG context
AB.output.jDMR.excent.CG <- buildPedigreeRegions(nodelist = nodeFile.exCG, edgelist = edgeFile, cytosine = "CG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.excent.CG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_excent_CG.RData")

nodeFile.exCHH <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_excent_CHH.fn"
#run and save CHH context
AB.output.jDMR.excent.CHH <- buildPedigreeRegions(nodelist = nodeFile.exCHH, edgelist = edgeFile, cytosine = "CHH", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.excent.CHH, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_excent_CHH.RData")

nodeFile.exCHG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_excent_CHG.fn"
#run and save CHG context
AB.output.jDMR.excent.CHG <- buildPedigreeRegions(nodelist = nodeFile.exCHG, edgelist = edgeFile, cytosine = "CHG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.excent.CHG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_excent_CHG.RData")


#H3K27
nodeFile.H3K27CG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_H3K27_CG.fn"
#run and save CG context
AB.output.jDMR.H3K27.CG <- buildPedigreeRegions(nodelist = nodeFile.H3K27CG, edgelist = edgeFile, cytosine = "CG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.H3K27.CG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_H3K27_CG.RData")

nodeFile.H3K27CHH <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_H3K27_CHH.fn"
#run and save CHH context
AB.output.jDMR.H3K27.CHH <- buildPedigreeRegions(nodelist = nodeFile.H3K27CHH, edgelist = edgeFile, cytosine = "CHH", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.H3K27.CHH, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_H3K27_CHH.RData")

nodeFile.H3K27CHG <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/nodelist_jDMR_H3K27_CHG.fn"
#run and save CHG context
AB.output.jDMR.H3K27.CHG <- buildPedigreeRegions(nodelist = nodeFile.H3K27CHG, edgelist = edgeFile, cytosine = "CHG", posteriorMaxFilter = 0.99)
save(AB.output.jDMR.H3K27.CHG, file = "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alphabeta/AB_output_jDMR_H3K27_CHG.RData")


