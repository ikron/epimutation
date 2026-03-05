##R script to calculate divergence data for different cytosince contexts
# Setting the paths
.libPaths(c("/projappl/project_2000350/rpackages", .libPaths()))

library(magrittr)
library(tidyverse)
library(rtracklayer)
library(data.table)
library(methimpute)
#library(jDMR)
source("/projappl/project_2000350/Genomics/jDMR_mod/jDMR_mod.R")

out.dir <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/jDMR_subsets"
reference <- "/scratch/project_2000350/genomics/methylation/mreference/fastafiles"
samplefile_H3K27 <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/list_file_H3K27.fn"
samplefile_euchromatin <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/list_file_euchromatin.fn"
samplefile_excent <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/list_files/list_file_excent.fn"

#matresults <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/jDMRmatrix"

#Call DMRs from methimpute files, using the grid approach
runMethimputeGrid(out.dir = out.dir, fasta = reference, samplefiles = samplefile_H3K27, win = 100, step = 100, nCytosines = 5, mincov = 5, context = c("CG", "CHG", "CHH"), genome = "Neurospora")

runMethimputeGrid(out.dir = out.dir, fasta = reference, samplefiles = samplefile_euchromatin, win = 100, step = 100, nCytosines = 5, mincov = 5, context = c("CG", "CHG", "CHH"), genome = "Neurospora")

runMethimputeGrid(out.dir = out.dir, fasta = reference, samplefiles = samplefile_excent, win = 100, step = 100, nCytosines = 5, mincov = 5, context = c("CG", "CHG", "CHH"), genome = "Neurospora")



#Makes matrix files for DMRs
#makeDMRmatrix(samplefiles = samplefile, input.dir = out.dir, out.dir = matresults, context = c("CG", "CHG", "CHH"))

#Filter DMRs 
#filterDMRmatrix(gridDMR = TRUE, data.dir = matresults)

#DMR annotation
#annotateDMRs(gff.files = c(gff.genes, gff.TE), annotation = c("gene", "TE"), input.dir = anno.dir, gff3.out = FALSE, out.dir = out.dir)
