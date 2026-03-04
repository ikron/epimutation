##R script to calculate divergence data for different cytosince contexts
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(data.table)

# Setting the paths
.libPaths(c("/projappl/project_2000350/rpackages", .libPaths()))
library(methimpute)
#library(jDMR)
source("/projappl/project_2000350/Genomics/jDMR_mod/jDMR_mod.R")

out.dir <- "/scratch/project_2000350/genomics/methylation/jDMR/nanopore/all"
reference <- "/scratch/project_2000350/genomics/methylation/mreference/fastafiles"
samplefile <- "/scratch/project_2000350/genomics/methylation/methimpute/control/listfiles_all_nanopore.fn"
matresults <- "/scratch/project_2000350/genomics/methylation/jDMRmatrix/nanopore/all"

#Call DMRs from methimpute files, using the grid approach
#runjDMRgrid(out.dir = out.dir, fasta.file = reference, samplefiles = samplefile, min.C = 5, genome = "Neurospora")
runMethimputeGrid(out.dir = out.dir, fasta = reference, samplefiles = samplefile, win = 100, step = 100, nCytosines = 5, mincov = 5, context = c("CG", "CHG", "CHH"), genome = "Neurospora")


#Makes matrix files for DMRs
makeDMRmatrix(samplefiles = samplefile, input.dir = out.dir, out.dir = matresults, context = c("CG", "CHG", "CHH"))

#Filter DMRs 
filterDMRmatrix(gridDMR = TRUE, data.dir = matresults)

#DMR annotation
#annotateDMRs(gff.files = c(gff.genes, gff.TE), annotation = c("gene", "TE"), input.dir = anno.dir, gff3.out = FALSE, out.dir = out.dir)


