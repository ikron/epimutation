library(data.table)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library("methimpute", lib = "/projappl/project_2000350/rpackages")

fasta.file <- "/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_lambaDNA_modified.fasta"
report.folder <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/methimpute"
input.folder <- "/scratch/project_2000350/genomics/nanopore/MAlines_5mC/bedMethyl_input"
samples <- list.files(input.folder) #Samples are the folders in the input folder
nsamples <- length(samples) #Number of samples

#Because most methylation extractor programs report only covered cytosines,we need to inflate the data to inlcude all cytosines (including non-covered sites
cytosine.positions <- extractCytosinesFromFASTA(fasta.file, contexts = c('CG','CHG','CHH'))

for(i in 1:nsamples) {
#Import modified data from Modkit 
sample <- samples[i]
path <- paste(input.folder, "/", sample, sep = "")   
data <- importBismark(path, chrom.lengths=NULL, skip = 0)
methylome <- inflateMethylome(data, cytosine.positions)
distcor <- distanceCorrelation(methylome)
fit <- estimateTransDist(distcor)
model <- callMethylation(data = methylome, transDist = fit$transDist)
exportMethylome(model, file = paste(report.folder, "/", sample, sep = ""))
}

