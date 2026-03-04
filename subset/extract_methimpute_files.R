#rm(list=ls())
library(data.table)
library(dplyr)
library(GenomicRanges)

#coded by Rashmi (with some modifications by Ilkka)
subset.methimpute <- function(methylomes, annotation, out.dir, filt.context){
  methylome.files <- list.files(methylomes, pattern="_All.txt", full.names=TRUE)
  annotation.files <- list.files(annotation, pattern=".bed", full.names=TRUE)
  for (i1 in seq_along(methylome.files)){
    meth.name <- gsub(".*methylome_|\\_All.txt","", basename(methylome.files[i1]))
    cat(paste0("\nRunning ", meth.name, " ...\n"), sep="")
    cat("\n")
    methylome <- fread(methylome.files[i1])
    #remove Mt and chloroplast coordinates etc.
    methylome <- methylome[seqnames %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"), ] #This drops mtDNA. lambda phage, and contigs not mapping to chromosomes
    # We only keep high confidence data and filter for context. But possibility to keep all contexts
    if(filt.context != "All") { methylome <- methylome[context==filt.context & posteriorMax>=0.99,] }
    if(filt.context == "All") { methylome <- methylome[posteriorMax>=0.99,] }
    m.gr <- GRanges(seqnames=methylome$seqnames, 
                    ranges=IRanges(start=methylome$start, width=1),
                    strand=methylome$strand,
                    context=methylome$context, 
                    counts.methylated=methylome$counts.methylated,
                    counts.total=methylome$counts.total,
                    posteriorMax=methylome$posteriorMax,
                    status=methylome$status,
                    rcmethlvl=methylome$rc.meth.lvl)
    
    for (i2 in seq_along(annotation.files)){
      #annotation.file.name <- gsub(".*At_segments_|\\.bed$", "", basename(annotation.file[i2]))
      annotation.file.name <- gsub(".bed$", "", basename(annotation.files[i2]))
      cat(paste0("Running for Annotation ", annotation.file.name, " ...\n"), sep="")
      annotation.file.in <- fread(annotation.files[i2])
      annotation.file.gr <- GRanges(seqnames=annotation.file.in$V1, ranges=IRanges(start=annotation.file.in$V2, end=annotation.file.in$V3))
      overlaps <- findOverlaps(m.gr, annotation.file.gr)
      overlaps.hits <- m.gr[queryHits(overlaps)]
      out.df <- data.frame(overlaps.hits)
      out.df <- out.df[-c(3,4)]
      fwrite(x=out.df, file=paste0(out.dir, "/", meth.name, "_", filt.context, "_", annotation.file.name, "_methimpute.txt"), 
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    }
    rm(m.gr)
  }
}


path.methylome <- "/scratch/project_2000350/genomics/methylation/methimpute/all"
path.annotation <- "/scratch/project_2000350/genomics/methylation/methimpute/annotation"
out.dir <- "/scratch/project_2000350/genomics/methylation/methimpute/subsets"
subset.methimpute(methylomes=path.methylome, 
                  annotation=path.annotation, 
                  out.dir=out.dir, 
                  filt.context="All")


### Some tests, not run ###
#path.methylome <- "~/Genomics/Neurospora/methylation/methimpute/all"
#path.annotation <- "~/Genomics/Neurospora/methylation/methimpute/annotation"
#out.dir <- "~/Genomics/Neurospora/methylation/methimpute/subsets"
#subset.methimpute(methylomes=path.methylome, 
#                  annotation=path.annotation, 
#                  out.dir=out.dir, 
#                  filt.context="All")

#methylome <- fread("~/Genomics/Neurospora/methylation/methimpute/ML7G40_bismark_pe_trim_CG_excent_methimpute.txt")
