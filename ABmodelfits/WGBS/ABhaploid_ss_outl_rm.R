##R script to calculate divergence data for different cytosince contexts
# Setting the paths
.libPaths(c("/projappl/project_2000350/rpackages", .libPaths()))

library(tidyverse)
library(data.table)
library(gtools)
library(igraph)
library(BiocParallel)
library(AlphaBeta)
#source("/scratch/project_2000350/genomics/methylation/ABmod.R") #Modified version of AlphaBeta (multicore setup needs to be different)
source("/scratch/project_2000350/genomics/methylation/ABhaploid.R") #Haploid models for AlphaBeta

#Load pedigree data, already processed
load(file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/all/MA_pedigree_CG_outl_rm.RData")
load(file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/all/MA_pedigree_CHG_outl_rm.RData")
load(file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/all/MA_pedigree_CHH_outl_rm.RData")

load(file = "/scratch/project_2000350/genomics/methylation/data/cytosine.prop.MA.anc.whole.RData") #Load cytosine proportions
#cyt.props.wgbs.whole

#Initial unmethylated proportions
p0uu_in.CG <- cyt.props.wgbs.whole[2,3]
p0uu_in.CHG <- cyt.props.wgbs.whole[4,3]
p0uu_in.CHH <- cyt.props.wgbs.whole[6,3]

### CG context ###
pedigree.all.ss.CG.outl.rm <- pedigree.all.ss.CG.outl.rm[,1:4]

#Run the models
neutral.haploid.ss.outl.rm.CG <- ABneutralHaploid(pedigree.data = pedigree.all.ss.CG.outl.rm, p0uu = p0uu_in.CG,
eqp = p0uu_in.CG, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CG_all_ss_outl_rm")
#In the actual analysis need to have Nstarts at least 50, here use 100 Nstarts

####################

### CHG context ###
pedigree.all.ss.CHG.outl.rm <- pedigree.all.ss.CHG.outl.rm[,1:4]

#Run the models
neutral.haploid.ss.outl.rm.CHG <- ABneutralHaploid(pedigree.data = pedigree.all.ss.CHG.outl.rm, p0uu = p0uu_in.CHG,
eqp = p0uu_in.CHG, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CHG_all_ss_outl_rm")

###################

### CHH context ###
pedigree.all.ss.CHH.outl.rm <- pedigree.all.ss.CHH.outl.rm[,1:4]

#Run the models
neutral.haploid.ss.outl.rm.CHH <- ABneutralHaploid(pedigree.data = pedigree.all.ss.CHH.outl.rm, p0uu = p0uu_in.CHH,
eqp = p0uu_in.CHH, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CHH_all_ss_outl_rm")

##################

### Save the results ###

save(neutral.haploid.ss.outl.rm.CG, neutral.haploid.ss.outl.rm.CHG, neutral.haploid.ss.outl.rm.CHH, file = "/scratch/project_2000350/genomics/methylation/ABresults/haploid/ABresults.haploid.ss.outl.rm.RData")




