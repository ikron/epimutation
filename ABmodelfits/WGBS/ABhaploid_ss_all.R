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

#Load pedigree data
load(file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/all/MA_CG_all.RData")
load(file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/all/MA_CHG_all.RData")
load(file = "/scratch/project_2000350/genomics/methylation/methimpute/pedigreedata/all/MA_CHH_all.RData")

load(file = "/scratch/project_2000350/genomics/methylation/data/cytosine.prop.MA.anc.whole.RData") #Load cytosine proportions
#cyt.props.wgbs.whole

### CG context ###
#Final pedigree CG
pedigree.CG <- output.all.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- cyt.props.wgbs.whole[2,3] #For unmethylated CG proportion in the MA ancestors
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Run the models
neutral.haploid.ss.all.CG <- ABneutralHaploid(pedigree.data = pedigree.final.CG, p0uu = p0uu_in.CG,
eqp = p0uu_in.CG, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CG_all_ss.RData")
#In the actual analysis need to have Nstarts at least 50, here use 100 Nstarts

####################

### CHG context ###
#Final pedigree CHG
pedigree.CHG <- output.all.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- cyt.props.wgbs.whole[4,3] #For unmethylated CHG proportion in the MA ancestors
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Run the models
neutral.haploid.ss.all.CHG <- ABneutralHaploid(pedigree.data = pedigree.final.CHG, p0uu = p0uu_in.CHG,
eqp = p0uu_in.CHG, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CHG_all_ss.RData")

###################

### CHH context ###
#Final pedigree CHH
pedigree.CHH <- output.all.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- cyt.props.wgbs.whole[6,3] #For unmethylated CHH proportion in the MA ancestors
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Run the models
neutral.haploid.ss.all.CHH <- ABneutralHaploid(pedigree.data = pedigree.final.CHH, p0uu = p0uu_in.CHH,
eqp = p0uu_in.CHH, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CHH_all_ss.RData")

##################

### Save the results ###

save(neutral.haploid.ss.all.CG, neutral.haploid.ss.all.CHG, neutral.haploid.ss.all.CHH, file = "/scratch/project_2000350/genomics/methylation/ABresults/haploid/ABresults.haploid.ss.all.RData")




