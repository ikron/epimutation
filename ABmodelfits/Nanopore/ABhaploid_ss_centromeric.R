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
source("/scratch/project_2000350/genomics/methylation/BOOTmodel_Haploid.R") #Bootstrap function for haploid model

#Load pedigree data
load(file = "/scratch/project_2000350/genomics/methylation/ABresults/input/AB_output_H3K9cent_CG.RData")
load(file = "/scratch/project_2000350/genomics/methylation/ABresults/input/AB_output_H3K9cent_CHG.RData")
load(file = "/scratch/project_2000350/genomics/methylation/ABresults/input/AB_output_H3K9cent_CHH.RData")

load(file = "/scratch/project_2000350/genomics/methylation/data/cytosine.prop.MA.anc.centromeric.RData") #Load cytosine proportions
#cyt.props.nano.centromeric

### CG context ###
#Final pedigree CG
pedigree.CG <- AB.output.H3K9cent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- cyt.props.nano.centromeric[2,3] #For unmethylated CG proportion in the MA ancestors
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Run the models
neutral.haploid.ss.centromeric.CG <- ABneutralHaploid(pedigree.data = pedigree.final.CG, p0uu = p0uu_in.CG,
eqp = p0uu_in.CG, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CG_centromeric_ss_nano")


boot.haploid.ss.centromeric.CG <- BOOTmodel_Haploid(model.fit = neutral.haploid.ss.centromeric.CG, Nboot = 1000, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "BOOT_ABneutralHaploid_CG_centromeric_ss_nano")

####################

### CHG context ###
#Final pedigree CHG
pedigree.CHG <- AB.output.H3K9cent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- cyt.props.nano.centromeric[4,3] #For unmethylated CHG proportion in the MA ancestors
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]


#Run the models
neutral.haploid.ss.centromeric.CHG <- ABneutralHaploid(pedigree.data = pedigree.final.CHG, p0uu = p0uu_in.CHG,
eqp = p0uu_in.CHG, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CHG_centromeric_ss_nano")


boot.haploid.ss.centromeric.CHG <- BOOTmodel_Haploid(model.fit = neutral.haploid.ss.centromeric.CHG, Nboot = 1000, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "BOOT_ABneutralHaploid_CHG_centromeric_ss_nano")


###################

### CHH context ###
#Final pedigree CHH
pedigree.CHH <- AB.output.H3K9cent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- cyt.props.nano.centromeric[6,3] #For unmethylated CHH proportion in the MA ancestors
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Run the models
neutral.haploid.ss.centromeric.CHH <- ABneutralHaploid(pedigree.data = pedigree.final.CHH, p0uu = p0uu_in.CHH,
eqp = p0uu_in.CHH, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CHH_centromeric_ss_nano")


boot.haploid.ss.centromeric.CHH <- BOOTmodel_Haploid(model.fit = neutral.haploid.ss.centromeric.CHH, Nboot = 1000, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "BOOT_ABneutralHaploid_CHH_centromeric_ss_nano")

##################

### Save the results ###

save(neutral.haploid.ss.centromeric.CG, boot.haploid.ss.centromeric.CG, neutral.haploid.ss.centromeric.CHG, boot.haploid.ss.centromeric.CHG, neutral.haploid.ss.centromeric.CHH, boot.haploid.ss.centromeric.CHH, file = "/scratch/project_2000350/genomics/methylation/ABresults/haploid/ABresults.haploid.ss.centromeric.RData")




