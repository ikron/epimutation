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
load(file = "/scratch/project_2000350/genomics/methylation/ABresults/input/AB_output_jDMR_H3K9cent_CG.RData")
load(file = "/scratch/project_2000350/genomics/methylation/ABresults/input/AB_output_jDMR_H3K9cent_CHG.RData")
load(file = "/scratch/project_2000350/genomics/methylation/ABresults/input/AB_output_jDMR_H3K9cent_CHH.RData")

load(file = "/scratch/project_2000350/genomics/methylation/data/DMR.prop.MA.anc.cent.RData") #Load cytosine proportions
#DMR.props.nano.cent

### CG context ###
#Final pedigree CG
pedigree.CG <- AB.output.jDMR.H3K9cent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- DMR.props.nano.cent[2,3] #For unmethylated CG proportion in the MA ancestors
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Run the models
neutral.haploid.DMR.centromeric.CG <- ABneutralHaploid(pedigree.data = pedigree.final.CG, p0uu = p0uu_in.CG,
eqp = p0uu_in.CG, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CG_centromeric_DMR_nano")


boot.haploid.DMR.centromeric.CG <- BOOTmodel_Haploid(model.fit = neutral.haploid.DMR.centromeric.CG, Nboot = 1000, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "BOOT_ABneutralHaploid_CG_centromeric_DMR_nano")

####################

### CHG context ###
#Final pedigree CHG
pedigree.CHG <- AB.output.jDMR.H3K9cent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- DMR.props.nano.cent[4,3] #For unmethylated CHG proportion in the MA ancestors
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]


#Run the models
neutral.haploid.DMR.centromeric.CHG <- ABneutralHaploid(pedigree.data = pedigree.final.CHG, p0uu = p0uu_in.CHG,
eqp = p0uu_in.CHG, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CHG_centromeric_DMR_nano")


boot.haploid.DMR.centromeric.CHG <- BOOTmodel_Haploid(model.fit = neutral.haploid.DMR.centromeric.CHG, Nboot = 1000, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "BOOT_ABneutralHaploid_CHG_centromeric_DMR_nano")


###################

### CHH context ###
#Final pedigree CHH
pedigree.CHH <- AB.output.jDMR.H3K9cent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- DMR.props.nano.cent[6,3] #For unmethylated CHH proportion in the MA ancestors
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Run the models
neutral.haploid.DMR.centromeric.CHH <- ABneutralHaploid(pedigree.data = pedigree.final.CHH, p0uu = p0uu_in.CHH,
eqp = p0uu_in.CHH, eqp.weight = 0.001, Nstarts = 100, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "ABneutralHaploid_CHH_centromeric_DMR_nano")


boot.haploid.DMR.centromeric.CHH <- BOOTmodel_Haploid(model.fit = neutral.haploid.DMR.centromeric.CHH, Nboot = 1000, out.dir = "/scratch/project_2000350/genomics/methylation/ABresults/haploid", out.name = "BOOT_ABneutralHaploid_CHH_centromeric_DMR_nano")

##################

### Save the results ###

save(neutral.haploid.DMR.centromeric.CG, boot.haploid.DMR.centromeric.CG, neutral.haploid.DMR.centromeric.CHG, boot.haploid.DMR.centromeric.CHG, neutral.haploid.DMR.centromeric.CHH, boot.haploid.DMR.centromeric.CHH, file = "/scratch/project_2000350/genomics/methylation/ABresults/haploid/ABresults.haploid.DMR.centromeric.RData")




