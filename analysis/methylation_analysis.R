### Analysis of DNA methylation data from the MA-experiment
library(tidyverse)
library(data.table)
#library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(latex2exp)
library(ape)
theme_set(theme_cowplot())

library(AlphaBeta)

source("~/Documents/tutkijatohtori/scripts/jDMR_functions.R") #Load functions dealing with DMRs
source("~/Documents/tutkijatohtori/scripts/ABhaploid.R") #Haploid versions of the AlphaBeta models
source("~/Documents/tutkijatohtori/scripts/BOOTmodel_Haploid.R") #Bootstrap function for haploid model

### * Testing that .bed files are OK

load("ForChrPlot.RData")

euchr <- read.table("./annotation/2489.euchromatin.bed", header = F, sep = "\t")
colnames(euchr) <- c("Chromosome", "start", "end")
euchr <- filter(euchr, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

excent.h3k9 <- read.table("./annotation/H3K9_ex_cent.bed", header = F, sep = "\t")
colnames(excent.h3k9)[1:3] <- c("Chromosome", "start", "end")
excent.h3k9 <- filter(excent.h3k9, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

excent2.h3k9 <- read.table("./annotation/H3K9_ex_cent2.bed", header = F, sep = "\t")
colnames(excent2.h3k9)[1:3] <- c("Chromosome", "start", "end")
excent2.h3k9 <- filter(excent2.h3k9, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

centromeric.h3k9 <- read.table("./annotation/centromericH3K9.bed", header = F, sep = "\t")
colnames(centromeric.h3k9)[1:3] <- c("Chromosome", "start", "end")

excent <- read.table("./annotation/excent.bed", header = F, sep = "\t")
colnames(excent)[1:3] <- c("Chromosome", "start", "end")
excent <- filter(excent, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

h3k27.exh3k9 <- read.table("./annotation/H3K27_exH3K9.bed", header = F, sep = "\t")
colnames(h3k27.exh3k9)[1:3] <- c("Chromosome", "start", "end")
h3k27.exh3k9 <- filter(h3k27.exh3k9, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Plotting Neurospora genes 
genes.plot <- genes
colnames(genes.plot)[1] <- "Chromosome"


Chr_test.plot <-  ggplot(chr.sizes, aes()) +
     geom_rect( aes(xmin = Start, xmax = End, ymin = 0, ymax = 5, fill = "Chromosome"), alpha = 0.5) +
     geom_rect(data = cent, aes(xmin = start, xmax = end, ymin = -1, ymax = -3, fill = "Centromeric")) +
     geom_rect(data = h3k27.plot, aes(xmin = Start, xmax = End, ymax = -4, ymin = -6, fill = "H3K27me3")) +
     geom_rect(data = h3k9.plot, aes(xmin = Start, xmax = End, ymax = -7, ymin = -9, fill = "H3K9me3")) +
     geom_rect(data = euchr, aes(xmin = start, xmax = end, ymax = -10, ymin = -12, fill = "Euchromatin")) +
     geom_rect(data = excent.h3k9, aes(xmin = start, xmax = end, ymax = -13, ymin = -15, fill = "H3K9_ex_cent")) +
     geom_rect(data = excent2.h3k9, aes(xmin = start, xmax = end, ymax = -16, ymin = -18, fill = "H3K9_ex_cent2")) +
     geom_rect(data = centromeric.h3k9, aes(xmin = start, xmax = end, ymax = -19, ymin = -21, fill = "Centromeric_H3K9")) +
      geom_rect(data = excent, aes(xmin = start, xmax = end, ymax = -23, ymin = -25, fill = "Excent")) +
     geom_rect(data = h3k27.exh3k9, aes(xmin = start, xmax = end, ymax = -26, ymin = -28, fill = "H3K27_exH3K9")) +      
     #geom_rect(data = h3k36.plot, aes(xmin = Start, xmax = End, ymax = -10, ymin = -13, fill = "H3K36me")) +
     geom_rect(data = genes.plot, aes(xmin = start, xmax = end, ymax = 5, ymin = 0, fill = "Gene")) +
     #geom_rect(data = h3k9.ex.cent.plot, aes(xmin = Start, xmax = End, ymax = -10, ymin = -12, fill = "hotpink")) +
     #geom_tile(data = genes.plot, aes(x = Position, y = 2.5, height = 5, width = 2), colour = "black", alpha = 0.5) +    
     xlab("Coordinate (bp)") +
     ylab("") +
     scale_fill_manual(breaks = c("Chromosome", "Centromeric", "H3K27me3", "H3K9me3", "Euchromatin", "H3K9_ex_cent", "H3K9_ex_cent2", "Centromeric_H3K9", "Excent", "H3K27_exH3K9"), values = c(Chromosome = "deepskyblue", Centromeric = "grey", H3K27me3 = "blue", H3K9me3 = "red", Euchromatin = "green", H3K9_ex_cent = "hotpink", H3K9_ex_cent2 = "darkred", Centromeric_H3K9 = "orange", Excent = "darkblue", H3K27_exH3K9 = "darkgreen")) +
     scale_x_continuous(expand = expansion(mult = c(0, 0.05)), breaks = c(0,2.5e6,5e6,7.5e6,1e7), labels = c(TeX("$0$"), TeX("$2.5 \\times 10^6$"), TeX("$5 \\times 10^6$"), TeX("$7.5 \\times 10^6$"), TeX("$1 \\times 10^7$")), limits = c(0, 1e7))  +
     guides(fill = guide_legend(nrow = 1)) +
     facet_grid(Chromosome ~ ., labeller = labeller(Chromosome = labels)) +
     theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "top", legend.justification = c(0.5, 0.5), legend.title = element_blank())

#Save the domains
save(cent, euchr, chr.sizes, excent, excent2.h3k9, h3k27.exh3k9, h3k27.plot, h3k36.plot, h3k9.plot, file = "./data/domains.RData")


### * Analysis of epimutation rates using AlphaBeta

### ** Build pedigrees and pairwise divergence matrices using AB

### NOTE ###
### Methylation divergences have been calculated using thr CSC cluster. Refer to the separate script files in ./pedigree_div ###

#Some tests to check that calculating pedigree divergence works
#This takes a long time
#output <- buildPedigree(nodelist = "nodelist.fn", edgelist = "edgelist.fn", cytosine = "CG", posteriorMaxFilter = 0.99)

#Build pedigrees for matA MA-lines
#output.matA <- buildPedigree(nodelist = "nodelist_matA.fn", edgelist = "edgelist_matA.fn", cytosine = "CG", posteriorMaxFilter = 0.99)
#output.matA.CHG <- buildPedigree(nodelist = "nodelist_matA.fn", edgelist = "edgelist_matA.fn", cytosine = "CHG", posteriorMaxFilter = 0.99)
#output.matA.CHH <- buildPedigree(nodelist = "nodelist_matA.fn", edgelist = "edgelist_matA.fn", cytosine = "CHH", posteriorMaxFilter = 0.99)
#output.mata.CG <- buildPedigree(nodelist = "nodelist_mata.fn", edgelist = "edgelist_mata.fn", cytosine = "CG", posteriorMaxFilter = 0.99)

### ** Running the AB models

## Load the data
#Output has been produced by cluster

#load(file = "./pedigree_data/MA_pedigree_data_CG_all.RData")
##load(file = "./pedigree_data/MA_pedigree_data_CHG_all.RData")
#load(file = "./pedigree_data/MA_pedigree_data_CHH_all.RData")
load(file = "./pedigree_data/all/MA_CG_all.RData")
load(file = "./pedigree_data/all/MA_CHG_all.RData")
load(file = "./pedigree_data/all/MA_CHH_all.RData")

#Final pedigree CG
pedigree.CG <- output.all.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.all.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.all.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.all.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.all.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.all.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]


### *** Run models for all samples for single cytosines

### Load the results for all different models run on the CSC cluster

load("./modeldata/haploid/all/ABresults.haploid.ss.all.RData")
#Variables containing the model results are:
#neutral.haploid.ss.all.CG
#neutral.haploid.ss.all.CHG
#neutral.haploid.ss.all.CHH

### Model comparisons ###
#Comparing the neutral models to the null model of no accumulation
comp.out.ss.CG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CG_all_ss.Rdata", pedigree.null = "./modeldata/ABnull_CG_global_estimates.Rdata")
comp.out.ss.CG$Ftest #p-value = 1.223804e-08

comp.out.ss.CHG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CHG_all_ss.Rdata", pedigree.null = "./modeldata/ABnull_CHG_global_estimates.Rdata")
comp.out.ss.CHG$Ftest #p-value = 4.178146e-11

comp.out.ss.CHH <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CHH_all_ss.Rdata", pedigree.null = "./modeldata/ABnull_CHH_global_estimates.Rdata")
comp.out.ss.CHH$Ftest #p-value = 7.157985e-13

### End model comparisons ###

### Plotting ###
#Make plot dataframe
plotdata <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

#Theoretical fits
  theory.fit.data <- neutral.haploid.ss.all.CG$for.fit.plot
  theory.fits <- c(neutral.haploid.ss.all.CG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.CG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.ss.all.CHG$for.fit.plot
  theory.fits <- c(neutral.haploid.ss.all.CHG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.CHG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.ss.all.CHH$for.fit.plot
  theory.fits <- c(neutral.haploid.ss.all.CHH$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.CHH <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHH", length(theory.fits)))

 theoryFits <- rbind(theory.fit.CG, theory.fit.CHG, theory.fit.CHH)


#Fits with outlier samples removed
#load the theory fits with outlier samples removed
load("./data/theoryFits.outlrm.RData")
#theoryFits.ss.outlrm contains the data

##Fitting some phenomenological models here
ss.CG <- filter(plotdata, context == "CG")

e1 <- D.value ~ a / (1+exp(-(b+c*dt))) #Fitting logistic growth
m1 <- nls(D.value ~ SSlogis(dt, Asym, xmid, scal), data = ss.CG)
ss.CG$pred <- predict(m1)

ones <- rep(1, nrow(ss.CG))
e2 <- D.value ~ a * ones

m0 <-  nls(e2, data = ss.CG, start = list(a = 1))
ss.CG$prednull <- predict(m0)

ss.CHG <- filter(plotdata, context == "CHG")

e1 <- D.value ~ a / (1+exp(-(b+c*dt)))
m1.2 <- nls(D.value ~ SSlogis(dt, Asym, xmid, scal), data = ss.CHG)
ss.CHG$pred <- predict(m1.2)

ones <- rep(1, nrow(ss.CHG))
e2 <- D.value ~ a * ones

m0.2 <-  nls(e2, data = ss.CHG, start = list(a = 1))
ss.CHG$prednull <- predict(m0.2)

ss.CHH <- filter(plotdata, context == "CHH")

e1 <- D.value ~ a / (1+exp(-(b+c*dt)))
m1.3 <- nls(D.value ~ SSlogis(dt, Asym, xmid, scal), data = ss.CHH)
ss.CHH$pred <- predict(m1.3)

ones <- rep(1, nrow(ss.CHH))
e2 <- D.value ~ a * ones

m0.3 <-  nls(e2, data = ss.CHH, start = list(a = 1))
ss.CHH$prednull <- predict(m0.3)

anova(m1, m0) #Compare DMR logistic and null models, CG # p = 2.121e-12

anova(m1.2, m0.2) #p = 2.995e-14

anova(m1.3, m0.3) #p = 4.463e-16

#Combine all data for plotting
plotdata <- cbind(plotdata, rbind(ss.CG[,7:8], ss.CHG[,7:8], ss.CHH[,7:8]))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methsites.plot <- ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    geom_line(data = theoryFits, aes(x = dt, y = divsim), col = "red" ) +
    geom_line(data = theoryFits.ss.outlrm, aes(x = dt, y = divsim), col = "orange" ) +    
    #geom_smooth() +
    geom_line(aes(x = dt, y = pred), col = "blue") +
    geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +        
    ggtitle("Single sites") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

save_plot(filename = "methsites_2026Feb.pdf", methsites.plot, base_height = 4, base_width = 4*1.5*1.618)



### *** Run models for all samples for DMRs

#Load data
#load(file = "./pedigree_data/MA_pedigree_data_regions_CG_all.RData")
#load(file = "./pedigree_data/MA_pedigree_data_regions_CHG_all.RData")
#load(file = "./pedigree_data/MA_pedigree_data_regions_CHH_all.RData")
load(file = "./pedigree_data/DMR_all/MA_DMR_all_CG.RData")
load(file = "./pedigree_data/DMR_all/MA_DMR_all_CHG.RData")
load(file = "./pedigree_data/DMR_all/MA_DMR_all_CHH.RData")

#Final pedigree CG
pedigree.CG <- output.DMR.all.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.DMR.all.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.DMR.all.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.DMR.all.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.DMR.all.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.DMR.all.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

## Load the data
#Output has been produced by cluster
load("./modeldata/haploid/all/ABresults.haploid.DMR.all.RData")
#Variables are
#neutral.haploid.DMR.all.CG
#neutral.haploid.DMR.all.CHG
#neutral.haploid.DMR.all.CHH

#Null models
outputABnull.CG.DMR <- ABnull(pedigree.data = pedigree.final.CG, out.dir = getwd(), out.name = "./modeldata/ABnull_DMR_CG_all")

outputABnull.CHG.DMR <- ABnull(pedigree.data = pedigree.final.CHG, out.dir = getwd(), out.name = "./modeldata/ABnull_DMR_CHG_all")

outputABnull.CHH.DMR <- ABnull(pedigree.data = pedigree.final.CHH, out.dir = getwd(), out.name = "./modeldata/ABnull_DMR_CHH_all")

### Model comparisons
comp.out.CG.DMR <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CG_all_DMR.Rdata", pedigree.null = "./modeldata/ABnull_DMR_CG_all.Rdata")
comp.out.CG.DMR$Ftest
#p < 2.2 * 10^-16

comp.out.CHG.DMR <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CHG_all_DMR.Rdata", pedigree.null = "./modeldata/ABnull_DMR_CHG_all.Rdata")
comp.out.CHG.DMR$Ftest
#p < 2.2 * 10^-16

comp.out.CHH.DMR <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CHH_all_DMR.Rdata", pedigree.null = "./modeldata/ABnull_DMR_CHH_all.Rdata")
comp.out.CHH.DMR$Ftest
#p < 2.2 * 10^-16

#Plotting
#Make plot dataframe
plotdata.DMR <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

#neutral.model.output.CG.DMR <- dget("./modeldata/ABneutralSOMA_CG_all_DMR.RData")
#neutral.model.output.CHG.DMR <- dget("./modeldata/ABneutralSOMA_CHG_all_DMR.RData")
#neutral.model.output.CHH.DMR <- dget("./modeldata/ABneutralSOMA_CHH_all_DMR.RData")


#Theoretical fits
  theory.fit.data.DMR <- neutral.haploid.DMR.all.CG$for.fit.plot
  theory.fits.DMR <- c(neutral.haploid.DMR.all.CG$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.haploid.DMR.all.CHG$for.fit.plot
  theory.fits.DMR <- c(neutral.haploid.DMR.all.CHG$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.haploid.DMR.all.CHH$for.fit.plot
  theory.fits.DMR <- c(neutral.haploid.DMR.all.CHH$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHH.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHH", length(theory.fits.DMR)))

 theoryFits.DMR <- rbind(theory.fit.CG.DMR, theory.fit.CHG.DMR, theory.fit.CHH.DMR)

##Fitting some phenomenological models here
DMR.CG <- filter(plotdata.DMR, context == "CG")

e1 <- D.value ~ a / (1+exp(-(b+c*dt))) #Fitting logistic growth
m1 <- nls(D.value ~ SSlogis(dt, Asym, xmid, scal), data = DMR.CG)
DMR.CG$pred <- predict(m1)

ones <- rep(1, nrow(DMR.CG))
e2 <- D.value ~ a * ones

m0 <-  nls(e2, data = DMR.CG, start = list(a = 1))
DMR.CG$prednull <- predict(m0)

DMR.CHG <- filter(plotdata.DMR, context == "CHG")

e1 <- D.value ~ a / (1+exp(-(b+c*dt)))
m1.2 <- nls(D.value ~ SSlogis(dt, Asym, xmid, scal), data = DMR.CHG)
DMR.CHG$pred <- predict(m1.2)

ones <- rep(1, nrow(DMR.CHG))
e2 <- D.value ~ a * ones

m0.2 <-  nls(e2, data = DMR.CHG, start = list(a = 1))
DMR.CHG$prednull <- predict(m0.2)

DMR.CHH <- filter(plotdata.DMR, context == "CHH")

e1 <- D.value ~ a / (1+exp(-(b+c*dt)))
m1.3 <- nls(D.value ~ SSlogis(dt, Asym, xmid, scal), data = DMR.CHH)
DMR.CHH$pred <- predict(m1.3)

ones <- rep(1, nrow(DMR.CHH))
e2 <- D.value ~ a * ones

m0.3 <-  nls(e2, data = DMR.CHH, start = list(a = 1))
DMR.CHH$prednull <- predict(m0.3)

anova(m1, m0) #Compare DMR logistic and null models, CG

anova(m1.2, m0.2)

anova(m1.3, m0.3)

#All p < 2.2e-16

#Combine all data for plotting
plotdata.DMR <- cbind(plotdata.DMR, rbind(DMR.CG[,7:8], DMR.CHG[,7:8], DMR.CHH[,7:8]))



myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.plot <- ggplot(plotdata.DMR, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "red" ) +
    geom_line(data = theoryFits.DMR.outlrm, aes(x = dt, y = divsim), col = "orange" ) +    
    geom_line(aes(x = dt, y = pred), col = "blue") +
    geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +    
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

save_plot(filename = "methDMR.pdf", methDMR.plot, base_height = 4, base_width = 4*1.5*1.618)

### Make a plot with both single sites and DMRs (for poster) ##
    
finalplot <- plot_grid(methsites.plot, methDMR.plot, nrow = 2, labels = c("A", "B"))
save_plot(filename = "./epigenet_ms/fig/methdiv_both.pdf", finalplot, base_height = 8, base_width = 4 * 1.5 * 1.618)



ggplot(koe, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    geom_line(aes(x = dt, y = pred), col = "red") +
    geom_line(aes(x = dt, y = prednull), col = "red", lty = "dashed") 

anova(m1, m0)

##############################################################


### *** Run models for all samples excluding centromeric regions

### **** Run models for single cytosines

#Load the pedigree data
load(file = "./pedigree_data/excent/MA_CG_excent.RData")
output.excent.CG <- output.CG
load(file = "./pedigree_data/excent/MA_CHG_excent.RData") #output.excent.CHG
load(file = "./pedigree_data/excent/MA_CHH_excent.RData") #output.excent.CHH

#Final pedigree CG
pedigree.CG <- output.excent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.excent.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.excent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.excent.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.excent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.excent.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

#Theoretical fits
  theory.fit.data <- neutral.model.output.CG$for.fit.plot
  theory.fits <- c(neutral.model.output.CG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.CG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CG", length(theory.fits)))

  theory.fit.data <- neutral.model.output.CHG$for.fit.plot
  theory.fits <- c(neutral.model.output.CHG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.CHG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHG", length(theory.fits)))

  theory.fit.data <- neutral.model.output.CHH$for.fit.plot
  theory.fits <- c(neutral.model.output.CHH$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.CHH <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHH", length(theory.fits)))

 theoryFits <- rbind(theory.fit.CG, theory.fit.CHG, theory.fit.CHH)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methsites.plot <- ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    geom_line(data = theory.fit.CG, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("Single sites") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

#Theoretical fit for single sites is very poor
### Final plot with DMRs included in the next section


### **** Run models for DMRs

#Load data
load(file = "./pedigree_data/DMR_excent/MA_DMR_excent_CG.RData")
DMR.excent.CG <- output.DMR.excent.CG
load(file = "./pedigree_data/DMR_excent/MA_DMR_excent_CHG.RData")
DMR.excent.CHG <- output.DMR.excent.CG
load(file = "./pedigree_data/DMR_excent/MA_DMR_excent_CHH.RData")
DMR.excent.CHH <- output.DMR.excent.CG

#Final pedigree CG
pedigree.CG <- DMR.excent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- DMR.excent.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- DMR.excent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- DMR.excent.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- DMR.excent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- DMR.excent.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]


#Plotting
#Make plot dataframe
plotdata.DMR <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

#Theoretical fits
  theory.fit.data.DMR <- neutral.model.output.CG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHH.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHH.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHH.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHH", length(theory.fits.DMR)))

 theoryFits.DMR <- rbind(theory.fit.CG.DMR, theory.fit.CHG.DMR, theory.fit.CHH.DMR)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.plot <- ggplot(plotdata.DMR, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

save_plot(filename = "methDMR_excent.pdf", methDMR.plot, base_height = 4, base_width = 4*1.5*1.618)

### Make a plot with both single sites and DMRs (for poster) ##
    
finalplot <- plot_grid(methsites.plot, methDMR.plot, nrow = 2, labels = c("A", "B"))
save_plot(filename = "./epigenet_ms/fig/methdiv_both_excent.pdf", finalplot, base_height = 8, base_width = 4 * 1.5 * 1.618)

##############################################################






### *** Run models for all samples for centromeric H3K9 regions

### **** Run models for single cytosines

#Load the pedigree data
load(file = "./pedigree_data/centromericH3K9/MA_CG_centromericH3K9.RData") #output.centromericH3K9.CG
load(file = "./pedigree_data/centromericH3K9/MA_CHG_centromericH3K9.RData")
load(file = "./pedigree_data/centromericH3K9/MA_CHH_centromericH3K9.RData") 

#Final pedigree CG
pedigree.CG <- output.centromericH3K9.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.centromericH3K9.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.centromericH3K9.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.centromericH3K9.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.centromericH3K9.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.centromericH3K9.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.centH3K9 <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))



#Theoretical fits
  theory.fit.data <- neutral.model.output.CG$for.fit.plot
  theory.fits <- c(neutral.model.output.CG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.CG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CG", length(theory.fits)))

  theory.fit.data <- neutral.model.output.CHG$for.fit.plot
  theory.fits <- c(neutral.model.output.CHG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.CHG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHG", length(theory.fits)))

  theory.fit.data <- neutral.model.output.CHH$for.fit.plot
  theory.fits <- c(neutral.model.output.CHH$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.CHH <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHH", length(theory.fits)))

 theoryFits <- rbind(theory.fit.CG, theory.fit.CHG, theory.fit.CHH)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methsites.plot <- ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    geom_line(data = theory.fit.CG, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("Single sites") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

#Theoretical fit for single sites is very poor
### Final plot with DMRs included in the next section



### **** Run models for DMRs

#Load data
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_CG.RData")
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_CHG.RData")
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_CHH.RData")


#Final pedigree CG
pedigree.CG <- output.DMR.centromericH3K9.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.DMR.centromericH3K9.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.DMR.centromericH3K9.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.DMR.centromericH3K9.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.DMR.centromericH3K9.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.DMR.centromericH3K9.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.DMR.centH3K9 <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))


#Theoretical fits
  theory.fit.data.DMR <- neutral.model.output.CG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHH.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHH.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHH.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHH", length(theory.fits.DMR)))

 theoryFits.DMR <- rbind(theory.fit.CG.DMR, theory.fit.CHG.DMR, theory.fit.CHH.DMR)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.centH3K9.plot <- ggplot(plotdata.DMR, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

#save_plot(filename = "methDMR_centromericH3K9.pdf", methDMR.plot, base_height = 4, base_width = 4*1.5*1.618)

### Make a plot with both single sites and DMRs (for poster) ##
    
finalplot <- plot_grid(methsites.plot, methDMR.plot, nrow = 2, labels = c("A", "B"))
save_plot(filename = "./epigenet_ms/fig/methdiv_both_centromericH3K9.pdf", finalplot, base_height = 8, base_width = 4 * 1.5 * 1.618)

##############################################################






### *** Run models for all samples for H3K9 excluding centromeric

### **** Run models for single cytosines
#Load the pedigree data
load(file = "./pedigree_data/H3K9_ex_cent/MA_CG_H3K9_ex_cent.RData") #output.centromericH3K9.CG
load(file = "./pedigree_data/H3K9_ex_cent/MA_CHG_H3K9_ex_cent.RData")
load(file = "./pedigree_data/H3K9_ex_cent/MA_CHH_H3K9_ex_cent.RData") 

#Final pedigree CG
pedigree.CG <- output.H3K9_excent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.H3K9_excent.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.H3K9_excent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.H3K9_excent.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.H3K9_excent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.H3K9_excent.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.H3K9 <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

### **** Run models for DMRs

#Load data
load(file = "./pedigree_data/DMR_H3K9_ex_cent/MA_DMR_H3K9_ex_cent_CG.RData")
load(file = "./pedigree_data/DMR_H3K9_ex_cent/MA_DMR_H3K9_ex_cent_CHG.RData")
load(file = "./pedigree_data/DMR_H3K9_ex_cent/MA_DMR_H3K9_ex_cent_CHH.RData")


#Final pedigree CG
pedigree.CG <- output.DMR.H3K9_ex_cent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.DMR.H3K9_ex_cent.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.DMR.H3K9_ex_cent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.DMR.H3K9_ex_cent.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.DMR.H3K9_ex_cent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.DMR.H3K9_ex_cent.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.DMR.H3K9 <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))




#Theoretical fits
  theory.fit.data.DMR <- neutral.model.output.CG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHH.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHH.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHH.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHH", length(theory.fits.DMR)))

 theoryFits.DMR <- rbind(theory.fit.CG.DMR, theory.fit.CHG.DMR, theory.fit.CHH.DMR)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.H3K9_ex_cent.plot <- ggplot(plotdata.DMR, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

#save_plot(filename = "methDMR_H3K9_ex_cent.pdf", methDMR.plot, base_height = 4, base_width = 4*1.5*1.618)



### *** Run models for all samples H3K27 regions

### **** Run models for single cytosines

load(file = "./pedigree_data/H3K27_ex_H3K9/MA_CG_H3K27_exH3K9.RData") #output.centromericH3K9.CG
load(file = "./pedigree_data/H3K27_ex_H3K9/MA_CHG_H3K27_exH3K9.RData")
load(file = "./pedigree_data/H3K27_ex_H3K9/MA_CHH_H3K27_exH3K9.RData") 

#Final pedigree CG
pedigree.CG <- output.H3K27_exH3K9.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.H3K27_exH3K9.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.H3K27_exH3K9.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.H3K27_exH3K9.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.H3K27_exH3K9.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.H3K27_exH3K9.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.H3K27 <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

### **** Run models for DMRs


#Load data
load(file = "./pedigree_data/DMR_H3K27_ex_H3K9/MA_DMR_H3K27_ex_H3K9_CG.RData")
load(file = "./pedigree_data/DMR_H3K27_ex_H3K9/MA_DMR_H3K27_ex_H3K9_CHG.RData")
load(file = "./pedigree_data/DMR_H3K27_ex_H3K9/MA_DMR_H3K27_ex_H3K9_CHH.RData")


#Final pedigree CG
pedigree.CG <- output.DMR.H3K27_ex_H3K9.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.DMR.H3K27_ex_H3K9.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.DMR.H3K27_ex_H3K9.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.DMR.H3K27_ex_H3K9.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.DMR.H3K27_ex_H3K9.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.DMR.H3K27_ex_H3K9.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.DMR.H3K27 <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))



#Theoretical fits
  theory.fit.data.DMR <- neutral.model.output.CG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHH.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHH.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHH.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHH", length(theory.fits.DMR)))

 theoryFits.DMR <- rbind(theory.fit.CG.DMR, theory.fit.CHG.DMR, theory.fit.CHH.DMR)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.H3K27_ex_H3K9.plot <- ggplot(plotdata.DMR, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)



### *** Run models for all samples for euchromatin

### **** Run models for single cytosines

#Load the pedigree data
load(file = "./pedigree_data/euchromatin/MA_CG_euchromatin.RData") 
load(file = "./pedigree_data/euchromatin/MA_CHG_euchromatin.RData")
load(file = "./pedigree_data/euchromatin/MA_CHH_euchromatin.RData") 

#Final pedigree CG
pedigree.CG <- output.euchromatin.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.euchromatin.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.euchromatin.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.euchromatin.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.euchromatin.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.euchromatin.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.euchromatin <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

### **** Run models for DMRs


#Load data
load(file = "./pedigree_data/DMR_euchromatin/MA_DMR_euchromatin_CG.RData")
load(file = "./pedigree_data/DMR_euchromatin/MA_DMR_euchromatin_CHG.RData")
load(file = "./pedigree_data/DMR_euchromatin/MA_DMR_euchromatin_CHH.RData")


#Final pedigree CG
pedigree.CG <- output.DMR.euchromatin.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.DMR.euchromatin.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- output.DMR.euchromatin.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.DMR.euchromatin.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.DMR.euchromatin.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.DMR.euchromatin.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.DMR.euchromatin <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))



#Theoretical fits
  theory.fit.data.DMR <- neutral.model.output.CG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHG.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHG.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.model.output.CHH.DMR$for.fit.plot
  theory.fits.DMR <- c(neutral.model.output.CHH.DMR$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHH.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHH", length(theory.fits.DMR)))

 theoryFits.DMR <- rbind(theory.fit.CG.DMR, theory.fit.CHG.DMR, theory.fit.CHH.DMR)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.euchromatin.plot <- ggplot(plotdata.DMR, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)


### *** Making divergence plots

### **** WGBS, single sites subsets divergence plots

#x label
myxlab <- TeX("$\\Delta$t (mitoses)")

plotdata.combined <- data.frame(rbind(plotdata.centH3K9, plotdata.H3K9, plotdata.H3K27, plotdata.euchromatin), domain = c( rep("centromeric H3K9", nrow(plotdata.centH3K9)), rep("H3K9 excl. cent.", nrow(plotdata.centH3K9)), rep("H3K27", nrow(plotdata.centH3K9)), rep("euchromatin", nrow(plotdata.centH3K9))))

combined.plot <- ggplot(plotdata.combined, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("Single sites, WGBS") +
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_grid(fct_relevel(domain, "centromeric H3K9", "euchromatin", "H3K9 excl. cent.", "H3K27") ~ context) +
    theme(text=element_text(size=16), axis.text = element_text(size = 16))

save_plot(filename = "./epigenet_ms/fig/divsubsetssingle.pdf", combined.plot, base_height=8.5, base_width = 11)

#Making a plot
meth.centH3K9.plot <- ggplot(plotdata.centH3K9, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("Single sites, centromeric H3K9") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)


#Making a plot
meth.H3K9.plot <- ggplot(plotdata.H3K9, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("Single sites, H3K9 excl. centromeric") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)


#Making a plot
meth.H3K27.plot <- ggplot(plotdata.H3K27, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("Single sites, H3K27") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)


#Making a plot
meth.euchromatin.plot <- ggplot(plotdata.euchromatin, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("Single sites, euchromatin") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)


div.subset.plot <- plot_grid(meth.centH3K9.plot, meth.H3K9.plot, meth.H3K27.plot, meth.euchromatin.plot, nrow = 4, labels = c("A", "B", "C", "D"))


### **** WGBS, DMR subsets divergence plots

#x label
myxlab <- TeX("$\\Delta$t (mitoses)")

plotdata.DMR.combined <- data.frame(rbind(plotdata.DMR.centH3K9, plotdata.DMR.H3K9, plotdata.DMR.H3K27, plotdata.DMR.euchromatin), domain = c( rep("centromeric H3K9", nrow(plotdata.DMR.centH3K9)), rep("H3K9 excl. cent.", nrow(plotdata.DMR.centH3K9)), rep("H3K27", nrow(plotdata.DMR.centH3K9)), rep("euchromatin", nrow(plotdata.DMR.centH3K9))))

combined.DMR.plot <- ggplot(plotdata.DMR.combined, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("DMRs, WGBS") +
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_grid(fct_relevel(domain, "centromeric H3K9", "euchromatin", "H3K9 excl. cent.", "H3K27") ~ context) +
    theme(text=element_text(size=16), axis.text = element_text(size = 16))

save_plot(filename = "./epigenet_ms/fig/divsubsetsDMR.pdf", combined.DMR.plot, base_height=8.5, base_width = 11)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.centH3K9.plot <- ggplot(plotdata.DMR.centH3K9, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.H3K9_ex_cent.plot <- ggplot(plotdata.DMR.H3K9, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.H3K27_ex_H3K9.plot <- ggplot(plotdata.DMR.H3K27, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.euchromatin.plot <- ggplot(plotdata.DMR.euchromatin, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)


div.DMR.subset.plot <- plot_grid(methDMR.centH3K9.plot, methDMR.H3K9_ex_cent.plot, methDMR.H3K27_ex_H3K9.plot, methDMR.euchromatin.plot, nrow = 4, labels = c("A", "B", "C", "D"))

### *** Models for whole genome, Nanopore data

### **** Models for single sites, Nanopore data

load(file = "./pedigree_data/nanopore_ABdata/AB_output_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_CHH.RData")

#AB.output.jDMR.CG is the name of the variable where divergence data is stored

#Final pedigree CG
pedigree.CG <- AB.output.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
meth.nanopore.plot <- ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("Single sites, Nanopore") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

### **** Models for DMRs, all regions, Nanopore data


load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_CHH.RData")

#AB.output.jDMR.CG is the name of the variable where divergence data is stored

#Final pedigree CG
pedigree.CG <- AB.output.jDMR.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.jDMR.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.jDMR.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.jDMR.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.jDMR.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.jDMR.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.DMR <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methDMR.nanopore.plot <- ggplot(plotdata.DMR, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("DMRs, Nanopore") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)


### *** Models for centromeric H3K9 regions, Nanopore data

### **** For single sites

load(file = "./pedigree_data/nanopore_ABdata/AB_output_H3K9cent_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_H3K9cent_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_H3K9cent_CHH.RData")

#AB.output.jDMR.CG is the name of the variable where divergence data is stored

load(file = "./data/cytosine.prop.MA.anc.centromeric.RData") #For proportion of methylated cytosines
#cyt.props.nano.centromeric

#Final pedigree CG
pedigree.CG <- AB.output.H3K9cent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- cyt.props.nano.centromeric[2,3] #AB.output.H3K9cent.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.H3K9cent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- cyt.props.nano.centromeric[4,3] #AB.output.H3K9cent.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.H3K9cent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- cyt.props.nano.centromeric[6,3] #AB.output.H3K9cent.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Load the data from models run on the cluster
load("./modeldata/haploid/centromericH3K9/ABresults.haploid.ss.centromeric.RData")
#Variables are
#neutral.haploid.ss.centromeric.CG
#neutral.haploid.ss.centromeric.CHG
#neutral.haploid.ss.centromeric.CHH
#boot.haploid.ss.centromeric.CG
#boot.haploid.ss.centromeric.CHG
#boot.haploid.ss.centromeric.CHH

### Plotting ###
#Make plot dataframe
plotdata.H3K9cent <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

#Predicted model divergences
theoryfit.ss <- data.frame(rbind(neutral.haploid.ss.centromeric.CG$for.fit.plot, neutral.haploid.ss.centromeric.CHG$for.fit.plot, neutral.haploid.ss.centromeric.CHH$for.fit.plot), context = c(rep("CG", nrow(neutral.haploid.ss.centromeric.CG$for.fit.plot)), rep("CHG", nrow(neutral.haploid.ss.centromeric.CHG$for.fit.plot)), rep("CHH", nrow(neutral.haploid.ss.centromeric.CHH$for.fit.plot))))

### Commands for plotting are in next sections of the code

### Bootstrap estimates
#CG
boot.haploid.ss.centromeric.CG$standard.errors #Intervals
boot.haploid.ss.centromeric.CG$boot.base #Estimates
boot.haploid.ss.centromeric.CG$boot.base[,2] / boot.haploid.ss.centromeric.CG$boot.base[,1] #B/a

#CHG
boot.haploid.ss.centromeric.CHG$standard.errors #Intervals
boot.haploid.ss.centromeric.CHG$boot.base #Estimates
boot.haploid.ss.centromeric.CHG$boot.base[,2] / boot.haploid.ss.centromeric.CHG$boot.base[,1] #B/a

#CHH
boot.haploid.ss.centromeric.CHH$standard.errors #Intervals
boot.haploid.ss.centromeric.CHH$boot.base #Estimates
boot.haploid.ss.centromeric.CHH$boot.base[,2] / boot.haploid.ss.centromeric.CHH$boot.base[,1] #B/a


#################################


### 1. Do model comparisons

#Null models

outputABnull <- ABnull(pedigree.data = pedigree.final.CG, out.dir = getwd(), out.name = "./modeldata/H3K9cent_nano/ABnull_nanopore_H3K9cent_CG")

outputABnull.CHG <- ABnull(pedigree.data = pedigree.final.CHG, out.dir = getwd(), out.name = "./modeldata/H3K9cent_nano/ABnull_nanopore_H3K9cent_CHG")

outputABnull.CHH <- ABnull(pedigree.data = pedigree.final.CHH, out.dir = getwd(), out.name = "./modeldata/H3K9cent_nano/ABnull_nanopore_H3K9cent_CHH")

#Model comparisons

## CG context
comp.out.CG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_CG_centromeric_ss_nano.Rdata", pedigree.null = "./modeldata/H3K9cent_nano/ABnull_nanopore_H3K9cent_CG.Rdata") #Neutral model is preferred over the null model
comp.out.CG$Ftest
# p = 1.087136e-09 

## CHG context
comp.out.CHG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_CHG_centromeric_ss_nano.Rdata", pedigree.null = "./modeldata/H3K9cent_nano/ABnull_nanopore_H3K9cent_CHG.Rdata")
comp.out.CHG$Ftest
#Neutral model is preferred over the null model
# p = 9.108510e-09


### CHH context
comp.out.CHH <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_CHH_centromeric_ss_nano.Rdata", pedigree.null = "./modeldata/H3K9cent_nano/ABnull_nanopore_H3K9cent_CHH.Rdata") #Neutral model is preferred over the null model
comp.out.CHH$Ftest
# p = 1.647544e-09


### **** For DMRs

load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K9cent_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K9cent_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K9cent_CHH.RData")

#AB.output.jDMR.H3K9cent.CG is the name of the variable where divergence data is stored

#Final pedigree CG
pedigree.CG <- AB.output.jDMR.H3K9cent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.jDMR.H3K9cent.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.jDMR.H3K9cent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.jDMR.H3K9cent.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.jDMR.H3K9cent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.jDMR.H3K9cent.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Load the data from models run on the cluster
load("./modeldata/haploid/centromericH3K9/ABresults.haploid.DMR.centromeric.RData")
#Variables are
#neutral.haploid.DMR.centromeric.CG
#neutral.haploid.DMR.centromeric.CHG
#neutral.haploid.DMR.centromeric.CHH
#boot.haploid.DMR.centromeric.CG
#boot.haploid.DMR.centromeric.CHG
#boot.haploid.DMR.centromeric.CHH


#Plotting
#Make plot dataframe
plotdata.DMR.H3K9cent <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

#Theoretical fits
  theory.fit.data.DMR <- neutral.haploid.DMR.centromeric.CG$for.fit.plot
  theory.fits.DMR <- c(neutral.haploid.DMR.centromeric.CG$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.haploid.DMR.centromeric.CHG$for.fit.plot
  theory.fits.DMR <- c(neutral.haploid.DMR.centromeric.CHG$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHG.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHG", length(theory.fits.DMR)))

  theory.fit.data.DMR <- neutral.haploid.DMR.centromeric.CHH$for.fit.plot
  theory.fits.DMR <- c(neutral.haploid.DMR.centromeric.CHH$estimates[1, "intercept"], theory.fit.data.DMR[,"div.sim"])
  theory.fit.t.DMR <- c(0, theory.fit.data.DMR[,"delta.t"])
  theory.fit.CHH.DMR <- data.frame("divsim" = theory.fits.DMR, "dt" = theory.fit.t.DMR, context = rep("CHH", length(theory.fits.DMR)))

 theoryFits.DMR <- rbind(theory.fit.CG.DMR, theory.fit.CHG.DMR, theory.fit.CHH.DMR)

### Bootstrap estimates
#CG
boot.haploid.DMR.centromeric.CG$standard.errors #Intervals
boot.haploid.DMR.centromeric.CG$boot.base #Estimates
boot.haploid.DMR.centromeric.CG$boot.base[,2] / boot.haploid.DMR.centromeric.CG$boot.base[,1] #B/a

#CHG
boot.haploid.DMR.centromeric.CHG$standard.errors #Intervals
boot.haploid.DMR.centromeric.CHG$boot.base #Estimates
boot.haploid.DMR.centromeric.CHG$boot.base[,2] / boot.haploid.DMR.centromeric.CHG$boot.base[,1] #B/a

#CHH
boot.haploid.DMR.centromeric.CHH$standard.errors #Intervals
boot.haploid.DMR.centromeric.CHH$boot.base #Estimates
boot.haploid.DMR.centromeric.CHH$boot.base[,2] / boot.haploid.DMR.centromeric.CHH$boot.base[,1] #B/a


### Model comparisons ###
### CG ###
comp.out.CG.DMR <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_CG_centromeric_DMR_nano.Rdata", pedigree.null = "./modeldata/DMR_H3K9cent_nano/ABnull_DMR_H3K9cent_CG.Rdata")
comp.out.CG.DMR$Ftest
#Neutral model is preferred over null model
#p < 2.2 * 10^-16 

### CHG ###
comp.out.CHG.DMR <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_CHG_centromeric_DMR_nano.Rdata", pedigree.null = "./modeldata/DMR_H3K9cent_nano/ABnull_DMR_H3K9cent_CHG.Rdata")
comp.out.CHG.DMR$Ftest
#Neutral model is preferred over null model
#p = 3.248477e-07

### CHH ###
comp.out.CHH.DMR <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_CHH_centromeric_DMR_nano.Rdata", pedigree.null = "./modeldata/DMR_H3K9cent_nano/ABnull_DMR_H3K9cent_CHH.Rdata")
comp.out.CHH.DMR$Ftest
#Neutral model is preferred over null model
#p = 3.008111e-12




### *** Models for H3K9 regions excluding centromeres, nanopore

### **** For single sites

load(file = "./pedigree_data/nanopore_ABdata/AB_output_H3K9excent_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_H3K9excent_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_H3K9excent_CHH.RData")

#Final pedigree CG
pedigree.CG <- AB.output.H3K9excent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.H3K9excent.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.H3K9excent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.H3K9excent.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.H3K9excent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.H3K9excent.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.H3K9excent <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))


### **** For DMRs

load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K9excent_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K9excent_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K9excent_CHH.RData")

#AB.output.jDMR.H3K9cent.CG is the name of the variable where divergence data is stored

#Final pedigree CG
pedigree.CG <- AB.output.jDMR.H3K9excent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.jDMR.H3K9excent.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.jDMR.H3K9excent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.jDMR.H3K9excent.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.jDMR.H3K9excent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.jDMR.H3K9excent.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.DMR.H3K9excent <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

### *** Models for H3K27 regions, nanopore

### **** For single sites

load(file = "./pedigree_data/nanopore_ABdata/AB_output_H3K27_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_H3K27_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_H3K27_CHH.RData")

#Final pedigree CG
pedigree.CG <- AB.output.H3K27.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.H3K27.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.H3K27.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.H3K27.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.H3K27.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.H3K27.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.H3K27 <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

### **** For DMRs


load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K27_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K27_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K27_CHH.RData")

#Final pedigree CG
pedigree.CG <- AB.output.jDMR.H3K27.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.jDMR.H3K27.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.jDMR.H3K27.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.jDMR.H3K27.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.jDMR.H3K27.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.jDMR.H3K27.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.DMR.H3K27 <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

### *** Models for euchromatic regions, nanopore

### **** For single sites

load(file = "./pedigree_data/nanopore_ABdata/AB_output_euchromatin_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_euchromatin_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_euchromatin_CHH.RData")

#Final pedigree CG
pedigree.CG <- AB.output.euchromatin.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.euchromatin.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.euchromatin.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.euchromatin.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.euchromatin.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.euchromatin.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.euchromatin <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))


### **** For DMRs

load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_euchromatin_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_euchromatin_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_euchromatin_CHH.RData")

#Final pedigree CG
pedigree.CG <- AB.output.jDMR.euchromatin.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.jDMR.euchromatin.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.jDMR.euchromatin.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.jDMR.euchromatin.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.jDMR.euchromatin.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.jDMR.euchromatin.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.DMR.euchromatin <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))


### *** Making divergence plots for Nanopore data

### **** Nanopore, single sites subsets divergence plots

myxlab <- TeX("$\\Delta$t (mitoses)")

plotdata.combined <- data.frame(rbind(plotdata.H3K9cent, plotdata.H3K9excent, plotdata.H3K27, plotdata.euchromatin), domain = c( rep("centromeric H3K9", nrow(plotdata.H3K9cent)), rep("H3K9 excl. cent.", nrow(plotdata.H3K9cent)), rep("H3K27", nrow(plotdata.H3K9cent)), rep("euchromatin", nrow(plotdata.H3K9cent))))

combined.plot <- ggplot(plotdata.combined, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("Single sites, Nanopore") +
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_grid(fct_relevel(domain, "centromeric H3K9", "euchromatin", "H3K9 excl. cent.", "H3K27") ~ context) +
    theme(text=element_text(size=16), axis.text = element_text(size = 16))

save_plot(filename = "./epigenet_ms/fig/Nano_divsingle.pdf", combined.plot, base_height=8.5, base_width = 11)

### **** Nanopore, DMRs subsets divergence plots

myxlab <- TeX("$\\Delta$t (mitoses)")

plotdata.DMR.combined <- data.frame(rbind(plotdata.DMR.H3K9cent, plotdata.DMR.H3K9excent, plotdata.DMR.H3K27, plotdata.DMR.euchromatin), domain = c( rep("centromeric H3K9", nrow(plotdata.DMR.H3K9cent)), rep("H3K9 excl. cent.", nrow(plotdata.DMR.H3K9cent)), rep("H3K27", nrow(plotdata.DMR.H3K9cent)), rep("euchromatin", nrow(plotdata.DMR.H3K9cent))))

combined.DMR.plot <- ggplot(plotdata.DMR.combined, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("DMRs, Nanopore") +
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_grid(fct_relevel(domain, "centromeric H3K9", "euchromatin", "H3K9 excl. cent.", "H3K27") ~ context) +
    theme(text=element_text(size=16), axis.text = element_text(size = 16))

save_plot(filename = "./epigenet_ms/fig/Nano_divDMR.pdf", combined.DMR.plot, base_height=8.5, base_width = 11)

### **** Nanopore, DMRs and single sites for centromeric H3K9

myxlab <- TeX("$\\Delta$t (mitoses)")

H3K9cent.plot  <- ggplot(plotdata.H3K9cent, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    geom_line(data = theoryfit.ss, aes(x = delta.t, y = div.sim), col = "blue" ) +    
    #geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    ggtitle("Single sites, centromeric H3K9, Nanopore") +
    scale_x_continuous(limits = c(0, 800)) +
    scale_y_continuous(limits = c(0, 0.2)) +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

H3K9cent.DMR.plot <- ggplot(plotdata.DMR.H3K9cent, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryfit.DMR, aes(x = delta.t, y = div.sim), col = "blue" ) +
    geom_line(data = theoryFits.DMR, aes(x = dt, y = divsim), col = "blue" ) +
    #geom_smooth() +
    #geom_abline(intercept = 0.07175384, slope = 0.0001920564) +
    ggtitle("DMRs, centromeric H3K9, Nanopore") +
    scale_x_continuous(limits = c(0, 800)) +
    scale_y_continuous(limits = c(0, 0.2)) +
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

plot.H3K9cent.combined <- plot_grid(H3K9cent.plot, H3K9cent.DMR.plot, nrow = 2, labels = c("A", "B"))
save_plot("./epigenet_ms/fig/H3K9cent_div.pdf", plot.H3K9cent.combined, base_height = 8, base_width = 4 * 1.5 * 1.618)

### ** Observed methylation frequencies in the MA ancestors

### *** WGBS data: DMRs, whole genome
datafolder <- "~/Genomics/Neurospora/methylation/methimpute/jDMR/"

#initialize results
DMR.props.wgbs.whole <- data.frame(context = c("CG", "CG", "CHG", "CHG", "CHH", "CHH"), status = c("M", "U", "M", "U", "M", "U"), prop = 0)

matA.r1.wgbs.CG <- fread(paste0(datafolder, "MANC1_bismark_pe_trim_CG.txt"), header = T)
matA.r1.wgbs.CG.table <- summarize(group_by(matA.r1.wgbs.CG, context, status), count = n())
matA.r1.wgbs.CHG <- fread(paste0(datafolder, "MANC1_bismark_pe_trim_CHG.txt"), header = T)
matA.r1.wgbs.CHG.table <- summarize(group_by(matA.r1.wgbs.CHG, context, status), count = n())
matA.r1.wgbs.CHH <- fread(paste0(datafolder, "MANC1_bismark_pe_trim_CHH.txt"), header = T)
matA.r1.wgbs.CHH.table <- summarize(group_by(matA.r1.wgbs.CHH, context, status), count = n())

matA.r2.wgbs.CG <- fread(paste0(datafolder, "MANC1_R2_bismark_pe_trim_CG.txt"), header = T)
matA.r2.wgbs.CG.table <- summarize(group_by(matA.r2.wgbs.CG, context, status), count = n())
matA.r2.wgbs.CHG <- fread(paste0(datafolder, "MANC1_R2_bismark_pe_trim_CHG.txt"), header = T)
matA.r2.wgbs.CHG.table <- summarize(group_by(matA.r2.wgbs.CHG, context, status), count = n())
matA.r2.wgbs.CHH <- fread(paste0(datafolder, "MANC1_R2_bismark_pe_trim_CHH.txt"), header = T)
matA.r2.wgbs.CHH.table <- summarize(group_by(matA.r2.wgbs.CHH, context, status), count = n())

matA.r3.wgbs.CG <- fread(paste0(datafolder, "MANC1_R3_bismark_pe_trim_CG.txt"), header = T)
matA.r3.wgbs.CG.table <- summarize(group_by(matA.r3.wgbs.CG, context, status), count = n())
matA.r3.wgbs.CHG <- fread(paste0(datafolder, "MANC1_R3_bismark_pe_trim_CHG.txt"), header = T)
matA.r3.wgbs.CHG.table <- summarize(group_by(matA.r3.wgbs.CHG, context, status), count = n())
matA.r3.wgbs.CHH <- fread(paste0(datafolder, "MANC1_R3_bismark_pe_trim_CHH.txt"), header = T)
matA.r3.wgbs.CHH.table <- summarize(group_by(matA.r3.wgbs.CHH, context, status), count = n())

mata.r1.wgbs.CG <- fread(paste0(datafolder, "MANC2a_bismark_pe_trim_CG.txt"), header = T)
mata.r1.wgbs.CG.table <- summarize(group_by(mata.r1.wgbs.CG, context, status), count = n())
mata.r1.wgbs.CHG <- fread(paste0(datafolder, "MANC2a_bismark_pe_trim_CHG.txt"), header = T)
mata.r1.wgbs.CHG.table <- summarize(group_by(mata.r1.wgbs.CHG, context, status), count = n())
mata.r1.wgbs.CHH <- fread(paste0(datafolder, "MANC2a_bismark_pe_trim_CHH.txt"), header = T)
mata.r1.wgbs.CHH.table <- summarize(group_by(mata.r1.wgbs.CHH, context, status), count = n())

mata.r2.wgbs.CG <- fread(paste0(datafolder, "MANC2a_R2_bismark_pe_trim_CG.txt"), header = T)
mata.r2.wgbs.CG.table <- summarize(group_by(mata.r2.wgbs.CG, context, status), count = n())
mata.r2.wgbs.CHG <- fread(paste0(datafolder, "MANC2a_R2_bismark_pe_trim_CHG.txt"), header = T)
mata.r2.wgbs.CHG.table <- summarize(group_by(mata.r2.wgbs.CHG, context, status), count = n())
mata.r2.wgbs.CHH <- fread(paste0(datafolder, "MANC2a_R2_bismark_pe_trim_CHH.txt"), header = T)
mata.r2.wgbs.CHH.table <- summarize(group_by(mata.r2.wgbs.CHH, context, status), count = n())

mata.r3.wgbs.CG <- fread(paste0(datafolder, "MANC2a_R3_bismark_pe_trim_CG.txt"), header = T)
mata.r3.wgbs.CG.table <- summarize(group_by(mata.r3.wgbs.CG, context, status), count = n())
mata.r3.wgbs.CHG <- fread(paste0(datafolder, "MANC2a_R3_bismark_pe_trim_CHG.txt"), header = T)
mata.r3.wgbs.CHG.table <- summarize(group_by(mata.r3.wgbs.CHG, context, status), count = n())
mata.r3.wgbs.CHH <- fread(paste0(datafolder, "MANC2a_R3_bismark_pe_trim_CHH.txt"), header = T)
mata.r3.wgbs.CHH.table <- summarize(group_by(mata.r3.wgbs.CHH, context, status), count = n())

CG.table <- matA.r1.wgbs.CG.table[,3] + matA.r2.wgbs.CG.table[,3] + matA.r3.wgbs.CG.table[,3] + mata.r1.wgbs.CG.table[,3] + mata.r2.wgbs.CG.table[,3] + mata.r3.wgbs.CG.table[,3]
DMR.props.wgbs.whole[1:2,3] <- c(CG.table[1,1] / sum(CG.table), CG.table[2,1] / sum(CG.table))

CHG.table <- matA.r1.wgbs.CHG.table[,3] + matA.r2.wgbs.CHG.table[,3] + matA.r3.wgbs.CHG.table[,3] + mata.r1.wgbs.CHG.table[,3] + mata.r2.wgbs.CHG.table[,3] + mata.r3.wgbs.CHG.table[,3]
DMR.props.wgbs.whole[3:4,3] <- c(CHG.table[1,1] / sum(CHG.table), CHG.table[2,1] / sum(CHG.table))

CHH.table <- matA.r1.wgbs.CHH.table[,3] + matA.r2.wgbs.CHH.table[,3] + matA.r3.wgbs.CHH.table[,3] + mata.r1.wgbs.CHH.table[,3] + mata.r2.wgbs.CHH.table[,3] + mata.r3.wgbs.CHH.table[,3]
DMR.props.wgbs.whole[5:6,3] <- c(CHH.table[1,1] / sum(CHH.table), CHH.table[2,1] / sum(CHH.table))

save(DMR.props.wgbs.whole, file = "./data/DMR.prop.MA.anc.whole.RData")

### *** WGBS data: single sites, whole genome

datafolder <- "~/Genomics/Neurospora/methylation/methimpute/all/"

#mat A ancestor replicates
matA.r1.wgbs <- fread(paste0(datafolder, "methylome_MANC1_bismark_pe_trim_All.txt"), header = T)
matA.r1.wgbs <- filter(matA.r1.wgbs, posteriorMax >= 0.99)
matA.r1.wgbs.table <- summarize(group_by(matA.r1.wgbs, context, status), count = n())
rm(matA.r1.wgbs) #Remove to free up memory

matA.r2.wgbs <- fread(paste0(datafolder, "methylome_MANC1_R2_bismark_pe_trim_All.txt"), header = T)
matA.r2.wgbs <- filter(matA.r2.wgbs, posteriorMax >= 0.99)
matA.r2.wgbs.table <- summarize(group_by(matA.r2.wgbs, context, status), count = n())
rm(matA.r2.wgbs) #Remove to free up memory

matA.r3.wgbs <- fread(paste0(datafolder, "methylome_MANC1_R3_bismark_pe_trim_All.txt"), header = T)
matA.r3.wgbs <- filter(matA.r3.wgbs, posteriorMax >= 0.99)
matA.r3.wgbs.table <- summarize(group_by(matA.r3.wgbs, context, status), count = n())
rm(matA.r3.wgbs) #Remove to free up memory

#mat a ancestor replicates
mata.r1.wgbs <- fread(paste0(datafolder, "methylome_MANC2a_bismark_pe_trim_All.txt"), header = T)
mata.r1.wgbs <- filter(mata.r1.wgbs, posteriorMax >= 0.99)
mata.r1.wgbs.table <- summarize(group_by(mata.r1.wgbs, context, status), count = n())
rm(mata.r1.wgbs) #Remove to free up memory

mata.r2.wgbs <- fread(paste0(datafolder, "methylome_MANC2a_R2_bismark_pe_trim_All.txt"), header = T)
mata.r2.wgbs <- filter(mata.r2.wgbs, posteriorMax >= 0.99)
mata.r2.wgbs.table <- summarize(group_by(mata.r2.wgbs, context, status), count = n())
rm(mata.r2.wgbs) #Remove to free up memory

mata.r3.wgbs <- fread(paste0(datafolder, "methylome_MANC2a_R3_bismark_pe_trim_All.txt"), header = T)
mata.r3.wgbs <- filter(mata.r3.wgbs, posteriorMax >= 0.99)
mata.r3.wgbs.table <- summarize(group_by(mata.r3.wgbs, context, status), count = n())
rm(mata.r3.wgbs) #Remove to free up memory

#Sum all replicates together
cyt.wgbs.whole <- matA.r1.wgbs.table[,1:2]
cyt.wgbs.whole[,3] <- matA.r1.wgbs.table[,3] + matA.r2.wgbs.table[,3] + matA.r3.wgbs.table[,3] + mata.r1.wgbs.table[,3] + mata.r2.wgbs.table[,3] + mata.r3.wgbs.table[,3]

cyt.props.wgbs.whole <- cyt.wgbs.whole[-c(1,4,7),1:2]
cyt.props.wgbs.whole$prop <- 0
cg.props <- unlist(c( (cyt.wgbs.whole[1,3] + cyt.wgbs.whole[2,3]) / (sum(cyt.wgbs.whole[1:3,3])), (cyt.wgbs.whole[3,3]) / (sum(cyt.wgbs.whole[1:3,3])) ))
chg.props <- unlist(c( (cyt.wgbs.whole[4,3] + cyt.wgbs.whole[5,3]) / (sum(cyt.wgbs.whole[4:6,3])), (cyt.wgbs.whole[6,3]) / (sum(cyt.wgbs.whole[4:6,3])) ))
chh.props <- unlist(c( (cyt.wgbs.whole[7,3] + cyt.wgbs.whole[8,3]) / (sum(cyt.wgbs.whole[7:9,3])), (cyt.wgbs.whole[9,3]) / (sum(cyt.wgbs.whole[7:9,3])) ))
cyt.props.wgbs.whole$prop <- c(cg.props, chg.props, chh.props)

cyt.props.wgbs.whole <- data.frame(cyt.props.wgbs.whole)

save(cyt.props.wgbs.whole, file = "./data/cytosine.prop.MA.anc.whole.RData")

#Load the data
datafolder <- "~/Genomics/Neurospora/methylation/methimpute/subsets/"

### *** WGBS data: single sites, centromeric regions
datafolder <- "~/Genomics/Neurospora/methylation/methimpute/subsets/"

anc.matA.r1.wgbs <- read.table(paste0(datafolder, "MANC1_bismark_pe_trim_All_centromericH3K9_methimpute.txt"), header = T)
all(anc.matA.r1.wgbs$posteriorMax >= 0.99)

anc.matA.r1.wgbs.cent <- summarize(group_by(anc.matA.r1.wgbs, context, status), count = n())

anc.matA.r2.wgbs <- read.table(paste0(datafolder, "MANC1_R2_bismark_pe_trim_All_centromericH3K9_methimpute.txt"), header = T)

anc.matA.r2.wgbs.cent <- summarize(group_by(anc.matA.r2.wgbs, context, status), count = n())

anc.matA.r3.wgbs <- read.table(paste0(datafolder, "MANC1_R3_bismark_pe_trim_All_centromericH3K9_methimpute.txt"), header = T)

anc.matA.r3.wgbs.cent <- summarize(group_by(anc.matA.r3.wgbs, context, status), count = n())

anc.mata.r1.wgbs <- read.table(paste0(datafolder, "MANC2a_bismark_pe_trim_All_centromericH3K9_methimpute.txt"), header = T)
all(anc.mata.r1.wgbs$posteriorMax >= 0.99)

anc.mata.r1.wgbs.cent <- summarize(group_by(anc.mata.r1.wgbs, context, status), count = n())

anc.mata.r2.wgbs <- read.table(paste0(datafolder, "MANC2a_R2_bismark_pe_trim_All_centromericH3K9_methimpute.txt"), header = T)

anc.mata.r2.wgbs.cent <- summarize(group_by(anc.mata.r2.wgbs, context, status), count = n())

anc.mata.r3.wgbs <- read.table(paste0(datafolder, "MANC2a_R3_bismark_pe_trim_All_centromericH3K9_methimpute.txt"), header = T)

anc.mata.r3.wgbs.cent <- summarize(group_by(anc.mata.r3.wgbs, context, status), count = n())

#Sum all replicates together
cyt.wgbs.cent <- anc.mata.r1.wgbs.cent[,1:2]
cyt.wgbs.cent[,3] <- anc.matA.r1.wgbs.cent[,3] + anc.matA.r2.wgbs.cent[,3] + anc.matA.r3.wgbs.cent[,3] + anc.mata.r1.wgbs.cent[,3] +  anc.mata.r2.wgbs.cent[,3] + anc.mata.r3.wgbs.cent[,3]

cyt.props.wgbs.centromeric <- cyt.wgbs.cent[-c(1,4,7),1:2]
cyt.props.wgbs.centromeric$prop <- 0
cg.props <- unlist(c( (cyt.wgbs.cent[1,3] + cyt.wgbs.cent[2,3]) / (sum(cyt.wgbs.cent[1:3,3])), (cyt.wgbs.cent[3,3]) / (sum(cyt.wgbs.cent[1:3,3])) ))
chg.props <- unlist(c( (cyt.wgbs.cent[4,3] + cyt.wgbs.cent[5,3]) / (sum(cyt.wgbs.cent[4:6,3])), (cyt.wgbs.cent[6,3]) / (sum(cyt.wgbs.cent[4:6,3])) ))
chh.props <- unlist(c( (cyt.wgbs.cent[7,3] + cyt.wgbs.cent[8,3]) / (sum(cyt.wgbs.cent[7:9,3])), (cyt.wgbs.cent[9,3]) / (sum(cyt.wgbs.cent[7:9,3])) ))
cyt.props.wgbs.centromeric$prop <- c(cg.props, chg.props, chh.props)

cyt.props.wgbs.centromeric <- data.frame(cyt.props.wgbs.centromeric)

### *** WGBS data: DMRs, centromeric regions

datafolder <- "~/Genomics/Neurospora/methylation/methimpute/jDMR/"

#initialize results
DMR.props.wgbs.cent <- data.frame(context = c("CG", "CG", "CHG", "CHG", "CHH", "CHH"), status = c("M", "U", "M", "U", "M", "U"), prop = 0)

matA.r1.wgbs.CG <- fread(paste0(datafolder, "MANC1_bismark_pe_trim_All_centromericH3K9_methimpute_CG.txt"), header = T)
matA.r1.wgbs.CG.table <- summarize(group_by(matA.r1.wgbs.CG, context, status), count = n())
matA.r1.wgbs.CHG <- fread(paste0(datafolder, "MANC1_bismark_pe_trim_All_centromericH3K9_methimpute_CHG.txt"), header = T)
matA.r1.wgbs.CHG.table <- summarize(group_by(matA.r1.wgbs.CHG, context, status), count = n())
matA.r1.wgbs.CHH <- fread(paste0(datafolder, "MANC1_bismark_pe_trim_All_centromericH3K9_methimpute_CHH.txt"), header = T)
matA.r1.wgbs.CHH.table <- summarize(group_by(matA.r1.wgbs.CHH, context, status), count = n())

matA.r2.wgbs.CG <- fread(paste0(datafolder, "MANC1_R2_bismark_pe_trim_All_centromericH3K9_methimpute_CG.txt"), header = T)
matA.r2.wgbs.CG.table <- summarize(group_by(matA.r2.wgbs.CG, context, status), count = n())
matA.r2.wgbs.CHG <- fread(paste0(datafolder, "MANC1_R2_bismark_pe_trim_All_centromericH3K9_methimpute_CHG.txt"), header = T)
matA.r2.wgbs.CHG.table <- summarize(group_by(matA.r2.wgbs.CHG, context, status), count = n())
matA.r2.wgbs.CHH <- fread(paste0(datafolder, "MANC1_R2_bismark_pe_trim_All_centromericH3K9_methimpute_CHH.txt"), header = T)
matA.r2.wgbs.CHH.table <- summarize(group_by(matA.r2.wgbs.CHH, context, status), count = n())

matA.r3.wgbs.CG <- fread(paste0(datafolder, "MANC1_R3_bismark_pe_trim_All_centromericH3K9_methimpute_CG.txt"), header = T)
matA.r3.wgbs.CG.table <- summarize(group_by(matA.r3.wgbs.CG, context, status), count = n())
matA.r3.wgbs.CHG <- fread(paste0(datafolder, "MANC1_R3_bismark_pe_trim_All_centromericH3K9_methimpute_CHG.txt"), header = T)
matA.r3.wgbs.CHG.table <- summarize(group_by(matA.r3.wgbs.CHG, context, status), count = n())
matA.r3.wgbs.CHH <- fread(paste0(datafolder, "MANC1_R3_bismark_pe_trim_All_centromericH3K9_methimpute_CHH.txt"), header = T)
matA.r3.wgbs.CHH.table <- summarize(group_by(matA.r3.wgbs.CHH, context, status), count = n())

mata.r1.wgbs.CG <- fread(paste0(datafolder, "MANC2a_bismark_pe_trim_All_centromericH3K9_methimpute_CG.txt"), header = T)
mata.r1.wgbs.CG.table <- summarize(group_by(mata.r1.wgbs.CG, context, status), count = n())
mata.r1.wgbs.CHG <- fread(paste0(datafolder, "MANC2a_bismark_pe_trim_All_centromericH3K9_methimpute_CHG.txt"), header = T)
mata.r1.wgbs.CHG.table <- summarize(group_by(mata.r1.wgbs.CHG, context, status), count = n())
mata.r1.wgbs.CHH <- fread(paste0(datafolder, "MANC2a_bismark_pe_trim_All_centromericH3K9_methimpute_CHH.txt"), header = T)
mata.r1.wgbs.CHH.table <- summarize(group_by(mata.r1.wgbs.CHH, context, status), count = n())

mata.r2.wgbs.CG <- fread(paste0(datafolder, "MANC2a_R2_bismark_pe_trim_All_centromericH3K9_methimpute_CG.txt"), header = T)
mata.r2.wgbs.CG.table <- summarize(group_by(mata.r2.wgbs.CG, context, status), count = n())
mata.r2.wgbs.CHG <- fread(paste0(datafolder, "MANC2a_R2_bismark_pe_trim_All_centromericH3K9_methimpute_CHG.txt"), header = T)
mata.r2.wgbs.CHG.table <- summarize(group_by(mata.r2.wgbs.CHG, context, status), count = n())
mata.r2.wgbs.CHH <- fread(paste0(datafolder, "MANC2a_R2_bismark_pe_trim_All_centromericH3K9_methimpute_CHH.txt"), header = T)
mata.r2.wgbs.CHH.table <- summarize(group_by(mata.r2.wgbs.CHH, context, status), count = n())

mata.r3.wgbs.CG <- fread(paste0(datafolder, "MANC2a_R3_bismark_pe_trim_All_centromericH3K9_methimpute_CG.txt"), header = T)
mata.r3.wgbs.CG.table <- summarize(group_by(mata.r3.wgbs.CG, context, status), count = n())
mata.r3.wgbs.CHG <- fread(paste0(datafolder, "MANC2a_R3_bismark_pe_trim_All_centromericH3K9_methimpute_CHG.txt"), header = T)
mata.r3.wgbs.CHG.table <- summarize(group_by(mata.r3.wgbs.CHG, context, status), count = n())
mata.r3.wgbs.CHH <- fread(paste0(datafolder, "MANC2a_R3_bismark_pe_trim_All_centromericH3K9_methimpute_CHH.txt"), header = T)
mata.r3.wgbs.CHH.table <- summarize(group_by(mata.r3.wgbs.CHH, context, status), count = n())

CG.table <- matA.r1.wgbs.CG.table[,3] + matA.r2.wgbs.CG.table[,3] + matA.r3.wgbs.CG.table[,3] + mata.r1.wgbs.CG.table[,3] + mata.r2.wgbs.CG.table[,3] + mata.r3.wgbs.CG.table[,3]
DMR.props.wgbs.cent[1:2,3] <- c(CG.table[1,1] / sum(CG.table), CG.table[2,1] / sum(CG.table))

CHG.table <- matA.r1.wgbs.CHG.table[,3] + matA.r2.wgbs.CHG.table[,3] + matA.r3.wgbs.CHG.table[,3] + mata.r1.wgbs.CHG.table[,3] + mata.r2.wgbs.CHG.table[,3] + mata.r3.wgbs.CHG.table[,3]
DMR.props.wgbs.cent[3:4,3] <- c(CHG.table[1,1] / sum(CHG.table), CHG.table[2,1] / sum(CHG.table))

CHH.table <- matA.r1.wgbs.CHH.table[,3] + matA.r2.wgbs.CHH.table[,3] + matA.r3.wgbs.CHH.table[,3] + mata.r1.wgbs.CHH.table[,3] + mata.r2.wgbs.CHH.table[,3] + mata.r3.wgbs.CHH.table[,3]
DMR.props.wgbs.cent[5:6,3] <- c(CHH.table[1,1] / sum(CHH.table), CHH.table[2,1] / sum(CHH.table))



### *** Nanopore data: single sites, centromeric regions
anc.matA.nano <- read.table(paste0(datafolder, "2489_MatA_All_centromericH3K9_methimpute.txt"), header = T)

all(anc.matA.nano$posteriorMax >= 0.99) #TRUE -> this data contains only high quality sites

anc.matA.nano.cent <- summarize(group_by(anc.matA.nano, context, status), count = n())

anc.mata.nano <- read.table(paste0(datafolder, "2489_Mata_All_centromericH3K9_methimpute.txt"), header = T)

#all(anc.mata.nano$posteriorMax >= 0.99) #TRUE

anc.mata.nano.cent <- summarize(group_by(anc.mata.nano, context, status), count = n())

#Calculating the average methylation frequencies for the two MA ancestors

cyt.props.nano.centromeric <- anc.matA.nano.cent[,1:2]

cyt.props.nano.centromeric$prop <- unlist( c( (anc.matA.nano.cent[1,3] + anc.mata.nano.cent[1,3]) / (anc.matA.nano.cent[1,3] + anc.mata.nano.cent[1,3] + anc.matA.nano.cent[2,3] + anc.mata.nano.cent[2,3]), (anc.matA.nano.cent[2,3] + anc.mata.nano.cent[2,3]) / (anc.matA.nano.cent[1,3] + anc.mata.nano.cent[1,3] + anc.matA.nano.cent[2,3] + anc.mata.nano.cent[2,3]), (anc.matA.nano.cent[3,3] + anc.mata.nano.cent[3,3]) / (anc.matA.nano.cent[3,3] + anc.mata.nano.cent[3,3] + anc.matA.nano.cent[4,3] + anc.mata.nano.cent[4,3]), (anc.matA.nano.cent[4,3] + anc.mata.nano.cent[4,3]) / (anc.matA.nano.cent[3,3] + anc.mata.nano.cent[3,3] + anc.matA.nano.cent[4,3] + anc.mata.nano.cent[4,3]), (anc.matA.nano.cent[5,3] + anc.mata.nano.cent[5,3]) / (anc.matA.nano.cent[5,3] + anc.mata.nano.cent[5,3] + anc.matA.nano.cent[6,3] + anc.mata.nano.cent[6,3]), (anc.matA.nano.cent[6,3] + anc.mata.nano.cent[6,3]) / (anc.matA.nano.cent[5,3] + anc.mata.nano.cent[5,3] + anc.matA.nano.cent[6,3] + anc.mata.nano.cent[6,3]) ))

cyt.props.nano.centromeric <- data.frame(cyt.props.nano.centromeric)
#This is methylation frequencies (mean of the two MA ancestors) for nanopore data, centrometic regions

save(cyt.props.wgbs.centromeric, cyt.props.nano.centromeric, file = "./data/cytosine.prop.MA.anc.centromeric.RData")

########################################################################################

### *** Nanopore data: DMRs, centromeric regions

datafolder <- "~/Genomics/Neurospora/methylation/methimpute/nanopore_centromeric/"

DMR.props.nano.cent <- data.frame(context = c("CG", "CG", "CHG", "CHG", "CHH", "CHH"), status = c("M", "U", "M", "U", "M", "U"), prop = 0)

anc.matA.nano.CG <- fread(paste0(datafolder, "2489_MatA_All_centromericH3K9_methimpute_CG.txt"), header = T)
anc.matA.nano.CG.table <- summarize(group_by(anc.matA.nano.CG, context, status), count = n())
anc.matA.nano.CHG <- fread(paste0(datafolder, "2489_MatA_All_centromericH3K9_methimpute_CHG.txt"), header = T)
anc.matA.nano.CHG.table <- summarize(group_by(anc.matA.nano.CHG, context, status), count = n())
anc.matA.nano.CHH <- fread(paste0(datafolder, "2489_MatA_All_centromericH3K9_methimpute_CHH.txt"), header = T)
anc.matA.nano.CHH.table <- summarize(group_by(anc.matA.nano.CHH, context, status), count = n())

anc.mata.nano.CG <- fread(paste0(datafolder, "2489_Mata_All_centromericH3K9_methimpute_CG.txt"), header = T)
anc.mata.nano.CG.table <- summarize(group_by(anc.mata.nano.CG, context, status), count = n())
anc.mata.nano.CHG <- fread(paste0(datafolder, "2489_Mata_All_centromericH3K9_methimpute_CHG.txt"), header = T)
anc.mata.nano.CHG.table <- summarize(group_by(anc.mata.nano.CHG, context, status), count = n())
anc.mata.nano.CHH <- fread(paste0(datafolder, "2489_Mata_All_centromericH3K9_methimpute_CHH.txt"), header = T)
anc.mata.nano.CHH.table <- summarize(group_by(anc.mata.nano.CHH, context, status), count = n())

CG.table <- anc.matA.nano.CG.table[,3] + anc.mata.nano.CG.table[,3]
DMR.props.nano.cent[1:2,3] <- c(CG.table[1,1] / sum(CG.table), CG.table[2,1] / sum(CG.table))

CHG.table <- anc.matA.nano.CHG.table[,3] + anc.mata.nano.CHG.table[,3]
DMR.props.nano.cent[3:4,3] <- c(CHG.table[1,1] / sum(CHG.table), CHG.table[2,1] / sum(CHG.table))

CHH.table <- anc.matA.nano.CHH.table[,3] + anc.mata.nano.CHH.table[,3]
DMR.props.nano.cent[5:6,3] <- c(CHH.table[1,1] / sum(CHH.table), CHH.table[2,1] / sum(CHH.table))

save(DMR.props.wgbs.cent, DMR.props.nano.cent, file = "./data/DMR.prop.MA.anc.cent.RData")


### * Exploring DMRs

### ** Analysis of DMRs

### *** Across the whole genome

### **** Load the data, combine different contexts, and annotate

DMRs.all.CG <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all/CG_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CG$context <- "CG"

DMRs.all.CHG <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all/CHG_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHG$context <- "CHG"

DMRs.all.CHH <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all/CHH_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHH$context <- "CHH"

DMRs.all <- rbind(DMRs.all.CG, DMRs.all.CHG, DMRs.all.CHH)
DMRs.all <- arrange(DMRs.all, seqnames, start)

#There are some bins that have been called in multiple contexts
bincheck <- paste0(DMRs.all[,1],DMRs.all[,2])
duplicatebins <- duplicated(bincheck)
#DMRs.all[duplicatebins,c(1:3,74)]

#Filter bins that are in different contexts (that is get only unique locations)
DMRs.all <- DMRs.all[!duplicatebins,]

### After combining different contexts there are some DMRs that should be combined into one (i.e. they are next to each other in the coordinates, e.g. CG and CHG DMRs that are next to each other) #TO DO!!!

#Filter everything that is not in the seven chromosomes
DMRs.all <- filter(DMRs.all, seqnames %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Load annotations
genes <- read.gff("~/Genomics/Neurospora/reference/Neurospora_crassa_NC12_genes.gff3")
genes <- filter(genes, type == "gene")
#We are only considering the seven chromosomes
genes <- filter(genes, seqid %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7") == T)
genes$seqid <- factor(genes$seqid) #Drop levels that are ot used

#Load TE annotations
TEs <- read.gff("~/Genomics/Neurospora/reference/Neurospora_crassa_NC12_TE.gff3")

#Promoter annotations, promdist the distance from TSS that still counts as promoter region
promoters <- annotate.promoters(genes, promdist = 500)

#This annotates all DMRs
DMRs.all <- annotate.DMRs(DMRs.all, genes, TEs, promoters) #Some nearby DMRs can occur in the same gene

## There are some warnings that occur
## Use options(warn=2), to convert warnings into errors and check what is the problem
## options(error = recover)
##

#Then we need information about the chromatin domains
#Load the data
#Euchromatin
euchr <- read.table("./annotation/2489.euchromatin.bed", header = F, sep = "\t")
colnames(euchr) <- c("Chromosome", "start", "end")
euchr <- filter(euchr, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Centromeric H3K9
centromeric.h3k9 <- read.table("./annotation/centromericH3K9.bed", header = F, sep = "\t")
colnames(centromeric.h3k9)[1:3] <- c("Chromosome", "start", "end")

#H3K9 not in centromeric regions
excent2.h3k9 <- read.table("./annotation/H3K9_ex_cent2.bed", header = F, sep = "\t")
colnames(excent2.h3k9)[1:3] <- c("Chromosome", "start", "end")
excent2.h3k9 <- filter(excent2.h3k9, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#H3K27 not in H3K9 regions
h3k27.exh3k9 <- read.table("./annotation/H3K27_exH3K9.bed", header = F, sep = "\t")
colnames(h3k27.exh3k9)[1:3] <- c("Chromosome", "start", "end")
h3k27.exh3k9 <- filter(h3k27.exh3k9, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Check DMR overlap for chromatin domains
DMRs.all <- DMR.by.domains(DMRs.all, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)

##Then need to check overlap for each gene and promoter in the different domains
genes <- annotation.by.domains(genes, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)
TEs <- annotation.by.domains(TEs, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)
promoters <- annotation.by.domains(promoters, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9, prom = TRUE)

#Then check whether methylation was gained or lost
#Load the methylation level files
DMRs.all.CG.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all/CG_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CG.meth$context <- "CG"

DMRs.all.CHG.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all/CHG_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHG.meth$context <- "CHG"

DMRs.all.CHH.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all/CHH_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHH.meth$context <- "CHH"

#For each DMR check whether methylation was gained or lost
meth <- classify.DMRs(DMRs.all, DMRs.all.CG.meth, DMRs.all.CHG.meth, DMRs.all.CHH.meth)
DMRs.all <- cbind(DMRs.all, meth) #Store methylation levels, and type

#Saving the data for easier loading the next time
save(DMRs.all, euchr, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9, genes, TEs, promoters, file = "./data/DMRs_all.RData")

### **** Analysis of all DMRs that occured across the MA experiment

load("./data/DMRs_all.RData") #Load the DMR data, domain, and annotations
colnames(DMRs.all)[1] <- "Chromosome" #Change seqname to chromosome

nrow(DMRs.all) #24208 DMRs observed in total

summary(DMRs.all$domain) #Breakdown of DMRs by domain

DMRs.all$anno.type <- factor(DMRs.all$anno.type)
summary(DMRs.all$anno.type) #Breakdown of DMRs by annotation

#Count DMRs for different domains and annotations
DMR.counts <- summarise(group_by(DMRs.all, anno.type, domain), count = n())

DMR.domains <- summarise(group_by(DMRs.all, domain), count = n())
#Then need the lenghts of these domains
domain.lengths <- c( sum(centromeric.h3k9$end - centromeric.h3k9$start) , sum(euchr$end - euchr$start), sum(h3k27.exh3k9$end - h3k27.exh3k9$start), sum(excent2.h3k9$end - excent2.h3k9$start) )
DMR.domains$sequence <- domain.lengths

#There are some genes in the centromeric regions, what are those?
filter(DMRs.all, anno.type == "gene", domain == "Centromeric")

#################################################################################3
##Poisson models for Annotation types considering different domains (that is: is it just that that TE's have so many DMRs so that is why centromeric regions have the most DMRs? Or do centromeric regions (and H3K9) have more DMRs than could be expected based on their TE content? ##T######

## Poisson model for DMR counts in the different domains...
#Reorder levels, so that euchromatin the one where everything else is compated to
DMR.domains$domain <- factor(DMR.domains$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27") )
#Expected counts, if DMRs would occur randomly across the genome
DMR.domains$expected <- sum(DMR.domains$count) * ( DMR.domains$sequence / sum(DMR.domains$sequence))

#Poisson model
#malli1 <- glm(count ~ offset(log(sequence)) + domain, family = poisson, data = DMR.domains)
malli1 <- glm(count ~ -1 + offset(log(expected)) + domain, family = poisson, data = DMR.domains)

#Format results
malli1.results <- cbind(coef(malli1), confint(malli1)) #Extract coefficients and conf int
malli1.results <- exp(malli1.results) #Change back to normal scale
malli1.results <- data.frame(malli1.results) #Change to dataframe
malli1.results$domain <- factor(c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
malli1.results$domain <- factor(malli1.results$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
colnames(malli1.results)[1:3] <- c("estimate", "lower", "upper")

##Then checking by domain and annotation
##Calculate lengths by annotation and domain
DMR.counts$sequence <- NA
DMR.counts[1,4] <- sum( filter(TEs, domain == "Centromeric")$end - filter(TEs, domain == "Centromeric")$start)
DMR.counts[2,4] <- sum( filter(TEs, domain == "Euchromatin")$end - filter(TEs, domain == "Euchromatin")$start)
DMR.counts[3,4] <- sum( filter(TEs, domain == "H3K27")$end - filter(TEs, domain == "H3K27")$start)
DMR.counts[4,4] <- sum( filter(TEs, domain == "H3K9")$end - filter(TEs, domain == "H3K9")$start)
#
DMR.counts[5,4] <- sum( filter(genes, domain == "Centromeric")$end - filter(genes, domain == "Centromeric")$start)
DMR.counts[6,4] <- sum( filter(genes, domain == "Euchromatin")$end - filter(genes, domain == "Euchromatin")$start)
DMR.counts[7,4] <- sum( filter(genes, domain == "H3K27")$end - filter(genes, domain == "H3K27")$start)
DMR.counts[8,4] <- sum( filter(genes, domain == "H3K9")$end - filter(genes, domain == "H3K9")$start)
#
DMR.counts[13,4] <- sum( filter(promoters, domain == "Centromeric")$p.end - filter(promoters, domain == "Centromeric")$p.start)
DMR.counts[14,4] <- sum( filter(promoters, domain == "Euchromatin")$p.end - filter(promoters, domain == "Euchromatin")$p.start)
DMR.counts[15,4] <- sum( filter(promoters, domain == "H3K27")$p.end - filter(promoters, domain == "H3K27")$p.start)
DMR.counts[16,4] <- sum( filter(promoters, domain == "H3K9")$p.end - filter(promoters, domain == "H3K9")$p.start)
#
DMR.counts[9,4] <- sum(centromeric.h3k9$end - centromeric.h3k9$start) - DMR.counts[1,4] - DMR.counts[5,4] - DMR.counts[13,4]
DMR.counts[10,4] <- sum(euchr$end - euchr$start) - DMR.counts[2,4] - DMR.counts[6,4] - DMR.counts[14,4]
DMR.counts[11,4] <- sum(h3k27.exh3k9$end - h3k27.exh3k9$start) - DMR.counts[3,4] - DMR.counts[7,4] - DMR.counts[15,4]
DMR.counts[12,4] <- sum(excent2.h3k9$end - excent2.h3k9$start) - DMR.counts[4,4] - DMR.counts[8,4] - DMR.counts[16,4]

#Releveling so that Euchromatin and intergenic are set as base comparisons
DMR.counts$domain <- factor(DMR.counts$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27") )
DMR.counts$anno.type <- factor(DMR.counts$anno.type, levels = c("intergenic", "TE", "gene", "promoter"))
#Expected counts, if DMRs would occur randomly across the genome
DMR.counts$expected <- sum(DMR.counts$count) * ( DMR.counts$sequence / sum(DMR.counts$sequence))

#Calculating the percentages
DMR.counts$percent <- rep(0,nrow(DMR.counts))
DMR.counts[c(1,5,9,13),]$percent <- round(DMR.counts[c(1,5,9,13),]$count / sum(DMR.counts[c(1,5,9,13),]$count) * 100,1) #Centromeric
DMR.counts[c(2,6,10,14),]$percent <- round(DMR.counts[c(2,6,10,14),]$count / sum(DMR.counts[c(2,6,10,14),]$count) * 100,1) #Euchromatin
DMR.counts[c(3,7,11,15),]$percent <- round(DMR.counts[c(3,7,11,15),]$count / sum(DMR.counts[c(3,7,11,15),]$count) * 100,1) #H3K27
DMR.counts[c(4,8,12,16),]$percent <- round(DMR.counts[c(4,8,12,16),]$count / sum(DMR.counts[c(4,8,12,16),]$count) * 100,1) #H3K9


#Poisson model
#malli2 <- glm(count ~ offset(log(sequence)) + domain + anno.type, family = poisson, data = DMR.counts)
#TE's do not seem to be especially enriched

#malli2 <- glm(count ~ offset(log(sequence)) + domain + anno.type + domain:anno.type, family = poisson, data = DMR.counts)

#malli2 <- glm(count ~ -1 + offset(log(expected)) + domain + anno.type + domain:anno.type, family = poisson, data = DMR.counts)

#Marginal effect of TEs
#margins(malli2, variable = 'anno.type')

#poissonmfx(count ~ offset(log(expected)) + domain + anno.type +  domain:anno.type, data = DMR.counts)

#m2 <- brm(count ~ offset(log(expected)) + domain + anno.type +  domain:anno.type, family = poisson, data = DMR.counts,
#    warmup = 1000, iter = 4000, chains = 4, cores = 4)
#me <- marginal_effects(m2, variable = 'anno.type')

malli3 <- glm(count ~ -1 + offset(log(expected)) +  domain:anno.type, family = poisson, data = DMR.counts)


#coef(malli3)

#Format results
malli3.results <- cbind(coef(malli3), confint(malli3)) #Extract coefficients and conf int
malli3.results <- exp(malli3.results) #Change back to normal scale
malli3.results <- data.frame(malli3.results) #Change to data frame

malli3.results$domain <- factor(rep(c("Euchromatin", "Centromeric", "H3K9", "H3K27"),4))
malli3.results$domain <- factor(malli3.results$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
malli3.results$annotation <- factor(c(rep("intergenic",4), rep("TE", 4), rep("gene", 4), rep("promoter", 4)))
colnames(malli3.results)[1:3] <- c("estimate", "lower", "upper")


#malli3.results.domains <- malli3.results[1:3,]
#malli3.results.domains$domain <- factor(c("Centromeric", "H3K9", "H3K27"))
#colnames(malli3.results.domains)[1:3] <- c("estimate", "lower", "upper")
#
#malli3.results.annotation <- malli3.results[-c(1:3),]
#malli3.results.annotation$domain <- factor(c("Euchromatin", "Euchromatin", "Euchromatin", rep(c("Centromeric", "H3K9", "H3K27"),3)))
#malli3.results.annotation$annotation <- factor(c("TE", "gene", "promoter", rep("TE", 3), rep("gene", 3), rep("promoter", 3)))
#colnames(malli3.results.annotation)[1:3] <- c("estimate", "lower", "upper")
#######################################################################################
 
### ***** Analysis of TEs and the DMRs that occur in them

##################################################################################
##Are the TEs in centromeric regions different to TEs elsewhere in the genome???##
##################################################################################

#Is a TE centromeric
TEs$is.cent <- TEs$domain == "Centromeric"
TEs$is.cent <- as.numeric(TEs$is.cent)
TEs$len <- TEs$end - TEs$start
#Extracting TE, type and class
TE.attr <- strsplit((TEs$attributes), ";")
TE.type <- unname(sapply(TE.attr, '[[', 1))
TE.type <- strsplit(TE.type, "'")
TEs$TEtype <- unname(sapply(TE.type, '[[', 2))
TE.class <- unname(sapply(TE.attr, '[[', 2))
TE.class <- strsplit(TE.class, "'")
TEs$TEclass <- unname(sapply(TE.class, '[[', 2))

TE.summary <- summarise(group_by(TEs, is.cent, TEclass), len = sum(len), count = n())
print(TE.summary, n = 50)

#Calculating the expected number of TE classes
TE.counts.cent <- sum(filter(TE.summary, is.cent == 1)$count)
TE.counts.noncent <- sum(filter(TE.summary, is.cent == 0)$count)
TE.counts.total <- sum(TE.summary$count)
TEperc <- round(( TE.summary$count / c( rep(TE.counts.noncent, 22), rep(TE.counts.cent, 16) )) * 100,1)
TE.summary$TEperc <- TEperc

#TE.summary2 <- summarise(group_by(TEs, is.cent, TEtype), len = sum(len), count = n())
#print(TE.summary2, n = 50)

##Make a figure of the TEs in centromeric and non-centromeric regions

TE.summary.plot <- ggplot(filter(TE.summary, TEperc >= 0.1), aes(x = TEclass, y = TEperc, fill = factor(is.cent))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(data = filter(TE.summary, TEperc >= 0.1), aes(x = TEclass, y = TEperc, label = paste(TEperc,"%",sep = "")), position = position_dodge(width = 1), colour = ifelse(filter(TE.summary, TEperc >= 0.1)$is.cent, "#E41A1C", "#377EB8"), vjust = -0.5 ) +    
    scale_y_continuous(expand = c(0,0), limits = c(0,60)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = c("1" = "#E41A1C", "0" = "#377EB8"), labels = c("centromeric regions", "non-centromeric regions")) +    
    ylab("% of TEs") +
    xlab("") +
    theme(legend.position = "top", legend.justification = c(0.5, 0.5), legend.title = element_blank())

### To do ###!!
#Check are certain TE classed enriched for DMRs

#Take all the DMRs that occur in TEs
DMRs.TE <- filter(DMRs.all, anno.type == "TE")
DMRs.TE$centromeric <- DMRs.TE$domain == "Centromeric"

#Count DMRs for different domains and annotations
DMR.TE.counts <- summarise(group_by(DMRs.TE, TE.class, centromeric, .drop = FALSE), count = n())
DMR.TE.counts <- as.data.frame(DMR.TE.counts)

TE2 <- summarise(group_by(TEs, TEclass, is.cent), len = sum(len), count = n())
TE2 <- as.data.frame(TE2)

#Need to the amount of sequence and expected values

#Some TE classes (that are rare) did not have any DMRs in them (0 DMRs)
#For some reason I don't get the behaviour I expect, with .drop = FALSE, have to do it manually
TE2$DMRcount <- 0
for(i in 1:nrow(DMR.TE.counts)) {
    index <- which(DMR.TE.counts$TE.class[i] == TE2$TEclass & DMR.TE.counts$centromeric[i] == TE2$is.cent) #Check the index
    TE2$DMRcount[index] <- DMR.TE.counts$count[i]
}

#Calculating expected numbers
TE2$DMRexp <- sum(TE2$DMRcount) * ( TE2$len / sum(TE2$len))

#Fit a model of TE class effects across the whole genome
malliTE <- glm(DMRcount ~ -1 + offset(log(DMRexp)) +  TEclass, family = poisson, data = TE2)

#Storing the results
malliTE.results <- cbind(coef(malliTE), confint(malliTE))
malliTE.results <- exp(malliTE.results) #Change back to normal scale
colnames(malliTE.results)[1:3] <- c("estimate", "lower", "upper")
malliTE.results <- data.frame(malliTE.results)
malliTE.results$TEclass <- filter(TE2, is.cent == 0)$TEclass

#Fit the model using sequence lenghts as offset, for prediction
malliTE2 <- glm(DMRcount ~ -1 + offset(log(len)) +  TEclass, family = poisson, data = TE2)

#Plot the results
#Two very small classes with 0 DMRs need to be dropped
malliTE.results <- malliTE.results[-c(7,21),]
    
DMR.TE.rel.rates <- ggplot(malliTE.results, aes( x = TEclass, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange(pch = 1) +
    geom_hline(yintercept = 1, lty = "dashed") +
    xlab("") +
    ylab("DMRs relative to expected") +
    scale_y_log10() +
    scale_x_discrete(guide = guide_axis(angle = 45))   

## Use the inferred DMR rates in different TE's to check if they can predict the amount of TE's in centromeric and non-centromeric regions
pred.data <- data.frame(TEclass = TE2$TEclass, len = TE2$len, is.cent = TE2$is.cent)

DMRpred <- predict(malliTE2, newdata = pred.data, type = 'response')
TE2$DMRpred <- DMRpred

TEpred.plot <- summarise(group_by(TE2, is.cent), observed = sum(DMRcount), predicted = sum(DMRpred))
TEpred.plot <- pivot_longer(TEpred.plot, c("observed", "predicted"))
TEpred.plot$is.cent <- factor(TEpred.plot$is.cent, labels = c("non-centromeric", "centromeric"))
colnames(TEpred.plot) <- c("region", "type", "count")

TEpred.fig <- ggplot(TEpred.plot, aes(x = region, y = count, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_y_continuous(expand = c(0,0)) +
    ylab("DMR count in TEs") +
    xlab("") +
    theme(legend.position = "bottom", legend.justification = c(0.5, 0.5), legend.title = element_blank())

#First align plots
aligned <- align_plots(TE.summary.plot, DMR.TE.rel.rates, align = "v", axis = "l")
bottom.row <- plot_grid(aligned[[2]], TEpred.fig, align = "h", axis = "b", ncol = 2, labels = c("B", "C"), rel_widths = c(1, 0.3))
#Make the final plot
TE.final <- plot_grid(aligned[[1]], bottom.row, nrow = 2, labels = c("A", ""))

#Save the TE summary plot
save_plot("./epigenet_ms/fig/TE_classes.pdf", TE.summary.plot,  base_width = 3.71*4)

save_plot("./epigenet_ms/fig/TE_DMR.pdf", TE.final,  base_width = 3.71*4, base_height = 3.71*2.4)

### ***** DMR gains and losses

#Check if gains of DMRs are more likely to be associated with expansion of a DMR into neighbouring regions, rather than a completely new DMR ##TO DO!

table(DMRs.all$DMR.type) #There are more gain events than losses
#Note that there must be so many more potential sites for gain, so that is why gains happen at slower rates

DMRs.cent <- filter(DMRs.all, domain == "Centromeric")
DMRs.h3k9 <- filter(DMRs.all, domain == "H3K9")

DMRs.cent[3:4,] #Here is a case where there was an expansion of the DMR in the same line (L36)
#And contraction of the DMR in line 13

#How to look at this?
#Split the data into lines (all time points, and the ancestor)
#Then could maybe more easily record contractions and expansions
#Would need to do this for the Nanopore data to have better resolution

##TO DO!###

##Check differential gains and losses in different regions, annotations
DMR.types <- summarize(group_by(DMRs.all, domain, anno.type, DMR.type), count = n())

DMR.gains <- filter(DMR.types, DMR.type == "gain")
DMR.gains <- arrange(DMR.gains, anno.type)
DMR.gains$sequence <- DMR.counts$sequence
DMR.gains$expected <- sum(DMR.gains$count) * (DMR.gains$sequence / sum(DMR.gains$sequence) )

malli.gain <- glm(count ~ -1 + offset(log(expected)) + domain:anno.type, family = poisson, data = DMR.gains)

malli.gain.results <- cbind(coef(malli.gain), confint(malli.gain)) #Extract coefficients and conf int
malli.gain.results <- exp(malli.gain.results) #Change back to normal scale
malli.gain.results <- data.frame(malli.gain.results) #Change to data frame

malli.gain.results$domain <- factor(rep(c( "Centromeric", "Euchromatin", "H3K9", "H3K27"),4))
malli.gain.results$domain <- factor(malli.gain.results$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
malli.gain.results$annotation <- factor(c(rep("TE", 4), rep("gene", 4), rep("intergenic",4), rep("promoter", 4)))
colnames(malli.gain.results)[1:3] <- c("estimate", "lower", "upper")

DMR.losses <- filter(DMR.types, DMR.type == "loss")
DMR.losses <- arrange(DMR.losses, anno.type)
DMR.losses$sequence <- DMR.counts$sequence
DMR.losses$expected <- sum(DMR.losses$count) * (DMR.losses$sequence / sum(DMR.losses$sequence) )

malli.loss <- glm(count ~ -1 + offset(log(expected)) + domain:anno.type, family = poisson, data = DMR.losses)

malli.loss.results <- cbind(coef(malli.loss), confint(malli.loss)) #Extract coefficients and conf int
malli.loss.results <- exp(malli.loss.results) #Change back to normal scale
malli.loss.results <- data.frame(malli.loss.results) #Change to data frame

malli.loss.results$domain <- factor(rep(c( "Centromeric", "Euchromatin", "H3K9", "H3K27"),4))
malli.loss.results$domain <- factor(malli.loss.results$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
malli.loss.results$annotation <- factor(c(rep("TE", 4), rep("gene", 4), rep("intergenic",4), rep("promoter", 4)))
colnames(malli.loss.results)[1:3] <- c("estimate", "lower", "upper")

combined.results <- data.frame( rbind(malli.gain.results, malli.loss.results), DMRtype = c( rep("DMR gain", nrow(malli.gain.results)), rep("DMR loss", nrow(malli.loss.results))))

#Calculating the percentages
DMR.gains$percent <- rep(0,nrow(DMR.gains))
DMR.gains[c(1,5,9,13),]$percent <- round(DMR.gains[c(1,5,9,13),]$count / sum(DMR.gains[c(1,5,9,13),]$count) * 100,1) #Centromeric
DMR.gains[c(2,6,10,14),]$percent <- round(DMR.gains[c(2,6,10,14),]$count / sum(DMR.gains[c(2,6,10,14),]$count) * 100,1) #Euchromatin
DMR.gains[c(3,7,11,15),]$percent <- round(DMR.gains[c(3,7,11,15),]$count / sum(DMR.gains[c(3,7,11,15),]$count) * 100,1) #H3K27
DMR.gains[c(4,8,12,16),]$percent <- round(DMR.gains[c(4,8,12,16),]$count / sum(DMR.gains[c(4,8,12,16),]$count) * 100,1) #H3K9

DMR.losses$percent <- rep(0,nrow(DMR.losses))
DMR.losses[c(1,5,9,13),]$percent <- round(DMR.losses[c(1,5,9,13),]$count / sum(DMR.losses[c(1,5,9,13),]$count) * 100,1) #Centromeric
DMR.losses[c(2,6,10,14),]$percent <- round(DMR.losses[c(2,6,10,14),]$count / sum(DMR.losses[c(2,6,10,14),]$count) * 100,1) #Euchromatin
DMR.losses[c(3,7,11,15),]$percent <- round(DMR.losses[c(3,7,11,15),]$count / sum(DMR.losses[c(3,7,11,15),]$count) * 100,1) #H3K27
DMR.losses[c(4,8,12,16),]$percent <- round(DMR.losses[c(4,8,12,16),]$count / sum(DMR.losses[c(4,8,12,16),]$count) * 100,1) #H3K9

DMR.combined <- rbind(DMR.gains, DMR.losses)

#Making plots
dmr.comb.rates <- ggplot(combined.results, aes(x = domain, colour = annotation, y = estimate,  ymin = lower, ymax = upper)) +
    geom_pointrange(position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 1, lty = "dashed", col = "black") +
    xlab("") +
    ylab("Relative to expected") +
    scale_y_log10() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +    
    #scale_colour_brewer(type = "qual", palette = "Set1") +
    scale_colour_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +
    facet_wrap(~ DMRtype, nrow = 2) +
    theme(legend.position = "none")

dmr.gain.rates <- ggplot(malli.gain.results, aes(x = domain, colour = annotation, y = estimate,  ymin = lower, ymax = upper)) +
    geom_pointrange(position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 1, lty = "dashed", col = "black") +
    xlab("") +
    ylab("Relative to expected") +
    scale_y_log10() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +    
    #scale_colour_brewer(type = "qual", palette = "Set1") +
    scale_colour_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +
    theme(legend.position = "none")

dmr.loss.rates <- ggplot(malli.loss.results, aes(x = domain, colour = annotation, y = estimate,  ymin = lower, ymax = upper)) +
    geom_pointrange(position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 1, lty = "dashed", col = "black") +
    xlab("") +
    ylab("Relative to expected") +
    scale_y_log10() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +    
    #scale_colour_brewer(type = "qual", palette = "Set1") +
    scale_colour_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +
    theme(legend.position = "none")

DMR.gains.plot <-  ggplot(DMR.gains, aes(fill = anno.type, x = anno.type, y = count)) +
    geom_bar(position = position_dodge(), stat = "identity") +
    geom_text(data = DMR.gains, aes(x = anno.type, y = count, label = paste(percent,"%",sep = "")), vjust = -0.5 ) +
    scale_y_continuous(expand = expansion(mult = c(0,.1))) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +    
    #scale_fill_brewer(type = "qual", palette = "Set1") +
    ylab("DMR gains") +
    xlab("") +
    facet_wrap( ~ domain,, nrow = 1, scales = "free") +
    theme(legend.title = element_blank())

DMR.losses.plot <-  ggplot(DMR.losses, aes(fill = anno.type, x = anno.type, y = count)) +
    geom_bar(position = position_dodge(), stat = "identity") +
    geom_text(data = DMR.losses, aes(x = anno.type, y = count, label = paste(percent,"%",sep = "")), vjust = -0.5 ) +
    scale_y_continuous(expand = expansion(mult = c(0,.1))) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +    
    #scale_fill_brewer(type = "qual", palette = "Set1") +
    ylab("DMR losses") +
    xlab("") +
    facet_wrap( ~ domain,, nrow = 1, scales = "free") +
    theme(legend.title = element_blank()) 

dmr.gain.loss.plot <- plot_grid(DMR.gains.plot + theme(legend.position = "none"), dmr.gain.rates, DMR.losses.plot + theme(legend.position = "bottom", legend.justification = c(0.7, 0.5)), dmr.loss.rates, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"), rel_widths = c(4,1.5,4,1.5), align = "hv", axis = "lb")

save_plot("./epigenet_ms/fig/DMR_gain_loss.pdf", dmr.gain.loss.plot, base_height = 3.71*2, base_width = 3.71*3.5)

### ***** DMRs in euchromatin and proximity to TEs

#I want to check if DMRs that happen in euchromatin (genes and promoters) are closer to transposable elements that would be expected by a random gene or a promoter

DMRs.euchromatin <- filter(DMRs.all, domain == "Euchromatin")
DMRs.eu.genes <- filter(DMRs.euchromatin, anno.type == "gene" | anno.type == "promoter")
dist.TEs <- distance.TE(DMRs.eu.genes, TEs) #Check distances to the nearest TE
DMRs.eu.genes <- cbind(DMRs.eu.genes, dist.TEs)

#Need to establish a random baseline
#Sample random genes and promoters from euchromatin
eu.genes.counts <- summarize(group_by(DMRs.eu.genes, anno.type), count = n())
#

DMR.dist.means <- summarize(group_by(DMRs.eu.genes, closest.TE), mean.dist = mean(dist.TE), count = n()) #Distance from DMRs in euchromatic genes and promoters is a bit shorter on average to LTR elements than to unknown elements, or SINES
#Other classes have too few cases to analyse meaningfully

ggplot(DMRs.eu.genes, aes(x = dist.TE)) +
    geom_density() +
    scale_x_continuous(expand = c(0,0))


#    scale_x_log10()    

ggplot(DMRs.eu.genes, aes(x = closest.TE, y = dist.TE)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(angle = 45))    

#Need to establish a random baseline
#Sample random genes and promoters from euchromatin
genes.eu <- filter(genes, domain == "Euchromatin")
promoters.eu <- filter(promoters, domain == "Euchromatin")
foo <- list(0)
reslist <- rep(foo, 100)
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 100, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
#This takes a long time, use a progress bar to monitor
for(i in 1:100) { #Make 100 simulations, sampling the same number of genes and promoters from euchromatin 
    gene.index <- sample(1:nrow(genes.eu), size = as.numeric(eu.genes.counts[1,2]), replace = FALSE)
    sampled.genes <- genes.eu[gene.index, c(1,4,5)]
    prom.index <- sample(1:nrow(promoters.eu), size = as.numeric(eu.genes.counts[2,2]), replace = FALSE)
    sampled.promoters <- promoters.eu[prom.index, c(1,4,5)]
    sampled.features <- rbind(sampled.genes, sampled.promoters)
    #
    dist.sample <- distance.TE.feat(sampled.features, TEs)
    #
    sampled.features <- cbind(sampled.features,dist.sample)
    reslist[[i]] <- sampled.features
    #Update progress bar, this should be inside the loop
    setTxtProgressBar(pb, i)
}

#Save the simulations results, as this seems to take some time
save(DMRs.eu.genes, reslist, file = "./data/TEdistances.RData")

TEdens.plot <- ggplot() +
    geom_density(data = DMRs.eu.genes, aes(x = dist.TE), lwd = 1.5) + 
    lapply(reslist, function(dat) {
        geom_density(data = dat, aes(x = dist.TE), colour = alpha("blue", 0.1)) } ) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    xlab("Distance to nearest TE (bp)")    
        
save_plot("./epigenet_ms/fig/TEdens.pdf", TEdens.plot, base_height = 3.71)

load("./data/domains.RData")


### ***** Make a new figure of the DMR locations
load("ForChrPlot.RData")


Chr_DMR.plot <-  ggplot(chr.sizes, aes()) +
     geom_rect( aes(xmin = Start, xmax = End, ymin = 0, ymax = 5, fill = "Chromosome"), alpha = 0.5) +
     geom_rect(data = cent, aes(xmin = start, xmax = end, ymin = -1, ymax = -3, fill = "Centromeric")) +
     geom_rect(data = h3k27.plot, aes(xmin = Start, xmax = End, ymax = -4, ymin = -6, fill = "H3K27me3")) +
     geom_rect(data = h3k9.plot, aes(xmin = Start, xmax = End, ymax = -7, ymin = -9, fill = "H3K9me3")) +
     #geom_rect(data = h3k36.plot, aes(xmin = Start, xmax = End, ymax = -10, ymin = -13, fill = "H3K36me")) +
     geom_rect(data = DMRs.all, aes(xmin = start, xmax = end, ymax = 5, ymin = 0, fill = "DMR")) +
     #geom_rect(data = h3k9.ex.cent.plot, aes(xmin = Start, xmax = End, ymax = -10, ymin = -12, fill = "hotpink")) +
     #geom_tile(data = aineisto, aes(x = Position, y = 2.5, height = 5, width = 2), colour = "black", alpha = 0.5) +    
     xlab("Position (bp)") +
     ylab("") +
     scale_fill_manual(breaks = c("Chromosome", "Centromeric", "H3K27me3", "H3K9me3", "DMR"), values = c(Chromosome = "deepskyblue", Centromeric = "grey", H3K27me3 = "blue", H3K9me3 = "red", DMR = "black")) +
     scale_x_continuous(expand = expansion(mult = c(0, 0.05)), breaks = c(0,2.5e6,5e6,7.5e6,1e7), labels = c(TeX("$0$"), TeX("$2.5 \\times 10^6$"), TeX("$5 \\times 10^6$"), TeX("$7.5 \\times 10^6$"), TeX("$1 \\times 10^7$")), limits = c(0, 1e7))  +
     guides(fill = guide_legend(nrow = 1)) +
     facet_grid(Chromosome ~ ., labeller = labeller(Chromosome = labels)) +
     theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "top", legend.justification = c(0.5, 0.5), legend.title = element_blank())

save_plot("./epigenet_ms/fig/chromosomes_DMRs.pdf", Chr_DMR.plot, base_height = 6, base_width = 9*1.618)

#Plot annotation types in the different domains
dmr.counts1 <- ggplot(DMR.counts, aes(fill = anno.type, x = domain, y = count)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    #scale_fill_brewer(type = "qual", palette = "Set1") +
    scale_fill_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +    
    ylab("DMR count") +
    xlab("") +    
    theme(legend.position = "none")

dmr.counts1.2 <- ggplot(DMR.domains, aes(x = domain, y = count)) +
    geom_bar(stat = "identity", colour = "black", fill = "lightgrey") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    #scale_fill_brewer(type = "qual", palette = "Set1") +
    ylab("DMR count") +
    xlab("") +    
    theme(legend.position = "none")

dmr.counts1 <- ggplot(DMR.counts, aes(fill = anno.type, x = domain, y = count)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_y_continuous(expand = c(0,0)) +
    #scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +
    #scale_fill_brewer(type = "qual", palette = "Set1") +
    ylab("DMR count") +
    xlab("") +
    facet_wrap( ~ domain, ncol = 1, scales = "free") +    
    coord_flip() +
    theme(legend.position = "none")


dmr.counts2.2 <- ggplot(DMR.counts, aes(fill = anno.type, x = anno.type, y = count)) +
    geom_bar(position = position_dodge(), stat = "identity") +
    geom_text(data = DMR.counts, aes(x = anno.type, y = count, label = paste(percent,"%",sep = "")), vjust = -0.5 ) +
    scale_y_continuous(expand = expansion(mult = c(0,.1))) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +    
    #scale_fill_brewer(type = "qual", palette = "Set1") +
    ylab("DMR count") +
    xlab("") +
    facet_wrap( ~ domain, nrow = 1, scales = "free") +
    theme(legend.title = element_blank()) 

#Same plot but with stacked proportions
dmr.counts2 <- ggplot(DMR.counts, aes(fill = anno.type, x = domain, y = count)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +        
    #scale_fill_brewer(type = "qual", palette = "Set1") +
    geom_text(aes(label = count), position = position_fill(vjust = 0.5)) +  # Show counts   
    ylab("Proportion") +
    xlab("") +
    theme(legend.title = element_blank())

legend <- get_legend(dmr.counts2)

legend_b <- get_legend(dmr.counts2 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.justification = c(0.5, 0.5))       
)

#Plotting
dmr.rel.rates <- ggplot(malli1.results, aes( x = domain, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange(pch = 1) +
    geom_hline(yintercept = 1, lty = "dashed") +
    xlab("") +
    ylab("Relative to expected") +
    scale_y_log10() +
    scale_x_discrete(guide = guide_axis(angle = 45))   

#brewer.pal(4, "Set1")  "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3"

dmr.rates.anno <- ggplot(malli3.results, aes( x = domain, colour = annotation, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange(position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 1, lty = "dashed", col = "black") +
    xlab("") +
    ylab("Relative to expected") +
    scale_y_log10() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +    
    #scale_colour_brewer(type = "qual", palette = "Set1") +
    scale_colour_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) + 
    theme(legend.position = "none")   

bottom.row <- plot_grid(dmr.counts1, dmr.counts2 + theme(legend.position = "none") , dmr.rel.rates, dmr.rates.anno, legend, nrow = 1, rel_widths = c(0.7,0.8,0.5, 0.7, 0.3), align = "v", axis = 'l', labels = c("B", "C", "D", "E"))

bottom.row2 <- plot_grid(dmr.counts1.2, dmr.rel.rates, dmr.counts2.2 + theme(legend.position = "none") ,  dmr.rates.anno, nrow = 1, rel_widths = c(0.25, 0.25, 1, 0.3), align = "hv", axis = 'bt', labels = c("B", "C", "D", "E"))

final.DMR.plot <- plot_grid(Chr_DMR.plot, bottom.row, nrow = 2, ncol = 1, labels = "A", rel_heights = c(1, 0.7))
save_plot("./epigenet_ms/fig/chromosomes_DMRs_final.pdf", final.DMR.plot, base_height = 9, base_width = 9*1.618)

final.DMR.plot2 <- plot_grid(Chr_DMR.plot, bottom.row2, legend_b, nrow = 3, ncol = 1, labels = "A", rel_heights = c(1, 0.6, 0.04))
save_plot("./epigenet_ms/fig/chromosomes_DMRs_final3.pdf", final.DMR.plot2, base_height = 10, base_width = 9.2*1.618)

save_plot("./epigenet_ms/fig/DMR_summary.pdf", plot_grid(bottom.row2, legend_b, nrow = 2, rel_heights = c(1,0.05)), base_height = 4, base_width = 9*1.618)


#For presentations etc.
save_plot("./epigenet_ms/fig/DMR_composition.pdf", bottom.row, base_height = 4, base_width = 9*1.618)


### ** DMRs that occurred in the Nanopore sample

### *** Load the data, combine different contexts, and annotate

DMRs.all.CG <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/nanopore_all/CG_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CG$context <- "CG"

DMRs.all.CHG <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/nanopore_all/CHG_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHG$context <- "CHG"

DMRs.all.CHH <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/nanopore_all/CHH_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHH$context <- "CHH"

DMRs.all <- rbind(DMRs.all.CG, DMRs.all.CHG, DMRs.all.CHH)
DMRs.all <- arrange(DMRs.all, seqnames, start)

#There are some bins that have been called in multiple contexts
bincheck <- paste0(DMRs.all[,1],DMRs.all[,2])
duplicatebins <- duplicated(bincheck)
#DMRs.all[duplicatebins,c(1:3,74)]

#Filter bins that are in different contexts (that is get only unique locations)
DMRs.all.nanopore <- DMRs.all[!duplicatebins,]

#Filter everything that is not in the seven chromosomes
DMRs.all.nanopore <- filter(DMRs.all.nanopore, seqnames %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Load annotations
genes <- read.gff("~/Genomics/Neurospora/reference/Neurospora_crassa_NC12_genes.gff3")
genes <- filter(genes, type == "gene")
#We are only considering the seven chromosomes
genes <- filter(genes, seqid %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7") == T)
genes$seqid <- factor(genes$seqid) #Drop levels that are ot used

#Load TE annotations
TEs <- read.gff("~/Genomics/Neurospora/reference/Neurospora_crassa_NC12_TE.gff3")

#Promoter annotations, promdist the distance from TSS that still counts as promoter region
promoters <- annotate.promoters(genes, promdist = 500)

#This annotates all DMRs
DMRs.all.nanopore <- annotate.DMRs(DMRs.all.nanopore, genes, TEs, promoters) #Some nearby DMRs can occur in the same gene

## There are some warnings that occur
## Use options(warn=2), to convert warnings into errors and check what is the problem
## options(error = recover)
##

#Then we need information about the chromatin domains
#Load the data
#Euchromatin
euchr <- read.table("./annotation/2489.euchromatin.bed", header = F, sep = "\t")
colnames(euchr) <- c("Chromosome", "start", "end")
euchr <- filter(euchr, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Centromeric H3K9
centromeric.h3k9 <- read.table("./annotation/centromericH3K9.bed", header = F, sep = "\t")
colnames(centromeric.h3k9)[1:3] <- c("Chromosome", "start", "end")

#H3K9 not in centromeric regions
excent2.h3k9 <- read.table("./annotation/H3K9_ex_cent2.bed", header = F, sep = "\t")
colnames(excent2.h3k9)[1:3] <- c("Chromosome", "start", "end")
excent2.h3k9 <- filter(excent2.h3k9, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#H3K27 not in H3K9 regions
h3k27.exh3k9 <- read.table("./annotation/H3K27_exH3K9.bed", header = F, sep = "\t")
colnames(h3k27.exh3k9)[1:3] <- c("Chromosome", "start", "end")
h3k27.exh3k9 <- filter(h3k27.exh3k9, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Check DMR overlap for chromatin domains
DMRs.all.nanopore <- DMR.by.domains(DMRs.all.nanopore, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)

##Then need to check overlap for each gene and promoter in the different domains
genes <- annotation.by.domains(genes, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)
TEs <- annotation.by.domains(TEs, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)
promoters <- annotation.by.domains(promoters, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9, prom = TRUE)

#Then check whether methylation was gained or lost
#Load the methylation level files
DMRs.all.CG.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/nanopore_all/CG_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CG.meth$context <- "CG"

DMRs.all.CHG.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/nanopore_all/CHG_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHG.meth$context <- "CHG"

DMRs.all.CHH.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/nanopore_all/CHH_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHH.meth$context <- "CHH"

colnames(DMRs.all.nanopore)[1] <- "Chromosome"
#For each DMR check whether methylation was gained or lost
meth <- classify.DMRs(DMRs.all.nanopore, ind.cols = 5:42, DMRs.all.CG.meth, DMRs.all.CHG.meth, DMRs.all.CHH.meth)
DMRs.all.nanopore <- cbind(DMRs.all.nanopore, meth) #Store methylation levels, and type

#For some reason 

#Saving the data for easier loading the next time
save(DMRs.all.nanopore, euchr, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9, genes, TEs, promoters, file = "./data/DMRs_all_nanopore.RData")


### * H3K9 chip-seq data

### ** Peaks called by Mariana

### Loading some of the peaks called by Mariana

#Set up folder, samples, and file extension names
peakdatafolder <- "~/Genomics/Neurospora/chipseq/peaks/"
peaksamples.matA <- c("G5_L5", "G8_L5", "G10_L5", "G20_L5", "G40_L5", "G5_L11", "G8_L11", "G10_L11", "G20_L11", "G40_L11", "anc1_matA_rep1", "anc1_matA_rep2", "anc2_matA_rep1", "anc2_matA_rep2")
fileext <- "_2broad.bdg"
#length(peaksamples.matA) #There are 14 samples

matA.peaks <- load.chipseq.peaks(peakdatafolder, peaksamples.matA, fileext)

#Load the mata peaks
peaksamples.mata <- c("G5_L25", "G8_L25", "G10_L25", "G20_L25", "G40_L25", "G5_L31", "G8_L31", "G10_L31", "G20_L25", "G40_L31", "anc1_mata_rep1", "anc2_mata_rep1")

mata.peaks <- load.chipseq.peaks(peakdatafolder, peaksamples.mata, fileext)


##Check the length of the called regions
peak.len <- matA.peaks$end - matA.peaks$start
hist(log10(peak.len)) #Histogram of peak lengths
quantile(log10(peak.len), probs = c(0.025, 0.5, 0.975))
#Smallest peak length is 500 bp, and maximum is 353 kb
#We need to set up some overlap threshold
#Since some peaks are very long, this threshold probably needs to depend on the peak lenght itself
#i.e. some percentage overlap

#This function checks for peak overlap among all peaks called from the samples
matA.peaks.merged <- merge.peaks(matA.peaks, peaksamples.matA)
#length(unique(koe$peak.group)) #392 unique peaks

mata.peaks.merged <- merge.peaks(mata.peaks, peaksamples.mata)
#length(unique(mata.peaks.merged$peak.group)) #401 unique peaks

#Next for each unique peak, get chr, median start, median end coordinates, and make a matrix of presence / absence for each sample.
matA.peaks.annotated <- annotate.peaks(matA.peaks.merged, peaksamples.matA)

mata.peaks.annotated <- annotate.peaks(mata.peaks.merged, peaksamples.mata)

#Then filter for variable peaks and and mark those that need to be checked in IGV
matA.variable <- filter(matA.peaks.annotated, is.polymorphic == T)
nrow(matA.variable) #There are 163 variable peaks, most should be checked in IGV

#Filter for mata peaks
mata.variable <- filter(mata.peaks.annotated, is.polymorphic == T)
nrow(mata.variable)

#Total number of variable peaks
nrow(matA.variable) + nrow(mata.variable)

##Save the annotated and variable peaks for manual examination in IGV
savefolder <- "~/Documents/tutkijatohtori/epimutation/MA_meth/data/"
write.table(matA.peaks.annotated, file = paste0(savefolder, "H3K9_matA_peaks_all.csv"), sep = ",", row.names = F, quote = F)
write.table(matA.variable, file = paste0(savefolder, "H3K9_matA_peaks_variable.csv"), sep = ",", row.names = F, quote = F)

write.table(mata.peaks.annotated, file = paste0(savefolder, "H3K9_mata_peaks_all.csv"), sep = ",", row.names = F, quote = F)
write.table(mata.variable, file = paste0(savefolder, "H3K9_mata_peaks_variable.csv"), sep = ",", row.names = F, quote = F)


### *** Analysis of peaks

## Load the curated peak data files

#Alternatetively do the analysis with all variable peaks

matA.peak.matrix <- peak.matrix(variable = matA.variable, merged = matA.peaks.merged, samples =  peaksamples.matA)

mata.peak.matrix <- peak.matrix(variable = mata.variable, merged = mata.peaks.merged, samples =  peaksamples.mata)


#Then need to calculate divergence
#Need information about number of mitoses separating different samples

#All pairwise sample combinations
comb.matA <- t(combn(peaksamples.matA,2))
#calc.dt(koe[1,][1], koe[1,][2]) #This seems to work
#Use the calc.dt function to calculate divergence times between pairs of samples
dt.matA <- mapply(calc.dt, comb.matA[,1], comb.matA[,2]) #For applying the calc.dt function to each row
#In order to convert transfers into mitoses, (25 mitoses per transfer on average)
#Multiple by 25
div.matA <- calc.h3k9.div(matA.peak.matrix, comb.matA)

comb.mata <- t(combn(peaksamples.mata,2))
dt.mata <- mapply(calc.dt, comb.mata[,1], comb.mata[,2])
div.mata <- calc.h3k9.div(mata.peak.matrix, comb.mata)

divergence.matA <- data.frame(sample1 = comb.matA[,1], sample2 = comb.matA[,2], dt = dt.matA, dm = dt.matA*25, div.h3k9 = div.matA)

divergence.mata <- data.frame(sample1 = comb.mata[,1], sample2 = comb.mata[,2], dt = dt.mata, dm = dt.mata*25, div.h3k9 = div.mata)

divergence.final <- rbind(divergence.matA, divergence.mata)

myxlab <- TeX("$\\Delta$t (mitoses)")
myanno1 <- TeX("$\\beta = -2.66 \\times 10^{-5}$, $p = 0.386$")
h3k9div.plot <- ggplot(divergence.final, aes(y = div.h3k9, x = dm)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_y_continuous(limits = c(0,1)) +
    ylab("H3K9me3 divergence") +
    xlab(myxlab) +
    annotate("text", x = 500, y = 0.75, label = myanno1)

malli.h3k9 <- lm(div.h3k9 ~ dm, data = divergence.final)

save_plot("./epigenet_ms/fig/h3k9div.pdf", h3k9div.plot)
### *** H3K9me3 peaks and DMRs

#Load the H3K9me3 peak data

datafolder <- "~/Documents/tutkijatohtori/epimutation/MA_meth/data/"

#Load the peaks from above
#matA.peaks
#matA.peaks.merged

#Load DMRs from the nanopore dataset
#Load DMRs from the nanopore dataset
load(file = "./data/DMRs_all_nanopore.RData") #In DMRs.all.nanopore

samples.matA <- c("MatA", "G5_L5", "G8_L5", "G10_L5", "G5_L11", "G8_L11", "G10_L11")
samples.mata <- c("Mata", "G5_L25", "G8_L25", "G10_L25", "G5_L31", "G8_L31", "G10_L31")

matA.by.sample <- DMR.H3K9.overlap.by.sample(DMRs.all.nanopore, samples.matA, matA.peaks.merged)
mata.by.sample <- DMR.H3K9.overlap.by.sample(DMRs.all.nanopore, samples.mata, mata.peaks.merged)

#Saving some intermediate results
#save(matA.by.sample, mata.by.sample, file = "./data/DMRs_H3K9_bysample.RData")
#load("./data/DMRs_H3K9_bysample.RData")

#Putting by sample datasets together
DMRs.matA <- rbind(matA.by.sample[[1]][,-5], matA.by.sample[[2]][,-5], matA.by.sample[[3]][,-5], matA.by.sample[[4]][,-5], matA.by.sample[[5]][,-5], matA.by.sample[[6]][,-5], matA.by.sample[[7]][,-5])
DMRs.mata <- rbind(mata.by.sample[[1]][,-5], mata.by.sample[[2]][,-5], mata.by.sample[[3]][,-5], mata.by.sample[[4]][,-5], mata.by.sample[[5]][,-5], mata.by.sample[[6]][,-5], mata.by.sample[[7]][,-5])

#Combining and removing duplicate positions
DMRs.unique <- rbind(DMRs.matA, DMRs.mata)
pos <- paste(DMRs.unique[,1], DMRs.unique[,2], sep = " ") #Combine chr and position
uniq.index <- !duplicated(pos) #Get only unique DMRs, duplicate rows are removed
DMRs.unique <- DMRs.unique[uniq.index,]

#Analysis of DMR overlap
#Making some tables of DMR and H3K9 overlap
nrow(DMRs.unique) # 14 503 DMRs that happened in these samples (duplicates are removed)

round((table(DMRs.unique$H3K9.overlap)[2] / nrow(DMRs.unique)) * 100,2) #90.64% of DMRs overlap with H3K9

res.table <- summarize(group_by(DMRs.unique, H3K9.overlap, DMR.type), total = n())
res.table$H3K9.overlap <- c("no H3K9me3\n overlap", "no H3K9me3\n overlap", "H3K9me3\n overlap", "H3K9me3\n overlap")
res.table$percent <- round(c(1290/(1290+67), 67/(1290+67), 8441/(8441+4705), 4705/(8441+4705))*100,1)

DMR.H3K9.p <- ggplot(res.table, aes(x = H3K9.overlap, fill = DMR.type, y = total)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(data = res.table, aes(x = H3K9.overlap, y = total, label = paste(percent,"%",sep = "")), vjust = -0.5, position = position_dodge(width = 1)) +
    ylab("DMR count") +
    xlab("") +
    scale_y_continuous(expand = expansion(mult = c(0,.10))) +
    scale_fill_manual(values = c(gain = "#E41A1C", loss = "#377EB8")) +
    labs(fill = "Methylation") +
    theme(legend.position = "top", legend.justification = "center")    

res.table2 <- summarize(group_by(DMRs.unique, H3K9.overlap, anno.type), total = n())
res.table2$percent <- round(ifelse(res.table2$H3K9.overlap == F, res.table2$total/sum(res.table2$total[1:4]), res.table2$total/sum(res.table2$total[5:8]))*100 ,1) #Calculate percentages within overlap class
res.table2$H3K9.overlap <- c(rep("no H3K9me3\n overlap", 4), rep("H3K9me3\n overlap", 4))

DMR.H3K9.anno.p <- ggplot(res.table2, aes(x = anno.type, fill = anno.type, y = total)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(data = res.table2, aes(x = anno.type, y = total, label = paste(percent,"%",sep = "")), vjust = -0.5) +
    ylab("DMR count") +
    xlab("") +
    scale_y_continuous(expand = expansion(mult = c(0,.1))) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +    
    scale_fill_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +
    facet_wrap( ~ H3K9.overlap, scales = "free") +
    theme(legend.position = "none", legend.justification = "center")  

DMR.H3K9.final <- plot_grid(DMR.H3K9.p, DMR.H3K9.anno.p, nrow = 1, labels = c("A", "B"), align = "h", axis = "b", rel_widths = c(0.65,1))
save_plot(filename = "./epigenet_ms/fig/H3K9_DMR_overlap.pdf", DMR.H3K9.final, base_height = 4, base_width = 8)

DMR.gains.plot <-  ggplot(DMR.gains, aes(fill = anno.type, x = anno.type, y = count)) +
    geom_bar(position = position_dodge(), stat = "identity") +
    geom_text(data = DMR.gains, aes(x = anno.type, y = count, label = paste(percent,"%",sep = "")), vjust = -0.5 ) +
    scale_y_continuous(expand = expansion(mult = c(0,.1))) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = c(TE = "#E41A1C", intergenic = "#4DAF4A", gene = "#377EB8", promoter = "#984EA3")) +    
    #scale_fill_brewer(type = "qual", palette = "Set1") +
    ylab("DMR gains") +
    xlab("") +
    facet_wrap( ~ domain,, nrow = 1, scales = "free") +
    theme(legend.title = element_blank())


#testing
koe1 <- matA.by.sample[[2]]
koe2 <- matA.by.sample[[3]]
koe3 <- rbind(koe1[,-5],koe2[,-5]) #BOOKMARK!, now filter for unique DMRs
pos <- paste(koe3[,1], koe3[,2], sep = " ")
uniq.index <- !duplicated(pos)
koe3 <- koe3[uniq.index,]

pos[!duplicated(pos)]

uniq <- duplicated(pos) | duplicated(pos, fromLast = TRUE) #

a <- c(rep("A", 3), rep("B", 3), rep("C",2))
b <- c(1,1,2,4,1,1,2,2)
df <-data.frame(a,b)

duplicated(df)
[1] FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE  TRUE

> df[duplicated(df), ]
  a b
2 A 1
6 B 1
8 C 2

> df[!duplicated(df), ]
  a b
1 A 1
3 A 2
4 B 4
5 B 1
7 C 2

#Need to make a function that checks whether each DMR overlaps with H3K9 by sample
DMR.H3K9.overlap.by.sample <- function(DMRs, samples, H3K9) {

    res.list <- rep(list(1), length(samples))
    #First loop over samples that have both DMR and H3K9 data
    for(i in 1:length(samples)) { #Loop over all samples
        sample.index <- which(grepl(samples[i], colnames(DMRs))) #sample index
        sample.DMR <- DMRs[,c(1:4,sample.index, 43:59)] #Take DMRs of current sample
        index <- sample.DMR[,5] != sample.DMR$DMR.common #Filter for variable DMRs
        sample.DMR <- sample.DMR[index,]
        sample.DMR$H3K9.overlap <- NA #Initialize overlap vector

        #Then take the H3K9 peaks for the current sample
        #Need to adjust sample names because they are different in H3K9 data
        if(samples[i] == "MatA") { cur.sample <- "matA" }
        if(samples[i] == "Mata") { cur.sample <- "mata" }
        if(samples[i] != "Mata" & samples[i] != "MatA") { cur.sample <- samples[i] }
        index2 <- grepl(cur.sample, H3K9$sampleID)
        sample.H3K9 <- H3K9[index2,]

        #Then check overlap for each DMR
        for(j in 1:nrow(sample.DMR)) { #Loop over each DMR
            current.DMR <- sample.DMR[j,] #Take the current DMR
            check <- apply(sample.H3K9[,2:3], MARGIN = 1, check.overlap, start.dmr = current.DMR$start, end.dmr = current.DMR$end)
            sample.DMR$H3K9.overlap[j] <- any(check)
        } #Done looping over all DMRs of a sample

        #Store the results
        res.list[[i]] <- sample.DMR
    } #Done looping over all samples

    #return results
    return(res.list)
    
    #Save results, sample by sample
}


#From the sample by sample results, need to filter unique DMRs, so that there is no double counting

#Mat a ancestor
mata.anc.DMR <- DMRs.all.nanopore[,c(1:4,5,43:59)]
#Filter for variable DMRs
index <- mata.anc.DMR$Mata != mata.anc.DMR$DMR.common
mata.anc.DMR <- mata.anc.DMR[index,]

index2 <- grepl("mata", mata.peaks.merged$sampleID)
mata.anc.H3K9 <- mata.peaks.merged[index2,]

#Check overlap for each DMR
#For each DMR
current.DMR <- mata.anc.DMR[1,]

check <- apply(mata.anc.H3K9[,2:3], MARGIN = 1, check.overlap, start.dmr = current.DMR$start, end.dmr = current.DMR$end)
any(check)
#Then store information about whether a DMR overlaps H3K9

#matA.peaks <- read.csv(file = paste(datafolder, "H3K9_matA_peaks_all.csv", sep = ""), sep = ",", header = T)

### * Genetic mutations and methylation changes

#Number of cytosine mutations in the MA-lines compared to number of single cytosine methylation changes

### ** Number of genetic mutations

genmutations <- read.table("~/Documents/tutkijatohtori/epimutation/MA_WGS/curated_mutations_final.csv", header = T, sep = ",")

mutations.lines <- group_by(genmutations, Line)
mutperline <- summarise(mutations.lines, nmut = n()) #Mutations per line

#Filtering for only lines that were bisulfite sequenced

mutperline <- mutperline[c(11,36,37,1,4,5,6,8,10,12,13,20,21,23,25,26,27,28,31,33),]
sum(mutperline[,2]) #725 genetic mutations in total


#Only point mutations that involved a cytosine as the ancestral of the derived base
Cmutations <- filter(genmutations, type == "point" & (anc.base == "C" | sample.base == "C"))

mutations.C.lines <- group_by(Cmutations, Line)
C.mutperline <- summarise(mutations.C.lines, nmut = n())

C.mutperline <- C.mutperline[c(11,36,37,1,4,5,6,8,10,12,13,20,21,23,25,26,27,28,31,33),]
sum(C.mutperline[,2]) #301 point mutations involving a cytosine

### ** Number of single cytosince changes in the lines

#Need to make a dataset that has all variable cytosines

### *** Testing some scripts

datafolder <- "~/Genomics/Neurospora/methylation/methimpute/all/"

#koe1 <- fread(paste0(datafolder,"methylome_MANC1_bismark_pe_trim_All.txt"), select = c("seqnames", "start", "status"))

#koe2 <- fread(paste0(datafolder,"methylome_ML1G40_bismark_pe_trim_All.txt"), select = c("seqnames", "start", "status"))

#any(koe1$start != koe2$start) #All positions are the same

#sum(koe$status != koe2$status)

#File list for testing my function

mylist <- c("methylome_MANC1_bismark_pe_trim_All.txt", "methylome_ML1G40_bismark_pe_trim_All.txt", "methylome_MG5L36_bismark_pe_trim_All.txt")

#Making a function that counts the number of variable cytosines

methylome.cytosines <- function(filelist, datafolder) {

    nfiles <- length(filelist) #Number of files to process

    #Set up the dataframe by loading the first sample
    if(grepl("R", filelist[1]) == T) {
        pieces <- unlist(strsplit(filelist[1], split = "_"))
        sample.name <- paste0(pieces[2], "_", pieces[3]) } else {
            sample.name <- unlist(strsplit(filelist[1], split = "_"))[2]
    } #Get the sample name from file name
    
    aineisto <- fread(paste0(datafolder,filelist[1]), select = c("seqnames", "start", "status"))
    setnames(aineisto, "status", sample.name) #Change the name of the variable   
          
    #Loop over all files
    for(i in 2:nfiles) {
            if(grepl("R", filelist[1]) == T) {
                pieces <- unlist(strsplit(filelist[1], split = "_"))
                sample.name <- paste0(pieces[2], "_", pieces[3]) } else {
                    sample.name <- unlist(strsplit(filelist[1], split = "_"))[2]
                } #Get the sample name from file name #Name of the current sample
        current <- fread(paste0(datafolder,filelist[i]), select = "status") #methylation status for sample i
        setnames(current, "status", sample.name) #Change the name to current sample
        #Combine to the original data table
        aineisto <- cbind(aineisto, current)
    } #Done looping over all files
   
    return(aineisto)
}

#This function check if a cytosine is variable
variable.cytosines <- function(data) {
    !(all(data == "U") | all(data == "M") | all(data == "I")) #TRUE site is variable, FALSE site is not variable
}

#Assemble the methylation status of each cytosine for all samples
methylomes <- methylome.cytosines(filelist = mylist, datafolder = datafolder)

#Then check whether each position is variable (first two columns contain chr and position: dropped)
variable <- apply(methylomes[,-c(1,2)], 1, variable.cytosines)
#Combine
methylomes <- cbind(methylomes, variable)

#Writing the data table to a file
fwrite(methylomes, file = paste0(datafolder, "methylome_cytosines_all_samples.txt"), sep = "\t")

#Putting some results to a table
res.mat.ss <- data.frame(total.C = nrow(methylomes), variable.C = sum(methylomes$variable))
res.mat.ss$varperc <- round((res.mat.ss[1,2] / res.mat.ss[1,1]) * 100 ,2) #Percentage of variable sites

#Save the data
save(res.mat.ss, file = "single.cytosine.variable.results.RData")

variable <- rep(FALSE, nrow(koe))
for(i in 1:nrow(koe)) { variable[i] <- variable.cytosines(koe[i,]) }

variable.cytosines <- function(data) {
    check <- !is.na(data) #Nedd to check for missing data
    if(any(check) == T) {
    data <- data[, ..check]
    return(!(all(data == "Unmethylated") | all(data == "Methylated"))) #TRUE site is variable, FALSE site is not variable
} else { return(NA) } #Return NA if all data is missing
}


### Complete R script was run on the CSC cluster
a <- "methylome_MANC1_R2_bismark_pe_trim_All.txt"

strsplit(a, "\\(R*\\)(*SKIP)(*F)|\\h*_\\h*", perl=T)

x <- "This is it, isn't it (well, yes), and (well, this, that, and this, too)"
strsplit(x, "\\([^()]*\\)(*SKIP)(*F)|\\h*,\\h*", perl=T)
### *** Analysing the full single cytosine methylomes

#Loading the data for the bisulfite sequencing set, all samples
load("~/Genomics/Neurospora/methylation/methimpute/all/single.cytosine.variable.results.RData")

methylomes <- fread("~/Genomics/Neurospora/methylation/methimpute/all/methylome_cytosines_all_samples.txt")

#Fixing the amount of variable sites, since need to use na.rm = T
res.mat.ss[2] <- sum(methylomes$variable, na.rm = T)
res.mat.ss$varperc <- round((res.mat.ss[1,2] / res.mat.ss[1,1]) * 100 ,2)

save(res.mat.ss, file = "~/Documents/tutkijatohtori/epimutation/MA_meth/data/single.cytosine.variable.results.RData")

load("~/Documents/tutkijatohtori/epimutation/MA_meth/data/single.cytosine.variable.results.RData")

### *** Single cytosine methylation of the nanopore dataset

datafolder <- "~/Genomics/Neurospora/methylation/methimpute/nanopore_all/"

#koe1 <- fread(paste0(datafolder,"2489_MatA.txt"), select = c("seqnames", "start", "status"))

mylist <- c("2489_MatA.txt")



### ** Number of genetic mutations in methylation changes in lines with nanopore data

### *** Single cytosine changes and genetic mutations

### **** Script for determining the number of cytosine changes in each transfer
##For a script to be run on the cluster
library(tidyverse)
library(data.table)

##Need to make a function that checks cytosine methylation changes in all intervals
intervals.cytosine.check <- function(aineisto, cytosines.all.nanopore) {
    ind.cols <- 3:40 #This is currently hard coded, should change this
    inds <- cytosines.all.nanopore[,..ind.cols]
    ind.names <- colnames(inds)
    mylist <- strsplit(ind.names, "_")
    gens <- unname(sapply(mylist, '[[', 1))
    lines <- unname(sapply(mylist, '[[', 2))

    #Loop over all lines and intervals
    intnum <- nrow(aineisto) #Number of intervals to check
    #Set up the results matrix
    res.mat <- matrix(rep(0, intnum*3), ncol = 3)
    colnames(res.mat) <- c("Total", "Gains", "Losses")
    
    for(i in 1:intnum) { #loop over all intervals
        current.line <- aineisto$Line[i]
        current.gen <- aineisto$Generation[i]

        #If current generation is G1 then we are comparing to the ancestor and methylation
        #comparison is relative to ancestor
        if(current.gen == "G1") {
            #Need to get the correct column index for inds
            cur.ind <- gens == current.gen & lines == current.line
            cur.ind <- which(cur.ind)
            current <- inds[,..cur.ind] #Current focal sample

            if(current.line %in% c("L2", "L5", "L11")) {
               changes <- current != inds[,2] #Comparing to mat A 
               types <- ifelse(inds[,2] == "Unmethylated" & current == "Methylated", "gain", "loss")  }
            
            if(current.line %in% c("L23", "L25", "L31")) {         
               changes <- current != inds[,1] #Comparing to mat a
               types <- ifelse(inds[,1] == "Unmethylated" & current == "Methylated", "gain", "loss")  }
                
            res.mat[i,1] <- sum(changes, na.rm = T) #Total number of changes in interval
            res.mat[i,2] <- sum(changes == TRUE & types == "gain", na.rm = T) #Number of DMR gains
            res.mat[i,3] <- sum(changes == TRUE & types == "loss", na.rm = T) #Number of DMR losses
        } #Done counting DMR changes

        #If not in first interval, need to compare with previous
        if(current.gen != "G1") {
            #Get the current ind
            cur.ind <- gens == current.gen & lines == current.line
            cur.ind <- which(cur.ind)
            current <- inds[,..cur.ind] #Current focal sample
            #Check which was the previous transfer
            if(current.gen == "G5") { previous.gen <- "G1" }
            if(current.gen == "G7") { previous.gen <- "G5" }
            if(current.gen == "G8") { previous.gen <- "G7" }
            if(current.gen == "G10") { previous.gen <- "G8" }
            if(current.gen == "G15") { previous.gen <- "G10" }
            gens == current.gen & lines == current.line
            prev.ind <- gens == previous.gen & lines == current.line
            prev.ind <- which(prev.ind)
            previous <- inds[,..prev.ind] #Sample from previous transfer

           
                changes <- current != previous
                types <- ifelse(previous == "Unmethylated" & current == "Methylated", "gain", "loss")
            res.mat[i,1] <- sum(changes, na.rm = T) #Total number of DMR changes
            res.mat[i,2] <- sum(changes == TRUE & types == "gain", na.rm = T) #Number of methylation gains
            res.mat[i,3] <- sum(changes == TRUE & types == "loss", na.rm = T) #Number of methylation losses
        } #Done counting DMR changes
    } #Done looping over all intervals

    #Return results
    return(cbind(aineisto, res.mat))
} #Done


#First load the data
datafolder <- "/scratch/project_2000350/genomics/methylation/methimpute/ss_methylome/"

#datafolder <- "~/Documents/tutkijatohtori/epimutation/MA_meth/data/"

aineisto <- read.csv(paste0(datafolder, "results_analysis_DMRS.csv"), header = T)
#cytosines.small <- fread(paste0(datafolder, "cytosines_nanopore_sample.txt"), header = T)
cytosines.small <- fread(paste0(datafolder, "cytosines_nanopore_small.txt"), header = T)

cytosines.all <- fread(paste0(datafolder, "methylome_cytosines_all_nanopore.txt")) #Load the data
#Need to filter for only the seven chromosomes
cytosines.all <- filter(cytosines.all, seqnames %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Filtering for centromeric cytosines only
cytosines.cent <- filter(cytosines.all, (seqnames == "Supercontig_12.1" & start >= 3678400 & start <= 3988600) | (seqnames == "Supercontig_12.2" & start >= 1091200 & start <= 1368400) | (seqnames == "Supercontig_12.3" & start >= 682000 & start <= 968000) | (seqnames == "Supercontig_12.4" & start >= 869000 & start <= 1091200) | (seqnames == "Supercontig_12.5" & start >= 888800 & start <= 1229800) | (seqnames == "Supercontig_12.6" & start >= 2783000 & start <= 3080000) | (seqnames == "Supercontig_12.7" & start >= 2065800 & start <= 2070200) | (seqnames == "Supercontig_12.7" & start >= 2074600 & start <= 2105400) | (seqnames == "Supercontig_12.7" & start >= 2123000 & start <= 2125200) | (seqnames == "Supercontig_12.7" & start >= 2131800 & start <= 2508000))  

aineisto.cytosines <- intervals.cytosine.check(aineisto, cytosines.all)

#Writing the results to a file
write.table(aineisto.cytosines, file = paste0(datafolder, "cytosines.intervals.csv") sep = ",", row.names = F)

   
### **** Analysis of cytosine methylation changes and genetic mutations

datafolder <- "~/Documents/tutkijatohtori/epimutation/MA_meth/data/"

aineisto.cytosines <- read.csv(paste0(datafolder, "cytosines.intervals.csv"), header = T)

#nomut <- filter(aineisto, Num_mutation == 0)
aineisto.cytosines$interval <- rep(c("Anc--T1", "T1--T5", "T5--T7", "T7--T8", "T8--T10", "T10--T15"), 6)
aineisto.cytosines$Line <- c(rep("Line 11", 6), rep("Line 25", 6), rep("Line 31", 6), rep("Line 5", 6), rep("Line 2", 6), rep("Line 23", 6))

#Making an annotation label
aineisto.cytosines$lab <- paste("g =", aineisto.cytosines$Num_mutation, sep = " ")


#Need to change the order of factors
aineisto.cytosines$interval <- factor(aineisto.cytosines$interval, levels = c("Anc--T1", "T1--T5", "T5--T7", "T7--T8", "T8--T10", "T10--T15"))
aineisto.cytosines$Line <- factor(aineisto.cytosines$Line, levels = c("Line 2", "Line 5", "Line 11", "Line 23", "Line 25", "Line 31"))

#Change into a long format for plotting
aineisto.cytosines.long <- pivot_longer(aineisto.cytosines, cols = c(7:9), names_to = "type", values_to = "epimutation")
#Add the type variable for plotting purposes
aineisto.cytosines$type <- NA

#Drawing some boxes to highlight transfers with zero mutations
library(grid)
rect <- rectGrob(
  x = unit(1, "in"),
  y = unit(1, "npc") - unit(1, "in"),
  width = unit(1, "in"),
  height = unit(1, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "skyblue2", alpha = 0.5)
)

#Making a plot with genetic mutations and methylation changes
gemut.cyt <- ggplot(aineisto.cytosines.long, aes(y = epimutation, x = interval, fill = type)) +
    geom_bar(stat = "identity", colour = "black", position = position_dodge()) +
    scale_fill_manual(values = c(Gains = "blue", Losses = "red", Total = "grey")) +
    scale_y_continuous(expand = c(0,0), breaks = c(seq(0,1e05, by = 2.5e04)), limits = c(0, 1.3e05)) +
    ylab("Number of cytosine methylation changes") +
    xlab("Transfer interval") +
    #annotate(geom = "text", x = 3, y = 5e05, label = "Some text") +
    geom_text(data = aineisto.cytosines, mapping = aes(x = interval, y = -Inf, label = lab), hjust = -0.1, vjust = -8) +
    facet_grid(Line ~ .) +
    theme(legend.position = "top", legend.justification = "center", legend.title = element_blank())

#Figure with highlighted boxes
gemut.cyt.final <- ggdraw(gemut.cyt) + 
  geom_rect(
    aes(xmin = 0.545, xmax = 0.67, ymin = 0.8, ymax = 0.94),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.41, xmax = 0.67, ymin = 0.65, ymax = 0.79),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.545, xmax = 0.67, ymin = 0.51, ymax = 0.645),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.68, xmax = 0.805, ymin = 0.36, ymax = 0.5),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.545, xmax = 0.67, ymin = 0.215, ymax = 0.355),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.14, xmax = 0.265, ymin = 0.215, ymax = 0.355),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  )  +
  geom_rect( #White box the cover the extra NA label... this is a hack... : (
    aes(xmin = 0.67, xmax = 0.75, ymin = 0.95, ymax = 0.99),
    fill = "white", alpha = 1,
    inherit.aes = FALSE
  )

#Saving the plot
save_plot(file = "./epigenet_ms/fig/genmut.pdf", gemut.cyt, base_height = 8)

save_plot(file = "./epigenet_ms/fig/genmut2.pdf", gemut.cyt.final, base_height = 7.8, base_width = 8)

#Fit regression of cytosine methylation changes against number of genetic mutations
summary(lm(Total ~ Num_mutation, data = aineisto.cytosines)) #p = 0.161, R2 = 0.03
summary(lm(Gains ~ Num_mutation, data = aineisto.cytosines)) #p = 0.346, R2 = 0
summary(lm(Losses ~ Num_mutation, data = aineisto.cytosines)) #p = 0.183, R2 = 0.02

ggplot(aineisto.cytosines, aes(y = Total, x = Num_mutation)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab("Number of genetic mutations") +
    ylab("Number of cytosine methylation changes")    

numbers <- ggplot(aineisto.cytosines.long, aes(y = epimutation, x = Num_mutation)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab("Number of genetic mutations") +
    ylab("Number of cytosine methylation changes") +
    facet_wrap( ~ type)

#Drawing the final plot with labels
numbers.final <- ggdraw(numbers) + draw_label("p = 0.346", x = 0.25, y = 0.55) + draw_label("p = 0.183", x = 0.55, y = 0.55) + draw_label("p = 0.161", x = 0.85, y = 0.35)

#Saving it
save_plot(file = "./epigenet_ms/fig/cytosines_genmut_numbers.pdf", numbers.final, base_height = 4, base_width = 8)


### *** DMR changes and genetic mutations

##Mariana had checked which genetic mutations were present in the MA lines sequenced with nanopore at different time points

datafolder <- "~/Documents/tutkijatohtori/epimutation/MA_meth/data/"

aineisto <- read.csv(paste0(datafolder, "results_analysis_DMRS.csv"), header = T)

#nomut <- filter(aineisto, Num_mutation == 0)
aineisto$interval <- rep(c("Anc--T1", "T1--T5", "T5--T7", "T7--T8", "T8--T10", "T10--T15"), 6)
#Need to have matching names in the genetic mutation data and DMRs
aineisto$Line <- fct_recode(aineisto$Line, "L02" = "L2")

#Load DMRs from the nanopore dataset
load(file = "./data/DMRs_all_nanopore.RData") #In DMRs.all.nanopore

#Getting all DMR that occurred in particular intervals
All.DMR.gen <- intervals.DMR.check(aineisto, DMRs.all.nanopore, ind.cols = 5:42)

#Filtering only for centromeric DMRs
DMRs.centromeric.nanopore <- filter(DMRs.all.nanopore, domain == "Centromeric")
cent.DMR.gen <-  intervals.DMR.check(aineisto, DMRs.centromeric.nanopore, ind.cols = 5:42)

#aineisto$Line <- c(rep("Line 11", 6), rep("Line 25", 6), rep("Line 31", 6), rep("Line 5", 6), rep("Line 2", 6), rep("Line 23", 6))

#Combining datasets for plotting
All.DMR.gen$domain <- "All"
cent.DMR.gen$domain <- "Centromeric"
All.DMR.gen$interval <- rep(c("Anc--T1", "T1--T5", "T5--T7", "T7--T8", "T8--T10", "T10--T15"), 6)
cent.DMR.gen$interval <- rep(c("Anc--T1", "T1--T5", "T5--T7", "T7--T8", "T8--T10", "T10--T15"), 6)
All.DMR.gen$Line <- c(rep("Line 11", 6), rep("Line 25", 6), rep("Line 31", 6), rep("Line 5", 6), rep("Line 2", 6), rep("Line 23", 6))
cent.DMR.gen$Line <- c(rep("Line 11", 6), rep("Line 25", 6), rep("Line 31", 6), rep("Line 5", 6), rep("Line 2", 6), rep("Line 23", 6))

#combine and make a longer version for plotting
DMR.gen <- rbind(All.DMR.gen, cent.DMR.gen)
aineisto.DMR.long <- pivot_longer(DMR.gen, cols = c(8:10), names_to = "type", values_to = "epimutation")

#Making an annotation label
All.DMR.gen$lab <- paste("g =", All.DMR.gen$Num_mutation, sep = " ")
All.DMR.gen$type <- NA #For plotting
All.DMR.gen$interval <- factor(All.DMR.gen$interval, levels = c("Anc--T1", "T1--T5", "T5--T7", "T7--T8", "T8--T10", "T10--T15"))
All.DMR.gen$Line <- factor(All.DMR.gen$Line, levels = c("Line 2", "Line 5", "Line 11", "Line 23", "Line 25", "Line 31"))

#Need to change the order of factors
aineisto.DMR.long$interval <- factor(aineisto.DMR.long$interval, levels = c("Anc--T1", "T1--T5", "T5--T7", "T7--T8", "T8--T10", "T10--T15"))
aineisto.DMR.long$Line <- factor(aineisto.DMR.long$Line, levels = c("Line 2", "Line 5", "Line 11", "Line 23", "Line 25", "Line 31"))

#Making the plot
#Making a plot with genetic mutations and methylation changes
gemut.DMR <- ggplot(aineisto.DMR.long, aes(y = epimutation, x = domain, fill = type)) +
    geom_bar(stat = "identity", colour = "black", position = position_dodge()) +
    scale_fill_manual(values = c(Gains = "blue", Losses = "red", Total = "grey")) +
    scale_y_continuous(expand = c(0,0), breaks = c(seq(0,4000, by = 2000)), limits = c(0, 6000)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +    
    ylab("Number of DMR changes") +
    xlab("") +
    #annotate(geom = "text", x = 3, y = 5e05, label = "Some text") +
    geom_text(data = All.DMR.gen, mapping = aes(x = domain, y = -Inf, label = lab), hjust = -1.0, vjust = -6.5) +
    facet_grid(Line ~ interval) +
    theme(legend.position = "top", legend.justification = "center", legend.title = element_blank())

#Figure with highlighted boxes
gemut.DMR.final <- ggdraw(gemut.DMR) + 
  geom_rect(
    aes(xmin = 0.535, xmax = 0.665, ymin = 0.785, ymax = 0.90),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.39, xmax = 0.665, ymin = 0.655, ymax = 0.78),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.535, xmax = 0.665, ymin = 0.525, ymax = 0.650),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.68, xmax = 0.810, ymin = 0.395, ymax = 0.525),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.535, xmax = 0.665, ymin = 0.270, ymax = 0.390),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = 0.105, xmax = 0.225, ymin = 0.270, ymax = 0.390),
    fill = "skyblue", alpha = 0.25,
    inherit.aes = FALSE
  )  +
  geom_rect( #White box the cover the extra NA label... this is a hack... : (
    aes(xmin = 0.67, xmax = 0.75, ymin = 0.95, ymax = 0.99),
    fill = "white", alpha = 1,
    inherit.aes = FALSE
  )


save_plot(file = "./epigenet_ms/fig/DMR_genmut.pdf", gemut.DMR.final, base_height = 7.8, base_width = 8)

#Calculating DMR changes and number of genetic mutations per transfer
#transfers per interval
All.DMR.gen$trans.int <- rep(c(1,4,2,1,2,5),6)
cent.DMR.gen$trans.int <- rep(c(1,4,2,1,2,5),6)
aineisto.DMR.long$trans.int <- rep(rep(c(1,1,1,4,4,4,2,2,2,1,1,1,2,2,2,5,5,5),6),2)


#Fit a regression of DMR changes against the number of genetic mutations
summary(lm(Total ~ Num_mutation, data = All.DMR.gen)) #p = 0.392, R2 = 0
summary(lm(Gains ~ Num_mutation, data = All.DMR.gen)) #p = 0.326, R2 = 0
summary(lm(Losses ~ Num_mutation, data = All.DMR.gen)) #p = 0.0505, R2 = 0.08

summary(lm(Total ~ Num_mutation, data = cent.DMR.gen)) #p = 0.194, R2= 0.02
summary(lm(Gains ~ Num_mutation, data = cent.DMR.gen)) #p = 0.579, R2 = 0
summary(lm(Losses ~ Num_mutation, data = cent.DMR.gen)) #p = 0.0466 R2 = 0.09

#Number of changes for with and without mutations
DMR.numbers <- ggplot(aineisto.DMR.long, aes(y = epimutation, x = Mutation)) +
    geom_boxplot() +
    geom_smooth(method = "lm") +
    xlab("Genetic mutations occurred?") +
    ylab("Number of DMR changes") +
    facet_grid(type ~ domain)

summary(lm(Total ~ Mutation, data = All.DMR.gen)) #p = 0.437, R2 = 0
summary(lm(Gains ~ Mutation, data = All.DMR.gen)) #p = 0.222, R2 = 0.02
summary(lm(Losses ~ Mutation, data = All.DMR.gen)) #p = 0.955, R2 = 0

summary(lm(Total ~ Mutation, data = cent.DMR.gen)) #p = 0.662, R2 = 0
summary(lm(Gains ~ Mutation, data = cent.DMR.gen)) #p = 0.340, R2 = 0
summary(lm(Losses ~ Mutation, data = cent.DMR.gen)) #p = 0.798, R2 = 0

#Drawing the final plot with labels
DMR.number.final <- ggdraw(DMR.numbers) + draw_label("p = 0.222", x = 0.35, y = 0.90) + draw_label("p = 0.955", x = 0.35, y = 0.625) + draw_label("p = 0.437", x = 0.35, y = 0.125) + draw_label("p = 0.340", x = 0.75, y = 0.90) + draw_label("p = 0.798", x = 0.75, y = 0.625) + draw_label("p = 0.662", x = 0.75, y = 0.325)

##Consider only intervals with a single transfer, so that DMR number is accurate
single.tra.long <- filter(aineisto.DMR.long, interval == "Anc--T1" | interval == "T7--T8")
single.tra.all <- filter(All.DMR.gen, interval == "Anc--T1" | interval == "T7--T8")
single.tra.cent <- filter(cent.DMR.gen, interval == "Anc--T1" | interval == "T7--T8")

DMR.number.single <- ggplot(single.tra.long, aes(y = epimutation, x = Mutation)) +
    geom_boxplot() +
    geom_smooth(method = "lm") +
    xlab("Genetic mutations occurred?") +
    ylab("Number of DMR changes") +
    facet_grid(type ~ domain)

summary(lm(Total ~ Mutation, data = single.tra.all)) #p = 0.0577
summary(lm(Gains ~ Mutation, data = single.tra.all)) #p = 0.143
summary(lm(Losses ~ Mutation, data = single.tra.all)) #p = 0.0541

summary(lm(Total ~ Mutation, data = single.tra.cent)) #p = 0.062
summary(lm(Gains ~ Mutation, data = single.tra.cent)) #p = 0.218
summary(lm(Losses ~ Mutation, data = single.tra.cent)) #p = 0.0956

DMR.number.single.final <- ggdraw(DMR.number.single) + draw_label("p = 0.143", x = 0.35, y = 0.90) + draw_label("p = 0.054", x = 0.35, y = 0.625) + draw_label("p = 0.058", x = 0.35, y = 0.125) + draw_label("p = 0.218", x = 0.75, y = 0.90) + draw_label("p = 0.096", x = 0.75, y = 0.625) + draw_label("p = 0.062", x = 0.75, y = 0.325)

plot.DMR.numbers <- plot_grid(DMR.number.final, DMR.number.single.final, labels = c("A", "B"),  ncol = 2)

save_plot(filename = "./epigenet_ms/fig/DMR_numbers_gen.pdf", plot.DMR.numbers, base_height = 8, base_width = 8)

#Number of changes per transfer
ggplot(aineisto.DMR.long, aes(y = epimutation/trans.int, x = Num_mutation/trans.int)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab("Number of genetic mutations") +
    ylab("Number of DMR changes") +
    facet_grid(type ~ domain)

#Number of changes and number of mutations
ggplot(aineisto.DMR.long, aes(y = epimutation, x = Num_mutation)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab("Number of genetic mutations") +
    ylab("Number of DMR changes") +
    facet_grid(type ~ domain)



#Saving it
save_plot(file = "./epigenet_ms/fig/cytosines_genmut_numbers.pdf", numbers.final, base_height = 4, base_width = 8)

### * Exploring divergence outliers

### ** Identifying the outliers

### Load the data

load(file = "./pedigree_data/all/MA_CG_all.RData")
load(file = "./pedigree_data/all/MA_CHG_all.RData")
load(file = "./pedigree_data/all/MA_CHH_all.RData")

#Need to do some hacks to get the samples for each pairwise comparison
#genTable <- fread("nodelist.fn")
#genTable <- genTable %>% filter(genTable$meth == "Y")
#pairs <- combn(genTable$filename, 2)
#pairs <- t(pairs) #transpose
#fix the sample names
#pairs[,1] <- gsub("/mnt/nas/ilkka/results/methimpute-out/methylome_", "", pairs[,1])
#pairs[,2] <- gsub("/mnt/nas/ilkka/results/methimpute-out/methylome_", "", pairs[,2])
#pairs[,1] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,1])
#pairs[,2] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,2])
##Okay, problem is that order of the samples in the pedigree divergence data and sample pairs is not the same


#Maybe load matA and mata pedigrees separately and figure out the order?
#mat A samples
genTable <- fread("nodelist_matA.fn")
genTable <- genTable %>% filter(genTable$meth == "Y")
pairs <- combn(genTable$filename, 2)
pairs <- t(pairs) #transpose
#fix the sample names
pairs[,1] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,1])
pairs[,2] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,2])
pairs[,1] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,1])
pairs[,2] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,2])
pairs.matA <- pairs

#mat a samples
genTable <- fread("nodelist_mata.fn")
genTable <- genTable %>% filter(genTable$meth == "Y")
pairs <- combn(genTable$filename, 2)
pairs <- t(pairs) #transpose
pairs[,1] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,1])
pairs[,2] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,2])
pairs[,1] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,1])
pairs[,2] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,2])
pairs.mata <- pairs

load(file = "./pedigree_data/MA_pedigree_CG_matA.RData")

#Final pedigree CG matA
pedigree.matA.CG <- output.matA$Pdata
pedigree.matA.CG[,1:3] <- pedigree.matA.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CG[,2] + pedigree.matA.CG[,3] - 2*pedigree.matA.CG[,1]
p0uu_in.CG <- output.matA$tmpp0
pedigree.matA.CG <- data.frame(pedigree.matA.CG, pairs.matA, dt)
colnames(pedigree.matA.CG)[5:6] <- c("sample1", "sample2")

load(file = "./pedigree_data/MA_pedigree_CG_mata.RData")

#Final pedigree CG mat a
pedigree.mata.CG <- output.CG$Pdata
pedigree.mata.CG[,1:3] <- pedigree.mata.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CG[,2] + pedigree.mata.CG[,3] - 2*pedigree.mata.CG[,1]
p0uu_in.CG <- output.CG$tmpp0
pedigree.mata.CG <- data.frame(pedigree.mata.CG, pairs.mata, dt)
colnames(pedigree.mata.CG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CG, pedigree.mata.CG)

#Plot of combined pedigrees divergence
myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
plot.combined.pedigrees <- ggplot(pedigree.combined, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits, aes(x = dt, y = divsim), col = "red" ) +
    #geom_smooth() +
    #geom_line(aes(x = dt, y = pred), col = "blue") +
    #geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +
    geom_hline(yintercept = 0.00175, lty = "dashed") +
    scale_y_continuous(limits = c(0.0008, 0.0026)) +
    ggtitle("All data") +    
    xlab(myxlab) +
    ylab("Methylation divergence")

#Check that loading all of the data directly gives the same output

#Check that pedigree is correct
#plotPedigree(nodelist = "nodelist.fn", edgelist = "edgelist.fn",
#sampling.design = "sibling", plot.width = 15, plot.height = 15,
#aspect.ratio = 2.5, vertex.size = 12, vertex.label = FALSE) #Should be OK

#Final pedigree CG
pedigree.CG <- output.all.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- output.all.CG$tmpp0
#pedigree.CG <- data.frame(pedigree.CG, pairs)
#colnames(pedigree.CG)[5:6] <- c("sample1", "sample2")
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]
pedigree.final.CG <- data.frame(pedigree.final.CG, dt = dt.final)

#Making a plot
plot.normal.CG.pedi <- ggplot(pedigree.final.CG, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits, aes(x = dt, y = divsim), col = "red" ) +
    #geom_smooth() +
    #geom_line(aes(x = dt, y = pred), col = "blue") +
    #geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +
    geom_hline(yintercept = 0.0015, lty = "dashed") +
    ggtitle("Single sites") +    
    xlab(myxlab) +
    ylab("Methylation divergence")

plot_grid(plot.combined.pedigrees, plot.normal.CG.pedi) #These are identical as they should be

### Figuring out which samples are the outliers

#Using the pedigree.combined data
#
pedigree.combined$comp <- paste(pedigree.combined$sample1, pedigree.combined$sample2, sep = "_")

#Using an arbitary threshold of 0.00175

#Number of different samples
samples <- unique(c(as.character(pedigree.combined$sample1), as.character(pedigree.combined$sample2)))
n.samples <- length(samples) #69 samples in total

#Counting the number of times a sample is involved in an outlier observation
outlier.res <- data.frame(sample = samples, outlier.count = 0)
for(i in 1:n.samples) {
    current.sample <- samples[i]
    #Filter for those cases where the current sample is either one of the two samples being compared
    current <- filter(pedigree.combined, sample1 == current.sample | sample2 == current.sample)
    outlier.res[i,2] <- sum(current$D.value > 0.00175) #Count number of outlier points
}

#Okay, this shows that there are four samples that have that seems to cause the majority of the outlier points: MG20L13, ML6G40, MG20L35, ML36G40


outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.outl.rm <- pedigree.combined[!outl.samp,]

##Highlight pairwise comparisons that involved an outlier sample
pedigree.combined$highl <- "other"
pedigree.combined$highl[(pedigree.combined$sample1 == "MG20L13" | pedigree.combined$sample2 == "MG20L13")] <- "L13G20"
pedigree.combined$highl[pedigree.combined$sample1 == "ML6G40" | pedigree.combined$sample2 == "ML6G40"] <- "L6G40"
pedigree.combined$highl[pedigree.combined$sample1 == "MG20L35" | pedigree.combined$sample2 == "MG20L35"] <- "L35G20"
pedigree.combined$highl[pedigree.combined$sample1 == "ML36G40" | pedigree.combined$sample2 == "ML36G40"] <- "L36G40"

plot.outl.high <- ggplot(pedigree.combined, aes(x = dt, y = D.value, colour = highl)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits, aes(x = dt, y = divsim), col = "red" ) +
    #geom_smooth() +
    #geom_line(aes(x = dt, y = pred), col = "blue") +
    #geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +
    geom_hline(yintercept = 0.00175, lty = "dashed") +
    scale_y_continuous(limits = c(0.0008, 0.0026)) +
    scale_colour_manual(values = c(other = "black", L13G20 = "red", L6G40 = "blue", L35G20 = "purple", L36G40 = "orange")) +
    ggtitle("Outlier samples highlighted") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    theme(legend.title = element_blank())

##Divergence without outlier samples
plot.outl.rm <- ggplot(pedigree.outl.rm, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits, aes(x = dt, y = divsim), col = "red" ) +
    #geom_smooth() +
    #geom_line(aes(x = dt, y = pred), col = "blue") +
    #geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +
    geom_hline(yintercept = 0.00175, lty = "dashed") +
    scale_y_continuous(limits = c(0.0008, 0.0026)) +    
    ggtitle("Outlier samples removed") +    
    xlab(myxlab) +
    ylab("Methylation divergence")

#Making the final outlier plot
outlier.plot <- plot_grid(plot.combined.pedigrees, plot.outl.high, plot.outl.rm, ncol = 3, rel_widths = c(1,1.3,1))
save_plot("./epigenet_ms/fig/outliers_ss.pdf",outlier.plot, base_width = 14, base_height = 4.25)


### Make a final dataset without outliers


#Final pedigree CHG
pedigree.CHG <- output.all.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- output.all.CHG$tmpp0
pedigree.CHG <- data.frame(pedigree.CHG, pairs)
colnames(pedigree.CHG)[5:6] <- c("sample1", "sample2")
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- output.all.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- output.all.CHH$tmpp0
pedigree.CHH <- data.frame(pedigree.CHH, pairs)
colnames(pedigree.CHH)[5:6] <- c("sample1", "sample2")
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
methsites.plot.outliers <- ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits, aes(x = dt, y = divsim), col = "red" ) +
    #geom_smooth() +
    #geom_line(aes(x = dt, y = pred), col = "blue") +
    #geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +
    geom_hline(data = data.frame(thr = c(0.0015, 0.0010, 0.00125), context = c("CG", "CHG", "CHH")), aes(yintercept = thr), lty = "dashed") +
    ggtitle("Single sites") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

#Divergence values larger than 0.0015 seem to be outliers
#What is the number of unique samples in pairwise sample combinations among the outlier points?
#So are all of the outlier points caused by one or two samples?
outliers <- filter(plotdata, (D.value > 0.0015 & context == "CG") | (D.value > 0.0010 & context == "CHG") | (D.value > 0.00125 & context == "CHH"))

#There are strains that particularly occur among the outlier points
table(outliers$sample1)
table(outliers$sample2)
#Counting the number of cases where both samples are different
outliers$sample1

 ggplot(outliers, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits, aes(x = dt, y = divsim), col = "red" ) +
    #geom_smooth() +
    #geom_line(aes(x = dt, y = pred), col = "blue") +
    #geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +
    geom_hline(data = data.frame(thr = c(0.0015, 0.0010, 0.00125), context = c("CG", "CHG", "CHH")), aes(yintercept = thr), lty = "dashed") +
    ggtitle("Single sites") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

##Plotting
plotdata <- cbind(pedigree.final.CG, dt.final)

ggplot(plotdata, aes(x = dt.final, y = D.value)) +
    geom_point(alpha = 0.2)





buildPedigree
function (nodelist, edgelist, cytosine = "CG", posteriorMaxFilter = 0.99) 
{
    mt <- startTime("constracting pedigree ...\n")
    dmatrix <- dMatrix(nodelist, cytosine, posteriorMaxFilter)
    message("generating Rc.Meth.lvl from sample file....\n")
    rclvl <- rc.meth.lvl(nodelist, cytosine, posteriorMaxFilter)
    props <- rclvl[which(as.character(rclvl[, 2]) == cytosine), 
        ]
    outliers <- "none"
    props <- rclvl[which(!is.element(rclvl[, 1], outliers) == 
        TRUE), ]
    tmpP0uu <- 1 - mean(as.numeric(as.character(props[, 3])))
    message("finilizing pedegree data...")
    edges <- fread(edgelist, header = TRUE, skip = 0)
    samples <- fread(nodelist, header = TRUE, skip = 0, select = c(2, 
        3, 4))
    if (!is.null(samples$Branchpoint_date)) {
        colnames(samples)[2] <- "gen"
    }
    tmpPedegree <- convertDMATRIX(samples, edges, dmatrix)
    cat(stopTime(mt))
    return(list(Pdata = tmpPedegree, tmpp0 = tmpP0uu))
}

buildPedigree(nodelist = "nodelist.fn", edgelist = "edgelist.fn", cytosine = "CG", posteriorMaxFilter = 0.99)

### ** Analysis of divergence without outlier samples

### *** Produce datasets that contain sample information and remove outliers

#Samples involved in pairwise comparisons

#mat A samples
genTable <- fread("nodelist_matA.fn")
genTable <- genTable %>% filter(genTable$meth == "Y")
pairs <- combn(genTable$filename, 2)
pairs <- t(pairs) #transpose
#fix the sample names
pairs[,1] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,1])
pairs[,2] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,2])
pairs[,1] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,1])
pairs[,2] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,2])
pairs.matA <- pairs

#mat a samples
genTable <- fread("nodelist_mata.fn")
genTable <- genTable %>% filter(genTable$meth == "Y")
pairs <- combn(genTable$filename, 2)
pairs <- t(pairs) #transpose
pairs[,1] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,1])
pairs[,2] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,2])
pairs[,1] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,1])
pairs[,2] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,2])
pairs.mata <- pairs

### **** Single cytosines, whole genome
load(file = "./pedigree_data/MA_pedigree_CG_matA.RData") #output.matA
load(file = "./pedigree_data/MA_pedigree_CHG_matA.RData") #output.matA.CHG
load(file = "./pedigree_data/MA_pedigree_CHH_matA.RData") #output.matA.CHH

load(file = "./pedigree_data/MA_pedigree_CG_mata.RData") #output.CG
load(file = "./pedigree_data/MA_pedigree_CHG_mata.RData") #output.CHG
load(file = "./pedigree_data/MA_pedigree_CHH_mata.RData") #output.mata.CHH

### CG ###
#Final pedigree CG matA
pedigree.matA.CG <- output.matA$Pdata
pedigree.matA.CG[,1:3] <- pedigree.matA.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CG[,2] + pedigree.matA.CG[,3] - 2*pedigree.matA.CG[,1]
p0uu_in.CG <- output.matA$tmpp0
pedigree.matA.CG <- data.frame(pedigree.matA.CG, pairs.matA, dt)
colnames(pedigree.matA.CG)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CG <- output.CG$Pdata
pedigree.mata.CG[,1:3] <- pedigree.mata.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CG[,2] + pedigree.mata.CG[,3] - 2*pedigree.mata.CG[,1]
p0uu_in.CG <- output.CG$tmpp0
pedigree.mata.CG <- data.frame(pedigree.mata.CG, pairs.mata, dt)
colnames(pedigree.mata.CG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CG, pedigree.mata.CG)

outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.all.ss.CG.outl.rm <- pedigree.combined[!outl.samp,]

save(pedigree.all.ss.CG.outl.rm, p0uu_in.CG, file = "./pedigree_data/all/MA_pedigree_CG_outl_rm.RData")
###########

### CHG ###
#Final pedigree CHG matA
pedigree.matA.CHG <- output.matA.CHG$Pdata
pedigree.matA.CHG[,1:3] <- pedigree.matA.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHG[,2] + pedigree.matA.CHG[,3] - 2*pedigree.matA.CHG[,1]
p0uu_in.CHG <- output.matA.CHG$tmpp0
pedigree.matA.CHG <- data.frame(pedigree.matA.CHG, pairs.matA, dt)
colnames(pedigree.matA.CHG)[5:6] <- c("sample1", "sample2")

#Final pedigree CHG mat a
pedigree.mata.CHG <- output.CHG$Pdata
pedigree.mata.CHG[,1:3] <- pedigree.mata.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHG[,2] + pedigree.mata.CHG[,3] - 2*pedigree.mata.CHG[,1]
p0uu_in.CHG <- output.CHG$tmpp0
pedigree.mata.CHG <- data.frame(pedigree.mata.CHG, pairs.mata, dt)
colnames(pedigree.mata.CHG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHG, pedigree.mata.CHG)

outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.all.ss.CHG.outl.rm <- pedigree.combined[!outl.samp,]

save(pedigree.all.ss.CHG.outl.rm, p0uu_in.CHG, file = "./pedigree_data/all/MA_pedigree_CHG_outl_rm.RData")
##########

### CHH ###
#Final pedigree CHH matA
pedigree.matA.CHH <- output.matA.CHH$Pdata
pedigree.matA.CHH[,1:3] <- pedigree.matA.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHH[,2] + pedigree.matA.CHH[,3] - 2*pedigree.matA.CHH[,1]
p0uu_in.CHH <- output.matA.CHH$tmpp0
pedigree.matA.CHH <- data.frame(pedigree.matA.CHH, pairs.matA, dt)
colnames(pedigree.matA.CHH)[5:6] <- c("sample1", "sample2")

#Final pedigree CHH mat a
pedigree.mata.CHH <- output.mata.CHH$Pdata
pedigree.mata.CHH[,1:3] <- pedigree.mata.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHH[,2] + pedigree.mata.CHH[,3] - 2*pedigree.mata.CHH[,1]
p0uu_in.CHH <- output.mata.CHH$tmpp0
pedigree.mata.CHH <- data.frame(pedigree.mata.CHH, pairs.mata, dt)
colnames(pedigree.mata.CHH)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHH, pedigree.mata.CHH)

outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.all.ss.CHH.outl.rm <- pedigree.combined[!outl.samp,]

save(pedigree.all.ss.CHH.outl.rm, p0uu_in.CHH, file = "./pedigree_data/all/MA_pedigree_CHH_outl_rm.RData")

#Check by plotting
#Plotting
#Make plot dataframe
plotdata <- data.frame(rbind(pedigree.all.ss.CG.outl.rm, pedigree.all.ss.CHG.outl.rm, pedigree.all.ss.CHH.outl.rm), context = c(rep("CG", nrow(pedigree.all.ss.CG.outl.rm)), rep("CHG", nrow(pedigree.all.ss.CHG.outl.rm)), rep("CHH", nrow(pedigree.all.ss.CHH.outl.rm))))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("Single sites") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)


### **** DMRs, whole genome

load(file = "./pedigree_data/DMR_all/MA_DMR_all_mat_CG.RData") #output.matA.DMR.all.CG, output.mata.DMR.all.CG
load(file = "./pedigree_data/DMR_all/MA_DMR_all_mat_CHG.RData") #output.matA.DMR.all.CHG ,output.mata.DMR.all.CHG
load(file = "./pedigree_data/DMR_all/MA_DMR_all_mat_CHH.RData") #output.matA.DMR.all.CHH, output.mata.DMR.all.CHH

### CG ###
#Final pedigree CG matA
pedigree.matA.CG <- output.matA.DMR.all.CG$Pdata
pedigree.matA.CG[,1:3] <- pedigree.matA.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CG[,2] + pedigree.matA.CG[,3] - 2*pedigree.matA.CG[,1]
p0uu_in.CG <- output.matA.DMR.all.CG$tmpp0
pedigree.matA.CG <- data.frame(pedigree.matA.CG, pairs.matA, dt)
colnames(pedigree.matA.CG)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CG <- output.mata.DMR.all.CG$Pdata
pedigree.mata.CG[,1:3] <- pedigree.mata.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CG[,2] + pedigree.mata.CG[,3] - 2*pedigree.mata.CG[,1]
p0uu_in.CG <- output.mata.DMR.all.CG$tmpp0
pedigree.mata.CG <- data.frame(pedigree.mata.CG, pairs.mata, dt)
colnames(pedigree.mata.CG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CG, pedigree.mata.CG)

outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.all.DMR.CG.outl.rm <- pedigree.combined[!outl.samp,]

save(pedigree.all.DMR.CG.outl.rm, p0uu_in.CG, file = "./pedigree_data/DMR_all/MA_DMR_all_CG_outl_rm.RData")
###########

### CHG ###
#Final pedigree CHG matA
pedigree.matA.CHG <- output.matA.DMR.all.CHG$Pdata
pedigree.matA.CHG[,1:3] <- pedigree.matA.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHG[,2] + pedigree.matA.CHG[,3] - 2*pedigree.matA.CHG[,1]
p0uu_in.CHG <- output.matA.DMR.all.CHG$tmpp0
pedigree.matA.CHG <- data.frame(pedigree.matA.CHG, pairs.matA, dt)
colnames(pedigree.matA.CHG)[5:6] <- c("sample1", "sample2")

#Final pedigree CHG mat a
pedigree.mata.CHG <- output.mata.DMR.all.CHG$Pdata
pedigree.mata.CHG[,1:3] <- pedigree.mata.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHG[,2] + pedigree.mata.CHG[,3] - 2*pedigree.mata.CHG[,1]
p0uu_in.CHG <- output.mata.DMR.all.CHG$tmpp0
pedigree.mata.CHG <- data.frame(pedigree.mata.CHG, pairs.mata, dt)
colnames(pedigree.mata.CHG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHG, pedigree.mata.CHG)

outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.all.DMR.CHG.outl.rm <- pedigree.combined[!outl.samp,]

save(pedigree.all.DMR.CHG.outl.rm, p0uu_in.CHG, file = "./pedigree_data/DMR_all/MA_DMR_all_CHG_outl_rm.RData")
###########

### CHH ###
#Final pedigree CHH matA
pedigree.matA.CHH <- output.matA.DMR.all.CHH$Pdata
pedigree.matA.CHH[,1:3] <- pedigree.matA.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHH[,2] + pedigree.matA.CHH[,3] - 2*pedigree.matA.CHH[,1]
p0uu_in.CHH <- output.matA.DMR.all.CHH$tmpp0
pedigree.matA.CHH <- data.frame(pedigree.matA.CHH, pairs.matA, dt)
colnames(pedigree.matA.CHH)[5:6] <- c("sample1", "sample2")

#Final pedigree CHH mat a
pedigree.mata.CHH <- output.mata.DMR.all.CHH$Pdata
pedigree.mata.CHH[,1:3] <- pedigree.mata.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHH[,2] + pedigree.mata.CHH[,3] - 2*pedigree.mata.CHH[,1]
p0uu_in.CHH <- output.mata.DMR.all.CHH$tmpp0
pedigree.mata.CHH <- data.frame(pedigree.mata.CHH, pairs.mata, dt)
colnames(pedigree.mata.CHH)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHH, pedigree.mata.CHH)

outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.all.DMR.CHH.outl.rm <- pedigree.combined[!outl.samp,]

save(pedigree.all.DMR.CHH.outl.rm, p0uu_in.CHH, file = "./pedigree_data/DMR_all/MA_DMR_all_CHH_outl_rm.RData")

#Check by plotting
#Plotting
#Make plot dataframe
plotdata <- data.frame(rbind(pedigree.all.DMR.CG.outl.rm, pedigree.all.DMR.CHG.outl.rm, pedigree.all.DMR.CHH.outl.rm), context = c(rep("CG", nrow(pedigree.all.DMR.CG.outl.rm)), rep("CHG", nrow(pedigree.all.DMR.CHG.outl.rm)), rep("CHH", nrow(pedigree.all.DMR.CHH.outl.rm))))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

### **** Single cytosines, centromeric H3K9

load(file = "./pedigree_data/centromericH3K9/MA_mat_CG_centromericH3K9.RData") #output.matA.centromericH3K9.CG, output.mata.centromericH3K9.CG
load(file = "./pedigree_data/centromericH3K9/MA_mat_CHG_centromericH3K9.RData") #output.matA.centromericH3K9.CHG, output.mata.centromericH3K9.CHG
load(file = "./pedigree_data/centromericH3K9/MA_mat_CHH_centromericH3K9.RData") #output.matA.CHH

### CG ###
#Final pedigree CG matA
pedigree.matA.CG <- output.matA.centromericH3K9.CG$Pdata
pedigree.matA.CG[,1:3] <- pedigree.matA.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CG[,2] + pedigree.matA.CG[,3] - 2*pedigree.matA.CG[,1]
p0uu_in.CG <- output.matA.centromericH3K9.CG$tmpp0
pedigree.matA.CG <- data.frame(pedigree.matA.CG, pairs.matA, dt)
colnames(pedigree.matA.CG)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CG <- output.mata.centromericH3K9.CG$Pdata
pedigree.mata.CG[,1:3] <- pedigree.mata.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CG[,2] + pedigree.mata.CG[,3] - 2*pedigree.mata.CG[,1]
p0uu_in.CG <- output.mata.centromericH3K9.CG$tmpp0
pedigree.mata.CG <- data.frame(pedigree.mata.CG, pairs.mata, dt)
colnames(pedigree.mata.CG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CG, pedigree.mata.CG)

outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.centromeric.ss.CG.outl.rm <- pedigree.combined[!outl.samp,]

save(pedigree.centromeric.ss.CG.outl.rm, p0uu_in.CG, file = "./pedigree_data/centromericH3K9/MA_centromeric_CG_outl_rm.RData")
###########

### CHG ###
#Final pedigree CHG matA
pedigree.matA.CHG <- output.matA.centromericH3K9.CHG$Pdata
pedigree.matA.CHG[,1:3] <- pedigree.matA.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHG[,2] + pedigree.matA.CHG[,3] - 2*pedigree.matA.CHG[,1]
p0uu_in.CHG <- output.matA.centromericH3K9.CHG$tmpp0
pedigree.matA.CHG <- data.frame(pedigree.matA.CHG, pairs.matA, dt)
colnames(pedigree.matA.CHG)[5:6] <- c("sample1", "sample2")

#Final pedigree CHG mat a
pedigree.mata.CHG <- output.mata.centromericH3K9.CHG$Pdata
pedigree.mata.CHG[,1:3] <- pedigree.mata.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHG[,2] + pedigree.mata.CHG[,3] - 2*pedigree.mata.CHG[,1]
p0uu_in.CHG <- output.mata.centromericH3K9.CHG$tmpp0
pedigree.mata.CHG <- data.frame(pedigree.mata.CHG, pairs.mata, dt)
colnames(pedigree.mata.CHG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHG, pedigree.mata.CHG)

outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.centromeric.ss.CHG.outl.rm <- pedigree.combined[!outl.samp,]

save(pedigree.centromeric.ss.CHG.outl.rm, p0uu_in.CHG, file = "./pedigree_data/centromericH3K9/MA_centromeric_CHG_outl_rm.RData")
##########

### CHH ###
#Final pedigree CHH matA
pedigree.matA.CHH <- output.matA.centromericH3K9.CHH$Pdata
pedigree.matA.CHH[,1:3] <- pedigree.matA.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHH[,2] + pedigree.matA.CHH[,3] - 2*pedigree.matA.CHH[,1]
p0uu_in.CHH <- output.matA.centromericH3K9.CHH$tmpp0
pedigree.matA.CHH <- data.frame(pedigree.matA.CHH, pairs.matA, dt)
colnames(pedigree.matA.CHH)[5:6] <- c("sample1", "sample2")

#Final pedigree CHH mat a
pedigree.mata.CHH <- output.mata.centromericH3K9.CHH$Pdata
pedigree.mata.CHH[,1:3] <- pedigree.mata.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHH[,2] + pedigree.mata.CHH[,3] - 2*pedigree.mata.CHH[,1]
p0uu_in.CHH <- output.mata.centromericH3K9.CHH$tmpp0
pedigree.mata.CHH <- data.frame(pedigree.mata.CHH, pairs.mata, dt)
colnames(pedigree.mata.CHH)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHH, pedigree.mata.CHH)

outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.centromeric.ss.CHH.outl.rm <- pedigree.combined[!outl.samp,]

save(pedigree.centromeric.ss.CHH.outl.rm, p0uu_in.CHH, file = "./pedigree_data/centromericH3K9/MA_centromeric_CHH_outl_rm.RData")


#Check by plotting
#Plotting
#Make plot dataframe
plotdata <- data.frame(rbind(pedigree.centromeric.ss.CG.outl.rm, pedigree.centromeric.ss.CHG.outl.rm, pedigree.centromeric.ss.CHH.outl.rm), context = c(rep("CG", nrow(pedigree.centromeric.ss.CG.outl.rm)), rep("CHG", nrow(pedigree.centromeric.ss.CHG.outl.rm)), rep("CHH", nrow(pedigree.centromeric.ss.CHH.outl.rm))))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("Single sites") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

### **** DMRs, centromeric H3K9

load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_mat_CG.RData") #output.matA.DMR.centromeric.CG, output.mata.DMR.centromeric.CG
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_mat_CHG.RData") #output.matA.DMR.centromeric.CHG ,output.mata.DMR.centromeric.CHG
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_mat_CHH.RData") #output.matA.DMR.centromeric.CHH, output.mata.DMR.centromeric.CHH



### CG ###
#Final pedigree CG matA
pedigree.matA.CG <- output.matA.DMR.centromeric.CG$Pdata
pedigree.matA.CG[,1:3] <- pedigree.matA.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CG[,2] + pedigree.matA.CG[,3] - 2*pedigree.matA.CG[,1]
p0uu_in.CG <- output.matA.DMR.centromeric.CG$tmpp0
pedigree.matA.CG <- data.frame(pedigree.matA.CG, pairs.matA, dt)
colnames(pedigree.matA.CG)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CG <- output.mata.DMR.centromeric.CG$Pdata
pedigree.mata.CG[,1:3] <- pedigree.mata.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CG[,2] + pedigree.mata.CG[,3] - 2*pedigree.mata.CG[,1]
p0uu_in.CG <- output.mata.DMR.centromeric.CG$tmpp0
pedigree.mata.CG <- data.frame(pedigree.mata.CG, pairs.mata, dt)
colnames(pedigree.mata.CG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CG, pedigree.mata.CG)

#Using different outliers samples here
outl.samp <- pedigree.combined$sample1 %in% c("MG5L20", "MG20L1", "MG20L17", "ML1G40", "ML6G40", "ML10G40", "MG5L30", "MG20L29", "MG20L30") | pedigree.combined$sample2 %in% c("MG5L20", "MG20L1", "MG20L17", "ML1G40", "ML6G40", "ML10G40", "MG5L30", "MG20L29", "MG20L30") 
#outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")

#Removing those samples from the data
pedigree.centromeric.DMR.CG.outl.rm2 <- pedigree.combined[!outl.samp,]

save(pedigree.centromeric.DMR.CG.outl.rm2, p0uu_in.CG, file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromeric_CG_outl_rm2.RData")
###########

### CHG ###
#Final pedigree CHG matA
pedigree.matA.CHG <- output.matA.DMR.centromeric.CHG$Pdata
pedigree.matA.CHG[,1:3] <- pedigree.matA.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHG[,2] + pedigree.matA.CHG[,3] - 2*pedigree.matA.CHG[,1]
p0uu_in.CHG <- output.matA.DMR.centromeric.CHG$tmpp0
pedigree.matA.CHG <- data.frame(pedigree.matA.CHG, pairs.matA, dt)
colnames(pedigree.matA.CHG)[5:6] <- c("sample1", "sample2")

#Final pedigree CHG mat a
pedigree.mata.CHG <- output.mata.DMR.centromeric.CHG$Pdata
pedigree.mata.CHG[,1:3] <- pedigree.mata.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHG[,2] + pedigree.mata.CHG[,3] - 2*pedigree.mata.CHG[,1]
p0uu_in.CHG <- output.mata.DMR.centromeric.CHG$tmpp0
pedigree.mata.CHG <- data.frame(pedigree.mata.CHG, pairs.mata, dt)
colnames(pedigree.mata.CHG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHG, pedigree.mata.CHG)

#Using different outliers samples here
outl.samp <- pedigree.combined$sample1 %in% c("MG5L20", "MG20L1", "MG20L17", "ML1G40", "ML6G40", "ML10G40", "MG5L30", "MG20L29", "MG20L30") | pedigree.combined$sample2 %in% c("MG5L20", "MG20L1", "MG20L17", "ML1G40", "ML6G40", "ML10G40", "MG5L30", "MG20L29", "MG20L30") 
#outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.centromeric.DMR.CHG.outl.rm2 <- pedigree.combined[!outl.samp,]

save(pedigree.centromeric.DMR.CHG.outl.rm2, p0uu_in.CHG, file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromeric_CHG_outl_rm2.RData")
###########

### CHH ###
#Final pedigree CHH matA
pedigree.matA.CHH <- output.matA.DMR.centromeric.CHH$Pdata
pedigree.matA.CHH[,1:3] <- pedigree.matA.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHH[,2] + pedigree.matA.CHH[,3] - 2*pedigree.matA.CHH[,1]
p0uu_in.CHH <- output.matA.DMR.centromeric.CHH$tmpp0
pedigree.matA.CHH <- data.frame(pedigree.matA.CHH, pairs.matA, dt)
colnames(pedigree.matA.CHH)[5:6] <- c("sample1", "sample2")

#Final pedigree CHH mat a
pedigree.mata.CHH <- output.mata.DMR.centromeric.CHH$Pdata
pedigree.mata.CHH[,1:3] <- pedigree.mata.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHH[,2] + pedigree.mata.CHH[,3] - 2*pedigree.mata.CHH[,1]
p0uu_in.CHH <- output.mata.DMR.centromeric.CHH$tmpp0
pedigree.mata.CHH <- data.frame(pedigree.mata.CHH, pairs.mata, dt)
colnames(pedigree.mata.CHH)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHH, pedigree.mata.CHH)

#Using different outliers samples here
outl.samp <- pedigree.combined$sample1 %in% c("MG5L20", "MG20L1", "MG20L17", "ML1G40", "ML6G40", "ML10G40", "MG5L30", "MG20L29", "MG20L30") | pedigree.combined$sample2 %in% c("MG5L20", "MG20L1", "MG20L17", "ML1G40", "ML6G40", "ML10G40", "MG5L30", "MG20L29", "MG20L30") 
#outl.samp <- pedigree.combined$sample1 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40") | pedigree.combined$sample2 %in% c("MG20L13", "ML6G40", "MG20L35", "ML36G40")
#Removing those samples from the data
pedigree.centromeric.DMR.CHH.outl.rm2 <- pedigree.combined[!outl.samp,]

save(pedigree.centromeric.DMR.CHH.outl.rm2, p0uu_in.CHH, file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromeric_CHH_outl_rm2.RData")

#Check by plotting
#Plotting
#Make plot dataframe
plotdata <- data.frame(rbind(pedigree.centromeric.DMR.CG.outl.rm, pedigree.centromeric.DMR.CHG.outl.rm, pedigree.centromeric.DMR.CHH.outl.rm), context = c(rep("CG", nrow(pedigree.centromeric.DMR.CG.outl.rm)), rep("CHG", nrow(pedigree.centromeric.DMR.CHG.outl.rm)), rep("CHH", nrow(pedigree.centromeric.DMR.CHH.outl.rm))))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

##With alternative outlier samples
plotdata2 <- data.frame(rbind(pedigree.centromeric.DMR.CG.outl.rm2, pedigree.centromeric.DMR.CHG.outl.rm2, pedigree.centromeric.DMR.CHH.outl.rm2), context = c(rep("CG", nrow(pedigree.centromeric.DMR.CG.outl.rm2)), rep("CHG", nrow(pedigree.centromeric.DMR.CHG.outl.rm2)), rep("CHH", nrow(pedigree.centromeric.DMR.CHH.outl.rm2))))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
ggplot(plotdata2, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)



### **** check centromeric H3K9 outliers

#By looking at plots for divergence in centromeric H3K9 it looks like removing outliers does not change the outlier points

#Are different samples causing the outlier point in centromeric H3K9?

### Identifying outliers for centromericH3K9 DMR samples, looks likely different than for whole genome

#Load the centromeric DMR data

load("./pedigree_data/DMR_centromericH3K9/MA_DMR_centromeric_CG_outl_rm.RData")
load("./pedigree_data/DMR_centromericH3K9/MA_DMR_centromeric_CHG_outl_rm.RData")
load("./pedigree_data/DMR_centromericH3K9/MA_DMR_centromeric_CHH_outl_rm.RData")

plotdata <- data.frame(rbind(pedigree.centromeric.DMR.CG.outl.rm, pedigree.centromeric.DMR.CHG.outl.rm, pedigree.centromeric.DMR.CHH.outl.rm), context = c(rep("CG", nrow(pedigree.centromeric.DMR.CG.outl.rm)), rep("CHG", nrow(pedigree.centromeric.DMR.CHG.outl.rm)), rep("CHH", nrow(pedigree.centromeric.DMR.CHH.outl.rm))))

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

##CHH DMR data shows the clearest outliers

#Loading the unfiltered dataset
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_mat_CG.RData") #output.matA.DMR.centromeric.CG, output.mata.DMR.centromeric.CG
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_mat_CHG.RData") #output.matA.DMR.centromeric.CHG ,output.mata.DMR.centromeric.CHG
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_mat_CHH.RData") #output.matA.DMR.centromeric.CHH, output.mata.DMR.centromeric.CHH

### CHH ###
#Final pedigree CHH matA
pedigree.matA.CHH <- output.matA.DMR.centromeric.CHH$Pdata
pedigree.matA.CHH[,1:3] <- pedigree.matA.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHH[,2] + pedigree.matA.CHH[,3] - 2*pedigree.matA.CHH[,1]
p0uu_in.CHH <- output.matA.DMR.centromeric.CHH$tmpp0
pedigree.matA.CHH <- data.frame(pedigree.matA.CHH, pairs.matA, dt)
colnames(pedigree.matA.CHH)[5:6] <- c("sample1", "sample2")

#Final pedigree CHH mat a
pedigree.mata.CHH <- output.mata.DMR.centromeric.CHH$Pdata
pedigree.mata.CHH[,1:3] <- pedigree.mata.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHH[,2] + pedigree.mata.CHH[,3] - 2*pedigree.mata.CHH[,1]
p0uu_in.CHH <- output.mata.DMR.centromeric.CHH$tmpp0
pedigree.mata.CHH <- data.frame(pedigree.mata.CHH, pairs.mata, dt)
colnames(pedigree.mata.CHH)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHH, pedigree.mata.CHH)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
plot.combined.pedigrees <- ggplot(pedigree.combined, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    geom_hline(yintercept = 0.03, lty = "dashed") +
    scale_y_continuous(limits = c(0,0.1)) +
    ggtitle("DMRs, CHH") +    
    xlab(myxlab) +
    ylab("Methylation divergence")

#Number of different samples
samples <- unique(c(as.character(pedigree.combined$sample1), as.character(pedigree.combined$sample2)))
n.samples <- length(samples) #69 samples in total

#Counting the number of times a sample is involved in an outlier observation
outlier.res <- data.frame(sample = samples, outlier.count = 0)
for(i in 1:n.samples) {
    current.sample <- samples[i]
    #Filter for those cases where the current sample is either one of the two samples being compared
    current <- filter(pedigree.combined, sample1 == current.sample | sample2 == current.sample)
    outlier.res[i,2] <- sum(current$D.value > 0.03) #Count number of outlier points
}

#Okay, this shows that there are nine samples that have that seems to cause the majority of the outlier points: MG5L20, MG5L20, MG20L17, ML1G40, ML6G40, ML10G40, MG5L30, MG20L29, MG20L30

outl.samp <- pedigree.combined$sample1 %in% c("MG5L20", "MG20L1", "MG20L17", "ML1G40", "ML6G40", "ML10G40", "MG5L30", "MG20L29", "MG20L30") | pedigree.combined$sample2 %in% c("MG5L20", "MG20L1", "MG20L17", "ML1G40", "ML6G40", "ML10G40", "MG5L30", "MG20L29", "MG20L30") 
pedigree.outl.rm <- pedigree.combined[!outl.samp,]

##Highlight pairwise comparisons that involved an outlier sample
pedigree.combined$highl <- "other"
pedigree.combined$highl[(pedigree.combined$sample1 == "MG5L20" | pedigree.combined$sample2 == "MG5L20")] <- "L20G5"
pedigree.combined$highl[pedigree.combined$sample1 == "MG20L1" | pedigree.combined$sample2 == "MG20L1"] <- "L1G20"
pedigree.combined$highl[pedigree.combined$sample1 == "MG20L17" | pedigree.combined$sample2 == "MG20L17"] <- "L17G20"
pedigree.combined$highl[pedigree.combined$sample1 == "ML1G40" | pedigree.combined$sample2 == "ML1G40"] <- "L1G40"
pedigree.combined$highl[(pedigree.combined$sample1 == "ML6G40" | pedigree.combined$sample2 == "ML6G40")] <- "L6G40"
pedigree.combined$highl[pedigree.combined$sample1 == "ML10G40" | pedigree.combined$sample2 == "ML10G40"] <- "L10G40"
pedigree.combined$highl[pedigree.combined$sample1 == "MG5L30" | pedigree.combined$sample2 == "MG5L30"] <- "L30G5"
pedigree.combined$highl[pedigree.combined$sample1 == "MG20L29" | pedigree.combined$sample2 == "MG20L29"] <- "L29G20"
pedigree.combined$highl[pedigree.combined$sample1 == "MG20L30" | pedigree.combined$sample2 == "MG20L30"] <- "L30G20"


plot.outl.high <- ggplot(pedigree.combined, aes(x = dt, y = D.value, colour = highl)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits, aes(x = dt, y = divsim), col = "red" ) +
    #geom_smooth() +
    #geom_line(aes(x = dt, y = pred), col = "blue") +
    #geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +
    geom_hline(yintercept = 0.03, lty = "dashed") +
    #scale_y_continuous(limits = c(0.0008, 0.0026)) +
    scale_colour_manual(values = c(other = "black", L20G5 = "red", L1G20 = "blue", L17G20 = "purple", L1G40 = "orange", L6G40 = "darkred", L10G40 = "darkblue", L30G5 = "darkcyan", L29G20 = "hotpink", L30G20 = "darkgoldenrod")) +
    scale_y_continuous(limits = c(0,0.1)) +    
    ggtitle("Outlier samples highlighted") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    theme(legend.title = element_blank())

##Divergence without outlier samples
plot.outl.rm <- ggplot(pedigree.outl.rm, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    #geom_line(data = theoryFits, aes(x = dt, y = divsim), col = "red" ) +
    #geom_smooth() +
    #geom_line(aes(x = dt, y = pred), col = "blue") +
    #geom_line(aes(x = dt, y = prednull), col = "blue", lty = "dotted") +
    scale_y_continuous(limits = c(0,0.1)) +    
    geom_hline(yintercept = 0.03, lty = "dashed") +
    #scale_y_continuous(limits = c(0.0008, 0.0026)) +    
    ggtitle("Outlier samples removed") +    
    xlab(myxlab) +
    ylab("Methylation divergence")

#Making the final outlier plot
outlier.plot <- plot_grid(plot.combined.pedigrees, plot.outl.high, plot.outl.rm, ncol = 3, rel_widths = c(1,1.3,1))
save_plot("./epigenet_ms/fig/outliers_DMR_centromeric.pdf",outlier.plot, base_width = 14, base_height = 4.25)


#How do the other contexts look like, once we have removed these outlier samples?
#Other context are not affected! Different outlier samples for the different contexts...
### *** Divergence across the whole genome without outliers

### Load the data for single cytosines and DMRs

#Single cytosines
load("./pedigree_data/all/MA_pedigree_CG_outl_rm.RData")
load("./pedigree_data/all/MA_pedigree_CHG_outl_rm.RData")
load("./pedigree_data/all/MA_pedigree_CHH_outl_rm.RData")

plotdata.ss <- data.frame(rbind(pedigree.all.ss.CG.outl.rm, pedigree.all.ss.CHG.outl.rm, pedigree.all.ss.CHH.outl.rm), context = c(rep("CG", nrow(pedigree.all.ss.CG.outl.rm)), rep("CHG", nrow(pedigree.all.ss.CHG.outl.rm)), rep("CHH", nrow(pedigree.all.ss.CHH.outl.rm))))

#Add AB fits
load("./modeldata/haploid/all/ABresults.haploid.ss.outl.rm.RData")
#Variables are
#neutral.haploid.ss.outl.rm.CG
#neutral.haploid.ss.outl.rm.CHG
#neutral.haploid.ss.outl.rm.CHH

#Theoretical fits
  theory.fit.data <- neutral.haploid.ss.outl.rm.CG$for.fit.plot
  theory.fits <- c(neutral.haploid.ss.outl.rm.CG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.ss.outlrm.CG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.ss.outl.rm.CHG$for.fit.plot
  theory.fits <- c(neutral.haploid.ss.outl.rm.CHG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.ss.outlrm.CHG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.ss.outl.rm.CHH$for.fit.plot
  theory.fits <- c(neutral.haploid.ss.outl.rm.CHH$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.ss.outlrm.CHH <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHH", length(theory.fits)))

 theoryFits.ss.outlrm <- rbind(theory.fit.ss.outlrm.CG, theory.fit.ss.outlrm.CHG, theory.fit.ss.outlrm.CHH)
                               
myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
plot.ss.div.all.outrm <- ggplot(plotdata.ss, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("Single sites") +
    geom_line(data = theoryFits.ss.outlrm, aes(x = dt, y = divsim), col = "orange" ) +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

### DMRs ###
load("./pedigree_data/DMR_all/MA_DMR_all_CG_outl_rm.RData")
load("./pedigree_data/DMR_all/MA_DMR_all_CHG_outl_rm.RData")
load("./pedigree_data/DMR_all/MA_DMR_all_CHH_outl_rm.RData")

#Add AB fits
load("./modeldata/haploid/all/ABresults.haploid.DMR.all.outl.rm.RData")
#Variables are
#neutral.haploid.DMR.all.outl.rm.CG
#neutral.haploid.DMR.all.outl.rm.CHG
#neutral.haploid.DMR.all.outl.rm.CHH

plotdata.DMR <- data.frame(rbind(pedigree.all.DMR.CG.outl.rm, pedigree.all.DMR.CHG.outl.rm, pedigree.all.DMR.CHH.outl.rm), context = c(rep("CG", nrow(pedigree.all.DMR.CG.outl.rm)), rep("CHG", nrow(pedigree.all.DMR.CHG.outl.rm)), rep("CHH", nrow(pedigree.all.DMR.CHH.outl.rm))))

#Theoretical fits
  theory.fit.data <- neutral.haploid.DMR.all.outl.rm.CG$for.fit.plot
  theory.fits <- c(neutral.haploid.DMR.all.outl.rm.CG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.DMR.outlrm.CG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.DMR.all.outl.rm.CHG$for.fit.plot
  theory.fits <- c(neutral.haploid.DMR.all.outl.rm.CHG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.DMR.outlrm.CHG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.DMR.all.outl.rm.CHH$for.fit.plot
  theory.fits <- c(neutral.haploid.DMR.all.outl.rm.CHH$estimates[1, "intercept"], theory.fit.data[,"div.sim"])
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
  theory.fit.DMR.outlrm.CHH <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHH", length(theory.fits)))

 theoryFits.DMR.outlrm <- rbind(theory.fit.DMR.outlrm.CG, theory.fit.DMR.outlrm.CHG, theory.fit.DMR.outlrm.CHH)

myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
plot.DMR.div.all.outrm <- ggplot(plotdata.DMR, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    geom_line(data = theoryFits.DMR.outlrm, aes(x = dt, y = divsim), col = "orange" ) +    
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

#Making the final plot
final.outlrm <- plot_grid(plot.ss.div.all.outrm, plot.DMR.div.all.outrm, labels = c("A", "B"), nrow = 2)
save_plot(filename = "methdiv_outlrm.pdf", final.outlrm, base_height = 8, base_width = 4 * 1.5 * 1.618)

#Saving the theory fits to include in other figure later
save(theoryFits.ss.outlrm, theoryFits.DMR.outlrm, file = "./data/theoryFits.outlrm.RData")


### Model comparisons against selection models for the outliers removed data set

### Single cytosines

#Run the null model

load("./pedigree_data/all/MA_pedigree_CG_outl_rm.RData")
load("./pedigree_data/all/MA_pedigree_CHG_outl_rm.RData")
load("./pedigree_data/all/MA_pedigree_CHH_outl_rm.RData")

outputABnull.CG.ss.outlrm <- ABnull(pedigree.data = pedigree.all.ss.CG.outl.rm, out.dir = getwd(), out.name = "./modeldata/ABnull_CG_all_ss_outlrm")

outputABnull.CHG.ss.outlrm <- ABnull(pedigree.data = pedigree.all.ss.CHG.outl.rm, out.dir = getwd(), out.name = "./modeldata/ABnull_CHG_all_ss_outlrm")

outputABnull.CHH.ss.outlrm <- ABnull(pedigree.data = pedigree.all.ss.CHH.outl.rm, out.dir = getwd(), out.name = "./modeldata/ABnull_CHH_all_ss_outlrm")


### 1. Do model comparisons

#Model comparisons for ss outliers removed data

## CG context
comp.out.CG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CG_all_ss_outl_rm.Rdata", pedigree.null = "./modeldata/ABnull_CG_all_ss_outlrm.Rdata") #Neutral model is preferred over the null model
comp.out.CG$Ftest # p < 2.2 * 10{ -16}

## CHG context
comp.out.CHG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CHG_all_ss_outl_rm.Rdata", pedigree.null = "./modeldata/ABnull_CHG_all_ss_outlrm.Rdata") #Neutral model is preferred over the null model
comp.out.CHG$Ftest
#p  < 2.2 * 10{ -16}

### CHH context
comp.out.CHH <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CHH_all_ss_outl_rm.Rdata", pedigree.null = "./modeldata/ABnull_CHH_all_ss_outlrm.Rdata") #Neutral model is preferred over the null model
comp.out.CHG$Ftest
#p  < 2.2 * 10{ -16}



### DMRs

#Run the null model

load("./pedigree_data/DMR_all/MA_DMR_all_CG_outl_rm.RData")
load("./pedigree_data/DMR_all/MA_DMR_all_CHG_outl_rm.RData")
load("./pedigree_data/DMR_all/MA_DMR_all_CHH_outl_rm.RData")

outputABnull.CG.DMR.outlrm <- ABnull(pedigree.data = pedigree.all.DMR.CG.outl.rm, out.dir = getwd(), out.name = "./modeldata/ABnull_CG_all_DMR_outlrm")

outputABnull.CHG.DMR.outlrm <- ABnull(pedigree.data = pedigree.all.DMR.CHG.outl.rm, out.dir = getwd(), out.name = "./modeldata/ABnull_CHG_all_DMR_outlrm")

outputABnull.CHH.DMR.outlrm <- ABnull(pedigree.data = pedigree.all.DMR.CHH.outl.rm, out.dir = getwd(), out.name = "./modeldata/ABnull_CHH_all_DMR_outlrm")

#Model comparisons for DMRs outliers removed data
## CG context
comp.out.CG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CG_all_DMR_outl_rm.Rdata", pedigree.null = "./modeldata/ABnull_CG_all_DMR_outlrm.Rdata") #Neutral model is preferred over the null model
comp.out.CG$Ftest
# p < 2.2 \times 10^{-16}

## CHG context
comp.out.CHG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CHG_all_DMR_outl_rm.Rdata", pedigree.null = "./modeldata/ABnull_CHG_all_DMR_outlrm.Rdata") #Neutral model is preferred over the null model
comp.out.CHG$Ftest
#p < 2.2 \times 10^{-16}

### CHH context
comp.out.CHH <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/all/ABneutralHaploid_CHH_all_DMR_outl_rm.Rdata", pedigree.null = "./modeldata/ABnull_CHH_all_DMR_outlrm.Rdata") #Neutral model is preferred over the null model
comp.out.CHH$Ftest
#p < 2.2 \times 10^{-16}


### Checking only initial divergence, deltat < 600 mitoses

ggplot(filter(plotdata.ss, dt < 600), aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("Single sites") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)

ggplot(filter(plotdata.DMR, dt < 600), aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("DMRs") +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_wrap( ~ context)
### * Initial linear divergence for WGBS centromeric H3K9me3

### Checking the initial linear divergence for WGBS data for centromeric H3K9me3
### Do we get consistent results with Nanopore data?

#Samples involved in pairwise comparisons

#mat A samples
genTable <- fread("nodelist_matA.fn")
genTable <- genTable %>% filter(genTable$meth == "Y")
pairs <- combn(genTable$filename, 2)
pairs <- t(pairs) #transpose
#fix the sample names
pairs[,1] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,1])
pairs[,2] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,2])
pairs[,1] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,1])
pairs[,2] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,2])
pairs.matA <- pairs

#mat a samples
genTable <- fread("nodelist_mata.fn")
genTable <- genTable %>% filter(genTable$meth == "Y")
pairs <- combn(genTable$filename, 2)
pairs <- t(pairs) #transpose
pairs[,1] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,1])
pairs[,2] <- gsub("~/Genomics/Neurospora/methylation/methimpute/methylome_", "", pairs[,2])
pairs[,1] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,1])
pairs[,2] <- gsub("_bismark_pe_trim_All.txt", "", pairs[,2])
pairs.mata <- pairs

### ** Single sites

#Data for centromericH3K9 single sites
load(file = "./pedigree_data/centromericH3K9/MA_mat_CG_centromericH3K9.RData") #output.matA.centromericH3K9.CG, output.mata.centromericH3K9.CG
load(file = "./pedigree_data/centromericH3K9/MA_mat_CHG_centromericH3K9.RData") #output.matA.centromericH3K9.CHG, output.mata.centromericH3K9.CHG
load(file = "./pedigree_data/centromericH3K9/MA_mat_CHH_centromericH3K9.RData") #output.matA.CHH

### CG ###
#Final pedigree CG matA
pedigree.matA.CG <- output.matA.centromericH3K9.CG$Pdata
pedigree.matA.CG[,1:3] <- pedigree.matA.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CG[,2] + pedigree.matA.CG[,3] - 2*pedigree.matA.CG[,1]
#p0uu_in.CG <- output.matA.centromericH3K9.CG$tmpp0
pedigree.matA.CG <- data.frame(pedigree.matA.CG, pairs.matA, dt)
colnames(pedigree.matA.CG)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CG <- output.mata.centromericH3K9.CG$Pdata
pedigree.mata.CG[,1:3] <- pedigree.mata.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CG[,2] + pedigree.mata.CG[,3] - 2*pedigree.mata.CG[,1]
#p0uu_in.CG <- output.mata.centromericH3K9.CG$tmpp0
pedigree.mata.CG <- data.frame(pedigree.mata.CG, pairs.mata, dt)
colnames(pedigree.mata.CG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CG, pedigree.mata.CG)

### Here can remove some outliers if that is necessary

### Filter for the early linear stages
pedigree.centromeric.ss.CG <- filter(pedigree.combined, dt < 400)

### CHG ###
#Final pedigree CG matA
pedigree.matA.CHG <- output.matA.centromericH3K9.CHG$Pdata
pedigree.matA.CHG[,1:3] <- pedigree.matA.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHG[,2] + pedigree.matA.CHG[,3] - 2*pedigree.matA.CHG[,1]
#p0uu_in.CHG <- output.matA.centromericH3K9.CHG$tmpp0
pedigree.matA.CHG <- data.frame(pedigree.matA.CHG, pairs.matA, dt)
colnames(pedigree.matA.CHG)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CHG <- output.mata.centromericH3K9.CHG$Pdata
pedigree.mata.CHG[,1:3] <- pedigree.mata.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHG[,2] + pedigree.mata.CHG[,3] - 2*pedigree.mata.CHG[,1]
#p0uu_in.CHG <- output.mata.centromericH3K9.CHG$tmpp0
pedigree.mata.CHG <- data.frame(pedigree.mata.CHG, pairs.mata, dt)
colnames(pedigree.mata.CHG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHG, pedigree.mata.CHG)

### Here can remove some outliers if that is necessary

### Filter for the early linear stages
pedigree.centromeric.ss.CHG <- filter(pedigree.combined, dt < 400)

### CHH ###
#Final pedigree CG matA
pedigree.matA.CHH <- output.matA.centromericH3K9.CHH$Pdata
pedigree.matA.CHH[,1:3] <- pedigree.matA.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHH[,2] + pedigree.matA.CHH[,3] - 2*pedigree.matA.CHH[,1]
p0uu_in.CHH <- output.matA.centromericH3K9.CHH$tmpp0
pedigree.matA.CHH <- data.frame(pedigree.matA.CHH, pairs.matA, dt)
colnames(pedigree.matA.CHH)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CHH <- output.mata.centromericH3K9.CHH$Pdata
pedigree.mata.CHH[,1:3] <- pedigree.mata.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHH[,2] + pedigree.mata.CHH[,3] - 2*pedigree.mata.CHH[,1]
p0uu_in.CHH <- output.mata.centromericH3K9.CHH$tmpp0
pedigree.mata.CHH <- data.frame(pedigree.mata.CHH, pairs.mata, dt)
colnames(pedigree.mata.CHH)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHH, pedigree.mata.CHH)

### Here can remove some outliers if that is necessary

### Filter for the early linear stages
pedigree.centromeric.ss.CHH <- filter(pedigree.combined, dt < 400)

### Load the propertions of unmethylated cytosines in the MA ancestor
load("./data/cytosine.prop.MA.anc.centromeric.RData")
p0uu_in.CG <- cyt.props.wgbs.centromeric[2,3] #Unmethylated CG
p0uu_in.CHG <- cyt.props.wgbs.centromeric[4,3] #Unmethylated CHG
p0uu_in.CHH <- cyt.props.wgbs.centromeric[6,3] #Unmethylated CHH


#Plotting
#Make plot dataframe
plotdata <- data.frame(rbind(pedigree.centromeric.ss.CG, pedigree.centromeric.ss.CHG, pedigree.centromeric.ss.CHH), context = c(rep("CG", nrow(pedigree.centromeric.ss.CG)), rep("CHG", nrow(pedigree.centromeric.ss.CHG)), rep("CHH", nrow(pedigree.centromeric.ss.CHH))))


### Simple tests for neutral and non-neutral models ###
#Run the models
### CG ###
neutral.haploid.ss.init.centromeric.CG <- ABneutralHaploid(pedigree.data = pedigree.centromeric.ss.CG[,1:4], p0uu = p0uu_in.CG,
eqp = p0uu_in.CG, eqp.weight = 0.001, Nstarts = 100, out.dir = "~/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_ss_CG")

neutral.haploid.ss.init.centromeric.CG$estimates[1,]
#alpha = 1.290712e-05
#beta = 6.156157e-05
#alpha is consistent with nanopore, beta slightly lower

#Null model
outputABnull <- ABnull(pedigree.data = pedigree.centromeric.ss.CG[,1:4], out.dir = "/home/ililkron/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "./modeldata/ABnull_initial_centromericH3K9_ss_CG_test")

##########

### CHG ###
neutral.haploid.ss.init.centromeric.CHG <- ABneutralHaploid(pedigree.data = pedigree.centromeric.ss.CHG[,1:4], p0uu = p0uu_in.CHG,
eqp = p0uu_in.CHG, eqp.weight = 0.001, Nstarts = 100, out.dir = "~/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_ss_CHG")

neutral.haploid.ss.init.centromeric.CHG$estimates[1,]
#alpha = 1.401567e-05
#beta = 7.786028e-05
#alpha and beta are consistent with nanopore

#Null model
outputABnull <- ABnull(pedigree.data = pedigree.centromeric.ss.CHG[,1:4], out.dir = "/home/ililkron/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "./modeldata/ABnull_initial_centromericH3K9_ss_CHG_test")

##########

### CHH ###
neutral.haploid.ss.init.centromeric.CHH <- ABneutralHaploid(pedigree.data = pedigree.centromeric.ss.CHH[,1:4], p0uu = p0uu_in.CHH,
eqp = p0uu_in.CHH, eqp.weight = 0.001, Nstarts = 100, out.dir = "~/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_ss_CHH")

neutral.haploid.ss.init.centromeric.CHH$estimates[1,]
#alpha = 1.027627e-05
#beta = 5.846393e-05
#alpha is slightly less than nanopore, beta is consistent with nanopore 

#Null model
outputABnull <- ABnull(pedigree.data = pedigree.centromeric.ss.CHH[,1:4], out.dir = "/home/ililkron/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "./modeldata/ABnull_initial_centromericH3K9_ss_CHH_test")


#Model comparisons 
comp.out.CG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_ss_CG.Rdata", pedigree.null = "./modeldata/ABnull_initial_centromericH3K9_ss_CG_test.Rdata")
#Neutral model preferred over null model
comp.out.CG$Ftest
#p-value = 3.285389e-06

comp.out.CHG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_ss_CHG.Rdata", pedigree.null = "./modeldata/ABnull_initial_centromericH3K9_ss_CHG_test.Rdata")
#Neutral model preferred over null model
comp.out.CHG$Ftest
#p-value = 8.963347e-11

comp.out.CHH <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_ss_CHH.Rdata", pedigree.null = "./modeldata/ABnull_initial_centromericH3K9_ss_CHH_test.Rdata")
#Neutral model preferred over null model
comp.out.CHH$Ftest
#p-value = 1.669695e-06

#################################################################

### Making a dataset for running the AB models on CSC cluster ###
#save(pedigree.centromeric.ss.CG, p0uu_in.CG, pedigree.centromeric.ss.CHG, p0uu_in.CHG, pedigree.centromeric.ss.CHH, p0uu_in.CHH, file = "./pedigree_data/centromeric_initial/MA_pedigree_centromeric_ss_initial.RData")


#Theoretical fits
  theory.fit.data <- neutral.haploid.ss.init.centromeric.CG$for.fit.plot
  theory.fits <- c(neutral.haploid.ss.init.centromeric.CG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])[1:300]
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"][1:299])
  theory.fit.CG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.ss.init.centromeric.CHG$for.fit.plot
  theory.fits <- c(neutral.haploid.ss.init.centromeric.CHG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])[1:300]
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"][1:299])
  theory.fit.CHG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.ss.init.centromeric.CHH$for.fit.plot
  theory.fits <- c(neutral.haploid.ss.init.centromeric.CHH$estimates[1, "intercept"], theory.fit.data[,"div.sim"])[1:300]
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"][1:299])
  theory.fit.CHH <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHH", length(theory.fits)))

theoryFits.ss.init <- rbind(theory.fit.CG, theory.fit.CHG, theory.fit.CHH)

#Saving the results
save(neutral.haploid.ss.init.centromeric.CG, neutral.haploid.ss.init.centromeric.CHG, neutral.haploid.ss.init.centromeric.CHH, plotdata, theoryFits.ss.init, file = "./modeldata/haploid/ABresults.haploid.ss.WGBS.init.RData")

#Load the data
load("./modeldata/haploid/ABresults.haploid.ss.WGBS.init.RData")

#Test
myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
init.ss.plot <- ggplot(plotdata, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("Single sites, centromeric H3K9, WGBS, initial divergence") +
    geom_line(data = theoryFits.ss.init, aes(x = dt, y = divsim), col = "blue" ) +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    #scale_y_continuous(limits = c(0, 0.1)) +    
    facet_wrap( ~ context)

### ** DMRs

#Data for centromericH3K9 single sites
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_mat_CG.RData") #output.matA.DMR.centromeric.CG, output.mata.DMR.centromeric.CG
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_mat_CHG.RData") #output.matA.DMR.centromeric.CHG, output.mata.DMR.centromeric.CHG
load(file = "./pedigree_data/DMR_centromericH3K9/MA_DMR_centromericH3K9_mat_CHH.RData") #output.matA.DMR.centromeric.CHH

### CG ###
#Final pedigree CG matA
pedigree.matA.CG <- output.matA.DMR.centromeric.CG$Pdata
pedigree.matA.CG[,1:3] <- pedigree.matA.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CG[,2] + pedigree.matA.CG[,3] - 2*pedigree.matA.CG[,1]
#p0uu_in.CG <- output.matA.centromericH3K9.CG$tmpp0
pedigree.matA.CG <- data.frame(pedigree.matA.CG, pairs.matA, dt)
colnames(pedigree.matA.CG)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CG <- output.mata.DMR.centromeric.CG$Pdata
pedigree.mata.CG[,1:3] <- pedigree.mata.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CG[,2] + pedigree.mata.CG[,3] - 2*pedigree.mata.CG[,1]
#p0uu_in.CG <- output.mata.centromericH3K9.CG$tmpp0
pedigree.mata.CG <- data.frame(pedigree.mata.CG, pairs.mata, dt)
colnames(pedigree.mata.CG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CG, pedigree.mata.CG)

### Here can remove some outliers if that is necessary

### Filter for the early linear stages
pedigree.centromeric.DMR.CG <- filter(pedigree.combined, dt < 400)

### CHG ###
#Final pedigree CG matA
pedigree.matA.CHG <- output.matA.DMR.centromeric.CHG$Pdata
pedigree.matA.CHG[,1:3] <- pedigree.matA.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHG[,2] + pedigree.matA.CHG[,3] - 2*pedigree.matA.CHG[,1]
#p0uu_in.CHG <- output.matA.centromericH3K9.CHG$tmpp0
pedigree.matA.CHG <- data.frame(pedigree.matA.CHG, pairs.matA, dt)
colnames(pedigree.matA.CHG)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CHG <- output.mata.DMR.centromeric.CHG$Pdata
pedigree.mata.CHG[,1:3] <- pedigree.mata.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHG[,2] + pedigree.mata.CHG[,3] - 2*pedigree.mata.CHG[,1]
#p0uu_in.CHG <- output.mata.centromericH3K9.CHG$tmpp0
pedigree.mata.CHG <- data.frame(pedigree.mata.CHG, pairs.mata, dt)
colnames(pedigree.mata.CHG)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHG, pedigree.mata.CHG)

### Here can remove some outliers if that is necessary

### Filter for the early linear stages
pedigree.centromeric.DMR.CHG <- filter(pedigree.combined, dt < 400)

### CHH ###
#Final pedigree CG matA
pedigree.matA.CHH <- output.matA.DMR.centromeric.CHH$Pdata
pedigree.matA.CHH[,1:3] <- pedigree.matA.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.matA.CHH[,2] + pedigree.matA.CHH[,3] - 2*pedigree.matA.CHH[,1]
#p0uu_in.CHH <- output.matA.centromericH3K9.CHH$tmpp0
pedigree.matA.CHH <- data.frame(pedigree.matA.CHH, pairs.matA, dt)
colnames(pedigree.matA.CHH)[5:6] <- c("sample1", "sample2")

#Final pedigree CG mat a
pedigree.mata.CHH <- output.mata.DMR.centromeric.CHH$Pdata
pedigree.mata.CHH[,1:3] <- pedigree.mata.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.mata.CHH[,2] + pedigree.mata.CHH[,3] - 2*pedigree.mata.CHH[,1]
#p0uu_in.CHH <- output.mata.centromericH3K9.CHH$tmpp0
pedigree.mata.CHH <- data.frame(pedigree.mata.CHH, pairs.mata, dt)
colnames(pedigree.mata.CHH)[5:6] <- c("sample1", "sample2")

#Combine data from the two pedigrees
pedigree.combined <- rbind(pedigree.matA.CHH, pedigree.mata.CHH)

### Here can remove some outliers if that is necessary

### Filter for the early linear stages
pedigree.centromeric.DMR.CHH <- filter(pedigree.combined, dt < 400)

### Load the proportions of unmethylated cytosines in the MA ancestor
load("./data/DMR.prop.MA.anc.cent.RData")
p0uu_in.CG <- DMR.props.wgbs.cent[2,3] #Unmethylated CG
p0uu_in.CHG <- DMR.props.wgbs.cent[4,3] #Unmethylated CHG
p0uu_in.CHH <- DMR.props.wgbs.cent[6,3] #Unmethylated CHH

#Plotting
#Make plot dataframe
plotdata.DMR <- data.frame(rbind(pedigree.centromeric.DMR.CG, pedigree.centromeric.DMR.CHG, pedigree.centromeric.DMR.CHH), context = c(rep("CG", nrow(pedigree.centromeric.DMR.CG)), rep("CHG", nrow(pedigree.centromeric.DMR.CHG)), rep("CHH", nrow(pedigree.centromeric.DMR.CHH))))


### Simple tests for neutral and non-neutral models ###
#Run the models
### CG ###
neutral.haploid.DMR.init.centromeric.CG <- ABneutralHaploid(pedigree.data = pedigree.centromeric.DMR.CG[,1:4], p0uu = p0uu_in.CG,
eqp = p0uu_in.CG, eqp.weight = 0.001, Nstarts = 100, out.dir = "~/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_DMR_CG")

neutral.haploid.DMR.init.centromeric.CG$estimates[1,]
#alpha = 6.169808e-05
#beta = 0.0002560612 = 25.60e-05
#alpha and beta are higher than the nanopore estimate (beta especially so)

#Null model
outputABnull <- ABnull(pedigree.data = pedigree.centromeric.DMR.CG[,1:4], out.dir = "/home/ililkron/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "./modeldata/ABnull_initial_centromericH3K9_DMR_CG_test")

### CHG ###
#Test with nanopore DMR methylation proportions
p0uu_in.CHG <- DMR.props.nano.cent[4,3]

#Run with wgbs DMR methylation proportions
p0uu_in.CHG <- DMR.props.wgbs.cent[4,3]
neutral.haploid.DMR.init.centromeric.CHG <- ABneutralHaploid(pedigree.data = pedigree.centromeric.DMR.CHG[,1:4], p0uu = p0uu_in.CHG,
eqp = p0uu_in.CHG, eqp.weight = 0.001, Nstarts = 100, out.dir = "~/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_DMR_CHG")

neutral.haploid.DMR.init.centromeric.CHG$estimates[1,]
#alpha = 3.982275e-05
#beta = 1.70e-03
#alpha is a bit higher and beta is significantly higher than nanopore estimate
#Model has numerical convergence issues


#By using nanopore unmethylated proportions we get
#alpha = 3.386421e-05
#beta = 24.17e-05
#Alpha and beta are quite a bit higher than for nanopore

#Null model
outputABnull <- ABnull(pedigree.data = pedigree.centromeric.DMR.CHG[,1:4], out.dir = "/home/ililkron/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "./modeldata/ABnull_initial_centromericH3K9_DMR_CHG_test")

### CHH ###
neutral.haploid.DMR.init.centromeric.CHH <- ABneutralHaploid(pedigree.data = pedigree.centromeric.DMR.CHH[,1:4], p0uu = p0uu_in.CHH,
eqp = p0uu_in.CHH, eqp.weight = 0.001, Nstarts = 100, out.dir = "~/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_DMR_CHH")

neutral.haploid.DMR.init.centromeric.CHH$estimates[1,]
#alpha = 4.605432e-05
#beta = 1.16e-03
#alpha is comparable with nanopore but beta is significantly higher. 
#Model has numerical convergence issues

#Null model
outputABnull <- ABnull(pedigree.data = pedigree.centromeric.DMR.CHH[,1:4], out.dir = "/home/ililkron/Documents/tutkijatohtori/epimutation/MA_meth", out.name = "./modeldata/ABnull_initial_centromericH3K9_DMR_CHH_test")

### Model comparisons ###
comp.out.CG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_DMR_CG.Rdata", pedigree.null = "./modeldata/ABnull_initial_centromericH3K9_DMR_CG_test.Rdata")
#Neutral model preferred over null model
comp.out.CG$Ftest
#p-value = 1.074212e-11

comp.out.CHG <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_DMR_CHG.Rdata", pedigree.null = "./modeldata/ABnull_initial_centromericH3K9_DMR_CHG_test.Rdata")
#Neutral model preferred over null model
comp.out.CHG$Ftest
#p-value = 1.174182e-03

comp.out.CHH <- FtestRSSHaploid(pedigree.select = "./modeldata/haploid/centromericH3K9/ABneutralHaploid_initial_centromeric_DMR_CHH.Rdata", pedigree.null = "./modeldata/ABnull_initial_centromericH3K9_DMR_CHH_test.Rdata")
#Neutral model preferred over null model
comp.out.CHH$Ftest
#p-value = 0.1897764

#Theoretical fits
  theory.fit.data <- neutral.haploid.DMR.init.centromeric.CG$for.fit.plot
  theory.fits <- c(neutral.haploid.DMR.init.centromeric.CG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])[1:300]
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"][1:299])
  theory.fit.CG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.DMR.init.centromeric.CHG$for.fit.plot
  theory.fits <- c(neutral.haploid.DMR.init.centromeric.CHG$estimates[1, "intercept"], theory.fit.data[,"div.sim"])[1:300]
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"][1:299])
  theory.fit.CHG <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHG", length(theory.fits)))

  theory.fit.data <- neutral.haploid.DMR.init.centromeric.CHH$for.fit.plot
  theory.fits <- c(neutral.haploid.DMR.init.centromeric.CHH$estimates[1, "intercept"], theory.fit.data[,"div.sim"])[1:300]
  theory.fit.t <- c(0, theory.fit.data[,"delta.t"][1:299])
  theory.fit.CHH <- data.frame("divsim" = theory.fits, "dt" = theory.fit.t, context = rep("CHH", length(theory.fits)))

theoryFits.DMR.init <- rbind(theory.fit.CG, theory.fit.CHG, theory.fit.CHH)

#Saving the results
save(neutral.haploid.DMR.init.centromeric.CG, neutral.haploid.DMR.init.centromeric.CHG, neutral.haploid.DMR.init.centromeric.CHH, plotdata.DMR, theoryFits.DMR.init, file = "./modeldata/haploid/ABresults.haploid.DMR.WGBS.init.RData")

#Test
myxlab <- TeX("$\\Delta$t (mitoses)")
#Making a plot
init.DMR.plot <- ggplot(plotdata.DMR, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    ggtitle("DMRs, centromeric H3K9, WGBS, initial divergence") +
    geom_line(data = theoryFits.DMR.init, aes(x = dt, y = divsim), col = "blue" ) +    
    xlab(myxlab) +
    ylab("Methylation divergence") +
    #scale_y_continuous(limits = c(0, 0.12)) +    
    facet_wrap( ~ context)


#making the final plot
init.final.plot <- plot_grid(init.ss.plot, init.DMR.plot, nrow = 2, labels = c("A", "B"))

#Saving the plot
save_plot("./epigenet_ms/fig/initial_wgbs_H3K9cent_div.pdf", init.final.plot, base_height = 8, base_width = 4*1.5*1.618)

#save_plot("./epigenet_ms/fig/H3K9cent_div.pdf", plot.H3K9cent.combined, base_height = 8, base_width = 4 * 1.5 * 1.618)
### * Summarize sequencing report statistics

datafolder <- "./reports/"

samples <- list.files(datafolder)

nsamples <- length(samples) #Number of samples

#Make the results matrix
res.mat <- data.frame(matrix(rep(0, nsamples*4), ncol = 4))
colnames(res.mat) <- c("sample", "reads", "mean.coverage", "mean.depth")

##Loop over all samples in the reports folder
for(i in 1:nsamples) {
    sample <- samples[i]
    work.dir <- datafolder

    ### Coverage information from samtools
    coverage.file <- paste(work.dir, sample, sep = "") #Set the coverage file
    coverage.table <- read.table(coverage.file, header = F, sep = "\t")
    colnames(coverage.table) <- c("chr", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")

    #Store the sample name
    sample <- gsub("trim.", "", sample)
    sample <- gsub("_1_bismark_bt2_pe_coverage.txt", "", sample)
    res.mat[i,1] <- sample

    #Storing the number of reads
    res.mat[i,2] <- sum(coverage.table[,4]) #All reads across chromosomes

    #Calculating coverage across the whole genome (weighted average)
    #weights average across all contigs
    weights <- coverage.table$endpos / sum(coverage.table$endpos)
    mean.coverage <- sum(weights*coverage.table$coverage) / sum(weights)
    #Store coverage data
    res.mat[i,3] <- round(mean.coverage,2)
    
    #Calculating mean depth across the whole genome
    #weighted average across all contigs
    weights <- coverage.table$covbases / sum(coverage.table$covbases)
    mean.depth.w <- sum(weights*coverage.table$meandepth) / sum(weights)
    #Store the mean depth data
    res.mat[i,4] <- round(mean.depth.w,2)
}

#save(res.mat, file = paste(datafolder, "seq_report.RData", sep = ""))
write.table(res.mat, file = paste(datafolder, "seq_report.csv", sep = ""), quote = F, sep = ",", row.names = F)

#Check bisulfite conversion efficiency data
load("./reports/bs_conv_eff.RData")

#Fix sample names
conv.eff[,1] <- gsub("methylome_", "", conv.eff[,1])
conv.eff[,1] <- gsub("_bismark_pe_trim_All.txt", "", conv.eff[,1])

write.table(conv.eff, file = paste(datafolder, "bs_conversion.csv", sep = ""), quote = F, sep = ",", row.names = F)
### * Methylation proportion in the MA lines over the course of the experiment

### ** Observed methylation frequencies in the ancestor and MA lines

## The idea is to look at observed methylation frequencies in the MA lines
## To address the criticism that methylation divergence could be just the lines gradually losing methylation


# Methylation proportions for centromeric regions were calculated in the cluster

#Load the data
load("~/Documents/tutkijatohtori/epimutation/MA_meth/data/WGBS.cyt.prop.sample.RData")

#Processing the sample string to extract line and transfer
#Doing some regular expression magick
WGBS.cent.cyt.prop.by.sample$Line <- str_extract(WGBS.cent.cyt.prop.by.sample$sample, "(?<=L)\\d+")
WGBS.cent.cyt.prop.by.sample$Line[1:6] <- "Anc"

WGBS.cent.cyt.prop.by.sample$Transfer <- str_extract(WGBS.cent.cyt.prop.by.sample$sample, "(?<=G)\\d+")
WGBS.cent.cyt.prop.by.sample$Transfer[1:6] <- "0"
WGBS.cent.cyt.prop.by.sample$Transfer <- as.integer(WGBS.cent.cyt.prop.by.sample$Transfer)

###Store mating type of each line
WGBS.cent.cyt.prop.by.sample$mat <- ifelse(as.numeric(WGBS.cent.cyt.prop.by.sample$Line) <= 20, "mat A", "mat a")
WGBS.cent.cyt.prop.by.sample$mat[1:6] <- c(rep("mat A",3), rep("mat a", 3))

colnames(WGBS.cent.cyt.prop.by.sample)[2:4] <- c("CG", "CHG", "CHH")
WGBS.cytprop.long <- pivot_longer(WGBS.cent.cyt.prop.by.sample, c(2,3,4), names_to = "context", values_to = "prop")

fit <-lm(lin$y~0 +lin$x,offset=rep(10,length(lin$x)))

#Next we need to fit regression for each MA line
lines <- WGBS.cent.cyt.prop.by.sample$Line[-c(1:6)]
mats <- WGBS.cent.cyt.prop.by.sample$mat[-c(1:6)]
preddata <- data.frame(Transfer = c(0,5,20,40))
WGBS.data <- NULL
for(i in 1:length(lines)) {
    #CG context
    cur.data <- filter(WGBS.cytprop.long, (Line == lines[i] | Line == "Anc") & mat == mats[i] & context == "CG")
    anc.prop <- mean(filter(cur.data, Line == "Anc")$prop)
    cur.lm <- lm(I(prop - anc.prop) ~ 0 + Transfer, data = cur.data)
    prop <- predict(cur.lm, newdata = preddata) + anc.prop
    cur.res <- data.frame(Line = lines[i], Transfer = c(0,5,20,40), prop = prop, context = "CG")
    WGBS.data <- rbind(WGBS.data, cur.res) #Storing the results

    #CHG context
    cur.data <- filter(WGBS.cytprop.long, (Line == lines[i] | Line == "Anc") & mat == mats[i] & context == "CHG")
    anc.prop <- mean(filter(cur.data, Line == "Anc")$prop)
    cur.lm <- lm(I(prop - anc.prop) ~ 0 + Transfer, data = cur.data)
    prop <- predict(cur.lm, newdata = preddata) + anc.prop
    cur.res <- data.frame(Line = lines[i], Transfer = c(0,5,20,40), prop = prop, context = "CHG")
    WGBS.data <- rbind(WGBS.data, cur.res) #Storing the results

    #CHH context
    cur.data <- filter(WGBS.cytprop.long, (Line == lines[i] | Line == "Anc") & mat == mats[i] & context == "CHG")
    anc.prop <- mean(filter(cur.data, Line == "Anc")$prop)
    cur.lm <- lm(I(prop - anc.prop) ~ 0 + Transfer, data = cur.data)
    prop <- predict(cur.lm, newdata = preddata) + anc.prop
    cur.res <- data.frame(Line = lines[i], Transfer = c(0,5,20,40), prop = prop, context = "CHH")
    WGBS.data <- rbind(WGBS.data, cur.res) #Storing the results
} #Done looping over all lines


#Load the Nanopore data
load("~/Documents/tutkijatohtori/epimutation/MA_meth/data/Nano.cyt.prop.sample.RData")

#Processing the sample string to extract line and transfer
#Doing some regular expression magick
Nano.cent.cyt.prop.by.sample$Line <- str_extract(Nano.cent.cyt.prop.by.sample$sample, "(?<=L)\\d+")
Nano.cent.cyt.prop.by.sample$Line[1:2] <- "Anc"

Nano.cent.cyt.prop.by.sample$Transfer <- str_extract(Nano.cent.cyt.prop.by.sample$sample, "(?<=G)\\d+")
Nano.cent.cyt.prop.by.sample$Transfer[1:2] <- "0"
Nano.cent.cyt.prop.by.sample$Transfer <- as.integer(Nano.cent.cyt.prop.by.sample$Transfer)

###Store mating type of each line
Nano.cent.cyt.prop.by.sample$mat <- ifelse(as.numeric(Nano.cent.cyt.prop.by.sample$Line) <= 20, "mat A", "mat a")
Nano.cent.cyt.prop.by.sample$mat[1:2] <- c("mat a", "mat A")

colnames(Nano.cent.cyt.prop.by.sample)[2:4] <- c("CG", "CHG", "CHH")
Nano.cytprop.long <- pivot_longer(Nano.cent.cyt.prop.by.sample, c(2,3,4), names_to = "context", values_to = "prop")



#Next we need to fit regression for each MA line
lines <- Nano.cent.cyt.prop.by.sample$Line[-c(1:2)]
mats <- Nano.cent.cyt.prop.by.sample$mat[-c(1:2)]
preddata <- data.frame(Transfer = c(0,1,5,7,8,10,15))
Nano.data <- NULL
for(i in 1:length(lines)) {
    #CG context
    cur.data <- filter(Nano.cytprop.long, (Line == lines[i] | Line == "Anc") & mat == mats[i] & context == "CG")
    anc.prop <- mean(filter(cur.data, Line == "Anc")$prop)
    cur.lm <- lm(I(prop - anc.prop) ~ 0 + Transfer, data = cur.data)
    prop <- predict(cur.lm, newdata = preddata) + anc.prop
    cur.res <- data.frame(Line = lines[i], Transfer =c(0,1,5,7,8,10,15), prop = prop, context = "CG")
    Nano.data <- rbind(Nano.data, cur.res) #Storing the results

    #CHG context
    cur.data <- filter(Nano.cytprop.long, (Line == lines[i] | Line == "Anc") & mat == mats[i] & context == "CHG")
    anc.prop <- mean(filter(cur.data, Line == "Anc")$prop)
    cur.lm <- lm(I(prop - anc.prop) ~ 0 + Transfer, data = cur.data)
    prop <- predict(cur.lm, newdata = preddata) + anc.prop
    cur.res <- data.frame(Line = lines[i], Transfer = c(0,1,5,7,8,10,15), prop = prop, context = "CHG")
    Nano.data <- rbind(Nano.data, cur.res) #Storing the results

    #CHH context
    cur.data <- filter(Nano.cytprop.long, (Line == lines[i] | Line == "Anc") & mat == mats[i] & context == "CHH")
    anc.prop <- mean(filter(cur.data, Line == "Anc")$prop)
    cur.lm <- lm(I(prop - anc.prop) ~ 0 + Transfer, data = cur.data)
    prop <- predict(cur.lm, newdata = preddata) + anc.prop
    cur.res <- data.frame(Line = lines[i], Transfer = c(0,1,5,7,8,10,15), prop = prop, context = "CHH")
    Nano.data <- rbind(Nano.data, cur.res) #Storing the results
} #Done looping over all lines


## Making figures

#Proportion of methylated cytosines in centromeric regions, WGBS data
p.WGBS <- ggplot(WGBS.data, aes(y = prop, x = Transfer, group = Line)) +
    geom_line(colour = "blue", alpha = 0.2) +
    #geom_point() +
    ylab("Prop. of methylated cytosines") +    
    facet_wrap(.~context)

p.nano <- ggplot(Nano.data, aes(y = prop, x = Transfer, group = Line)) +
    geom_line(colour = "blue", alpha = 0.2) +
    #geom_point() +
    ylab("Prop. of methylated cytosines") +    
    facet_wrap(.~context)


## Get the number of gains and losses for centromeric regions for the Nanopore data

### Cytosines ###
datafolder <- "~/Documents/tutkijatohtori/epimutation/MA_meth/data/"
aineisto.cytosines <- read.csv(paste0(datafolder, "cytosines.intervals.centromeric.csv"), header = T)

### DMRs ###
aineisto.DMR <- read.csv(paste0(datafolder, "results_analysis_DMRS.csv"), header = T)
#nomut <- filter(aineisto, Num_mutation == 0)
aineisto.DMR$interval <- rep(c("Anc--T1", "T1--T5", "T5--T7", "T7--T8", "T8--T10", "T10--T15"), 6)
#Need to have matching names in the genetic mutation data and DMRs
aineisto.DMR$Line <- fct_recode(aineisto.DMR$Line, "L02" = "L2")

#Load DMRs from the nanopore dataset
load(file = "./data/DMRs_all_nanopore.RData") #In DMRs.all.nanopore

#Filtering only for centromeric DMRs
DMRs.centromeric.nanopore <- filter(DMRs.all.nanopore, domain == "Centromeric")
cent.DMR.gen <-  intervals.DMR.check(aineisto.DMR, DMRs.centromeric.nanopore, ind.cols = 5:42)

### Making a figure

#Gains and losses of centromeric cytosines
cent.cyt.events <- summarize(aineisto.cytosines, Gains = sum(Gains), Losses = sum(Losses))
#Gains and losses of centromeric DMRs
cent.DMR.events <- summarize(cent.DMR.gen, Gains = sum(Gains), Losses = sum(Losses))

events <- data.frame( count = unlist(c(cent.cyt.events, cent.DMR.events)), type = c("Gain", "Loss", "Gain", "Loss"), scale = c("Single cytosines", "Single cytosines", "DMRs", "DMRs"))

methylation.events <- ggplot(events, aes(y = count, x = type, fill = type)) +
    geom_bar(stat = "identity", colour = "black") +
    xlab("") +
    ylab("Count") +
    labs(fill = "") +    
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    scale_fill_discrete(labels = c("Gain" = "Gain of methylation", "Loss" = "Loss of methylation")) +
    facet_wrap(. ~ scale, scales = "free")

#Making the final plot
final.p <- plot_grid(p.WGBS, p.nano, methylation.events, ncol = 1, labels = c("A", "B", "C"), label_x = -0.03) + theme(plot.margin = unit(c(0,0,0,1), "cm"))

save_plot(filename = "./epigenet_ms/fig/methylation_prop.pdf", final.p, base_height = 9, base_width = 8)
### * Checking robustness of results to DMR spesifications

### ** Loading the WGBS DMR data with 100 bp, 200 bp, and 300 bp DMRs

### *** Generate the set of 200 bp DMRs

DMRs.all.CG <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_200bp/CG_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CG$context <- "CG"

DMRs.all.CHG <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_200bp/CHG_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHG$context <- "CHG"

DMRs.all.CHH <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_200bp/CHH_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHH$context <- "CHH"

DMRs.all <- rbind(DMRs.all.CG, DMRs.all.CHG, DMRs.all.CHH)
DMRs.all <- arrange(DMRs.all, seqnames, start)

#There are some bins that have been called in multiple contexts
bincheck <- paste0(DMRs.all[,1],DMRs.all[,2])
duplicatebins <- duplicated(bincheck)
#DMRs.all[duplicatebins,c(1:3,74)]

#Filter bins that are in different contexts (that is get only unique locations)
DMRs.all <- DMRs.all[!duplicatebins,]

### After combining different contexts there are some DMRs that should be combined into one (i.e. they are next to each other in the coordinates, e.g. CG and CHG DMRs that are next to each other) #TO DO!!!

#Filter everything that is not in the seven chromosomes
DMRs.all <- filter(DMRs.all, seqnames %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Load annotations
genes <- read.gff("~/Genomics/Neurospora/reference/Neurospora_crassa_NC12_genes.gff3")
genes <- filter(genes, type == "gene")
#We are only considering the seven chromosomes
genes <- filter(genes, seqid %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7") == T)
genes$seqid <- factor(genes$seqid) #Drop levels that are ot used

#Load TE annotations
TEs <- read.gff("~/Genomics/Neurospora/reference/Neurospora_crassa_NC12_TE.gff3")

#Promoter annotations, promdist the distance from TSS that still counts as promoter region
promoters <- annotate.promoters(genes, promdist = 500)

#This annotates all DMRs
DMRs.all <- annotate.DMRs(DMRs.all, genes, TEs, promoters) #Some nearby DMRs can occur in the same gene

## There are some warnings that occur
## Use options(warn=2), to convert warnings into errors and check what is the problem
## options(error = recover)
##

#Then we need information about the chromatin domains
#Load the data
#Euchromatin
euchr <- read.table("./annotation/2489.euchromatin.bed", header = F, sep = "\t")
colnames(euchr) <- c("Chromosome", "start", "end")
euchr <- filter(euchr, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Centromeric H3K9
centromeric.h3k9 <- read.table("./annotation/centromericH3K9.bed", header = F, sep = "\t")
colnames(centromeric.h3k9)[1:3] <- c("Chromosome", "start", "end")

#H3K9 not in centromeric regions
excent2.h3k9 <- read.table("./annotation/H3K9_ex_cent2.bed", header = F, sep = "\t")
colnames(excent2.h3k9)[1:3] <- c("Chromosome", "start", "end")
excent2.h3k9 <- filter(excent2.h3k9, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#H3K27 not in H3K9 regions
h3k27.exh3k9 <- read.table("./annotation/H3K27_exH3K9.bed", header = F, sep = "\t")
colnames(h3k27.exh3k9)[1:3] <- c("Chromosome", "start", "end")
h3k27.exh3k9 <- filter(h3k27.exh3k9, Chromosome %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#Check DMR overlap for chromatin domains
DMRs.all <- DMR.by.domains(DMRs.all, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)

##Then need to check overlap for each gene and promoter in the different domains
genes <- annotation.by.domains(genes, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)
TEs <- annotation.by.domains(TEs, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)
promoters <- annotation.by.domains(promoters, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9, prom = TRUE)

#Then check whether methylation was gained or lost
#Load the methylation level files
DMRs.all.CG.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_200bp/CG_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CG.meth$context <- "CG"

DMRs.all.CHG.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_200bp/CHG_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHG.meth$context <- "CHG"

DMRs.all.CHH.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_200bp/CHH_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHH.meth$context <- "CHH"

colnames(DMRs.all)[1] <- "Chromosome" #Need to change seqnames into "Chromosome"
#For each DMR check whether methylation was gained or lost
meth <- classify.DMRs(DMRs.all, ind.cols = 5:73,  DMRs.all.CG.meth, DMRs.all.CHG.meth, DMRs.all.CHH.meth)
DMRs.all <- cbind(DMRs.all, meth) #Store methylation levels, and type
DMRs.all.200bp <- DMRs.all

#Saving the data for easier loading the next time
save(DMRs.all.200bp, file = "./data/DMRs_all_200bp.RData")

### *** Generate the set of 300 bp DMRs

DMRs.all.CG <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_300bp/CG_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CG$context <- "CG"

DMRs.all.CHG <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_300bp/CHG_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHG$context <- "CHG"

DMRs.all.CHH <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_300bp/CHH_StateCalls-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHH$context <- "CHH"

DMRs.all <- rbind(DMRs.all.CG, DMRs.all.CHG, DMRs.all.CHH)
DMRs.all <- arrange(DMRs.all, seqnames, start)

#There are some bins that have been called in multiple contexts
bincheck <- paste0(DMRs.all[,1],DMRs.all[,2])
duplicatebins <- duplicated(bincheck)
#DMRs.all[duplicatebins,c(1:3,74)]

#Filter bins that are in different contexts (that is get only unique locations)
DMRs.all <- DMRs.all[!duplicatebins,]

### After combining different contexts there are some DMRs that should be combined into one (i.e. they are next to each other in the coordinates, e.g. CG and CHG DMRs that are next to each other) #TO DO!!!

#Filter everything that is not in the seven chromosomes
DMRs.all <- filter(DMRs.all, seqnames %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))

#This annotates all DMRs
DMRs.all <- annotate.DMRs(DMRs.all, genes, TEs, promoters) #Some nearby DMRs can occur in the same gene

#Check DMR overlap for chromatin domains
DMRs.all <- DMR.by.domains(DMRs.all, centromeric.h3k9, h3k27.exh3k9, excent2.h3k9)


#Then check whether methylation was gained or lost
#Load the methylation level files
DMRs.all.CG.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_300bp/CG_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CG.meth$context <- "CG"

DMRs.all.CHG.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_300bp/CHG_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHG.meth$context <- "CHG"

DMRs.all.CHH.meth <- read.table("~/Genomics/Neurospora/methylation/methimpute/jDMRmatrix/all_300bp/CHH_rcMethlvl-filtered-merged.txt", sep = "\t", header = T)
DMRs.all.CHH.meth$context <- "CHH"

colnames(DMRs.all)[1] <- "Chromosome" #Need to change seqnames into "Chromosome"
#For each DMR check whether methylation was gained or lost
meth <- classify.DMRs(DMRs.all, ind.cols = 5:73,  DMRs.all.CG.meth, DMRs.all.CHG.meth, DMRs.all.CHH.meth)
DMRs.all <- cbind(DMRs.all, meth) #Store methylation levels, and type
DMRs.all.300bp <- DMRs.all

#Saving the data for easier loading the next time
save(DMRs.all.300bp, file = "./data/DMRs_all_300bp.RData")


### *** Load the complete data
### 100 bp DMRs can be loaded with
load("./data/DMRs_all.RData") #DMRs.all
## Set of 200 bp
load("./data/DMRs_all_200bp.RData") #DMRs.all.200bp
## Set of 300 bp
load("./data/DMRs_all_300bp.RData") #DMRs.all.300bp


### ** Checking for DMR enrichment

#Then need the lenghts of these domains
domain.lengths <- c( sum(centromeric.h3k9$end - centromeric.h3k9$start) , sum(euchr$end - euchr$start), sum(h3k27.exh3k9$end - h3k27.exh3k9$start), sum(excent2.h3k9$end - excent2.h3k9$start) )

### 100bp DMRs
DMR.domains <- summarise(group_by(DMRs.all, domain), count = n())
DMR.domains$sequence <- domain.lengths

## Poisson model for DMR counts in the different domains...
#Reorder levels, so that euchromatin the one where everything else is compated to
DMR.domains$domain <- factor(DMR.domains$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27") )
#Expected counts, if DMRs would occur randomly across the genome
DMR.domains$expected <- sum(DMR.domains$count) * ( DMR.domains$sequence / sum(DMR.domains$sequence))

#Poisson model
#malli1 <- glm(count ~ offset(log(sequence)) + domain, family = poisson, data = DMR.domains)
malli1 <- glm(count ~ -1 + offset(log(expected)) + domain, family = poisson, data = DMR.domains)

#Format results
malli1.results <- cbind(coef(malli1), confint(malli1)) #Extract coefficients and conf int
malli1.results <- exp(malli1.results) #Change back to normal scale
malli1.results <- data.frame(malli1.results) #Change to dataframe
malli1.results$domain <- factor(c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
malli1.results$domain <- factor(malli1.results$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
colnames(malli1.results)[1:3] <- c("estimate", "lower", "upper")


#Plotting
ggplot(malli1.results, aes( x = domain, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange(pch = 1) +
    geom_hline(yintercept = 1, lty = "dashed") +
    xlab("") +
    ylab("Relative to expected") +
    scale_y_log10() +
    scale_x_discrete(guide = guide_axis(angle = 45))   


### 200bp DMRs
DMR.domains.200 <- summarise(group_by(DMRs.all.200bp, domain), count = n())
DMR.domains.200$sequence <- domain.lengths

#Reorder levels, so that euchromatin the one where everything else is compated to
DMR.domains.200$domain <- factor(DMR.domains.200$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27") )
#Expected counts, if DMRs would occur randomly across the genome
DMR.domains.200$expected <- sum(DMR.domains.200$count) * ( DMR.domains.200$sequence / sum(DMR.domains.200$sequence))

#Poisson model
#malli1 <- glm(count ~ offset(log(sequence)) + domain, family = poisson, data = DMR.domains)
malli1.200 <- glm(count ~ -1 + offset(log(expected)) + domain, family = poisson, data = DMR.domains.200)

#Format results
malli1.results.200 <- cbind(coef(malli1.200), confint(malli1.200)) #Extract coefficients and conf int
malli1.results.200 <- exp(malli1.results.200) #Change back to normal scale
malli1.results.200 <- data.frame(malli1.results.200) #Change to dataframe
malli1.results.200$domain <- factor(c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
malli1.results.200$domain <- factor(malli1.results.200$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
colnames(malli1.results.200)[1:3] <- c("estimate", "lower", "upper")

### 300bp DMRs
DMR.domains.300 <- summarise(group_by(DMRs.all.300bp, domain), count = n())
DMR.domains.300$sequence <- domain.lengths

#Reorder levels, so that euchromatin the one where everything else is compated to
DMR.domains.300$domain <- factor(DMR.domains.300$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27") )
#Expected counts, if DMRs would occur randomly across the genome
DMR.domains.300$expected <- sum(DMR.domains.300$count) * ( DMR.domains.300$sequence / sum(DMR.domains.300$sequence))

#Poisson model
#malli1 <- glm(count ~ offset(log(sequence)) + domain, family = poisson, data = DMR.domains)
malli1.300 <- glm(count ~ -1 + offset(log(expected)) + domain, family = poisson, data = DMR.domains.300)

#Format results
malli1.results.300 <- cbind(coef(malli1.300), confint(malli1.300)) #Extract coefficients and conf int
malli1.results.300 <- exp(malli1.results.300) #Change back to normal scale
malli1.results.300 <- data.frame(malli1.results.300) #Change to dataframe
malli1.results.300$domain <- factor(c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
malli1.results.300$domain <- factor(malli1.results.300$domain, levels = c("Euchromatin", "Centromeric", "H3K9", "H3K27"))
colnames(malli1.results.300)[1:3] <- c("estimate", "lower", "upper")

#Saving the results

DMR.window.results <- rbind(malli1.results, malli1.results.200, malli1.results.300)
DMR.window.results$specification <- c(rep("100 bp, 5 mC", 4), rep("200 bp, 10 mC", 4), rep("300 bp, 12 mC",4))

save(DMR.window.results, file = "./data/DMR_enrich_windows.RData")



### ** Divergence in centromeric regions

### Centromeric 100 bp DMRs
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K9cent_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K9cent_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/AB_output_jDMR_H3K9cent_CHH.RData")

#AB.output.jDMR.H3K9cent.CG is the name of the variable where divergence data is stored

#Final pedigree CG
pedigree.CG <- AB.output.jDMR.H3K9cent.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.jDMR.H3K9cent.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.jDMR.H3K9cent.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.jDMR.H3K9cent.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.jDMR.H3K9cent.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.jDMR.H3K9cent.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Plotting
#Make plot dataframe
plotdata.DMR.cent.100bp <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))


### Centromeric 200 bp DMRs
load(file = "./pedigree_data/nanopore_ABdata/centromeric_200bp/AB_output_jDMR_centromeric_200bp_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/centromeric_200bp/AB_output_jDMR_centromeric_200bp_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/centromeric_200bp/AB_output_jDMR_centromeric_200bp_CHH.RData")

#AB.output.jDMR.centromeric.200bp.CG is the name of the variable where divergence data is stored

#Final pedigree CG
pedigree.CG <- AB.output.jDMR.centromeric.200bp.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.jDMR.centromeric.200bp.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.jDMR.centromeric.200bp.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.jDMR.centromeric.200bp.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.jDMR.centromeric.200bp.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.jDMR.centromeric.200bp.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Make plot dataframe
plotdata.DMR.cent.200bp <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))


### Centromeric 300 bp DMRs
load(file = "./pedigree_data/nanopore_ABdata/centromeric_300bp/AB_output_jDMR_centromeric_300bp_CG.RData")
load(file = "./pedigree_data/nanopore_ABdata/centromeric_300bp/AB_output_jDMR_centromeric_300bp_CHG.RData")
load(file = "./pedigree_data/nanopore_ABdata/centromeric_300bp/AB_output_jDMR_centromeric_300bp_CHH.RData")

#AB.output.jDMR.centromeric.200bp.CG is the name of the variable where divergence data is stored

#Final pedigree CG
pedigree.CG <- AB.output.jDMR.centromeric.300bp.CG$Pdata
pedigree.CG[,1:3] <- pedigree.CG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CG[,2] + pedigree.CG[,3] - 2*pedigree.CG[,1]
p0uu_in.CG <- AB.output.jDMR.centromeric.300bp.CG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CG[,1])
pedigree.final.CG <- pedigree.CG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHG
pedigree.CHG <- AB.output.jDMR.centromeric.300bp.CHG$Pdata
pedigree.CHG[,1:3] <- pedigree.CHG[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHG[,2] + pedigree.CHG[,3] - 2*pedigree.CHG[,1]
p0uu_in.CHG <- AB.output.jDMR.centromeric.300bp.CHG$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHG[,1])
pedigree.final.CHG <- pedigree.CHG[tokeep,]
dt.final <- dt[tokeep]

#Final pedigree CHH
pedigree.CHH <- AB.output.jDMR.centromeric.300bp.CHH$Pdata
pedigree.CHH[,1:3] <- pedigree.CHH[,1:3]*25 #Convert transfers to mitoses
dt <- pedigree.CHH[,2] + pedigree.CHH[,3] - 2*pedigree.CHH[,1]
p0uu_in.CHH <- AB.output.jDMR.centromeric.300bp.CHH$tmpp0
#Need to remove comparisons to the other mating type
tokeep <- !is.infinite(pedigree.CHH[,1])
pedigree.final.CHH <- pedigree.CHH[tokeep,]
dt.final <- dt[tokeep]

#Make plot dataframe
plotdata.DMR.cent.300bp <- data.frame(rbind(pedigree.final.CG, pedigree.final.CHG, pedigree.final.CHH), dt = rep(dt.final, 3), context = c(rep("CG", nrow(pedigree.final.CG)), rep("CHG", nrow(pedigree.final.CHG)), rep("CHH", nrow(pedigree.final.CHH))))


DMR.window.results$specification <- c(rep("100 bp, 5 mC", 4), rep("200 bp, 10 mC", 4), rep("300 bp, 12 mC",4))

#Combining the divergence results
plotdata.DMR.cent.windows <- rbind(plotdata.DMR.cent.100bp, plotdata.DMR.cent.200bp, plotdata.DMR.cent.300bp)
plotdata.DMR.cent.windows$specification <- c(rep("100 bp, 5 mC", nrow(plotdata.DMR.cent.100bp)), rep("200 bp, 10 mC", nrow(plotdata.DMR.cent.200bp)), rep("300 bp, 12 mC",nrow(plotdata.DMR.cent.300bp)))

#Save the results for plotting
save(plotdata.DMR.cent.windows, file = "./data/DMR_divergence_windows.RData")

### ** Plotting the final results

#Loading the WGBS DMR results
load(file = "./data/DMR_enrich_windows.RData")

#Loading the Nanopore DMR centromeric divergence results
load(file = "./data/DMR_divergence_windows.RData")

#Plotting
p.enrich <- ggplot(DMR.window.results, aes( x = domain, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange(pch = 1) +
    geom_hline(yintercept = 1, lty = "dashed") +
    xlab("") +
    ylab("Relative to expected") +
    scale_y_log10() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    facet_wrap(. ~ specification)


myxlab <- TeX("$\\Delta$t (mitoses)")
p.cent.div <- ggplot(plotdata.DMR.cent.windows, aes(x = dt, y = D.value)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "lm") +
    #ggtitle("DMRs, Nanopore") +
    xlab(myxlab) +
    ylab("Methylation divergence") +
    facet_grid(context ~ specification) #+
    #theme(text=element_text(size=16), axis.text = element_text(size = 16))

p.DMR.window <- plot_grid(p.enrich, p.cent.div, nrow = 2, labels = c("A", "B"), rel_heights = c(1.5,3))

save_plot(filename = "./epigenet_ms/fig/DMR_windows.pdf", p.DMR.window, base_height = 10, base_width = 7)
### * Checking DMRs in Nanopore dataset, regarding nucleosome drift

### Load the data
load("./data/DMRs_all_nanopore.RData")

#Folder for the methimpute data
datafolder <- "~/Genomics/Neurospora/methylation/methimpute/nanopore_centromeric/"
samples <- list.files(datafolder) #These are the MA ancestor samples
samples <- samples[-c(2:4,6:8)] #Drop samples that are not going to used for the analysis

#Need to write a function that loads a particular region for a set of samples
load.methylome.region.MA <- function(chr, reg.start, reg.end, samples) {
    nsamples <- length(samples) #Number of samples
    region <- NULL
    mycols <- c("seqnames", "start", "posteriorMax", "rcmethlvl")
    for(i in 1:nsamples) { #Loop over all samples
        sample <- fread(file = paste0(datafolder, samples[i]))
        sample <- sample[, ..mycols] #Select only the columns that are needed
        sample <- sample[seqnames == chr & start > reg.start & start < reg.end] #Filter for the region
        sample$sample <- gsub(".txt", "", samples[i]) #Store sample name
        #Making variables transfer and line
        if(grepl("MatA", samples[i]) == T) { sample$transfer <- "Transfer 0"; sample$line <- "Ancestor mat A" }
        if(grepl("Mata", samples[i]) == T) { sample$transfer <- "Transfer 0"; sample$line <- "Ancestor mat a" }
        if(grepl("MatA", samples[i]) == F & grepl("Mata", samples[i]) ==F) {
            sample$line <- paste("Line", str_extract(samples[i], "(?<=L)\\d+"), sep = " ")
            sample$transfer <- paste("Transfer", str_extract(samples[i], "(?<=G)\\d+"), sep = " ")
            }                                                                   
        region <- rbind(region, sample) #
    } #Done looping over all samples
    return(region)
    } #Done


#Filtering for the centromeric DMRs
DMRs.cent.nano <- filter(DMRs.all.nanopore, domain == "Centromeric") #Some 10 000 DMRs left

DMRs.cent.nano[9,1:42]

#Should probably do some sort of classification of of DMRs per line basis?
test <- grepl("L11|MatA", colnames(DMRs.cent.nano))
test[1:4] <- TRUE

samples[grepl("L11|MatA", samples)]

DMRs.cent.nano[13,1:42]

tail(filter(DMRs.cent.nano, Chromosome == "Supercontig_12.2", DMR.type == "gain"), n = 15)

### Locus chr 1, 3687501
loc1 <- load.methylome.region.MA(chr = "Supercontig_12.1", reg.start = 3687400, reg.end = 3687700, samples = samples[grepl("L2|MatA", samples)])
loc1$transfer <- factor(loc1$transfer, levels = c("Transfer 0", "Transfer 1", "Transfer 5", "Transfer 7", "Transfer 8", "Transfer 10", "Transfer 15"))

ggplot(loc1, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        ylab("Methylation (%)") +
        xlab("Position")  +
        scale_y_continuous( expand = c(0.0,0)) +
        facet_grid(transfer ~ .)

### This is a potential demonstration
loc2 <- load.methylome.region.MA(chr = "Supercontig_12.1", reg.start = 3687400, reg.end = 3688000, samples = samples[grepl("L23|Mata", samples)])
loc2$transfer <- fct_recode(loc2$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc2$transfer <- factor(loc2$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

ggplot(loc2, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        annotate("rect", x = NULL, xmin = 3687600, xmax = 3687800, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)

### Another locus (I guess this could be compatible with nucleosome drift?)
loc3 <- load.methylome.region.MA(chr = "Supercontig_12.1", reg.start = 3689200, reg.end = 3689700, samples = samples[grepl("L11|MatA", samples)])
loc3$transfer <- fct_recode(loc3$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc3$transfer <- factor(loc3$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

ggplot(loc3, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        annotate("rect", x = NULL, xmin = 3689455, xmax = 3689526, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)

### This locus looks quite interesting
loc4 <- load.methylome.region.MA(chr = "Supercontig_12.1", reg.start = 3714100, reg.end = 3715000, samples = samples[grepl("L5|MatA", samples)])
loc4$transfer <- fct_recode(loc4$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc4$transfer <- factor(loc4$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

p.loc4 <- ggplot(loc4, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        annotate("rect", x = NULL, xmin = 3714312, xmax = 3714525, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        annotate("rect", x = NULL, xmin = 3714565, xmax = 3714645, y = NULL, ymin = 0, ymax = 20, fill = "orange", alpha = 0.2) +    
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        ggtitle("Line 5, chr 1: 3714100 - 3715000") +     
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)


##Supercontig_12.7 2495101 2495200
loc5 <- load.methylome.region.MA(chr = "Supercontig_12.7", reg.start = 2495100, reg.end = 2495500, samples = samples[grepl("L11|MatA", samples)])
loc5$transfer <- fct_recode(loc5$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc5$transfer <- factor(loc5$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

p.loc5 <- ggplot(loc5, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        annotate("rect", x = NULL, xmin = 2495185, xmax = 2495242, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        annotate("rect", x = NULL, xmin = 2495250, xmax = 2495341, y = NULL, ymin = 0, ymax = 20, fill = "orange", alpha = 0.2) +    
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        ggtitle("Line 11, chr 7: 2495100 - 2495500") +     
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)


##Supercontig_12.3 967301 968000 #This gain of methylation is not compatible with nucleosome drift!
loc6 <- load.methylome.region.MA(chr = "Supercontig_12.3", reg.start = 966500, reg.end = 968100, samples = samples[grepl("L5|MatA", samples)])
loc6$transfer <- fct_recode(loc6$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc6$transfer <- factor(loc6$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

p.loc6 <- ggplot(loc6, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        annotate("rect", x = NULL, xmin = 966800, xmax = 967870, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        #annotate("rect", x = NULL, xmin = 2495250, xmax = 2495341, y = NULL, ymin = 0, ymax = 20, fill = "orange", alpha = 0.2) +    
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        ggtitle("Line 5, chr 3: 966500 - 968100") +     
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)



head(filter(DMRs.cent.nano, Chromosome == "Supercontig_12.2", DMR.type == "gain"), n = 35)


loc7 <- load.methylome.region.MA(chr = "Supercontig_12.2", reg.start = 1355500, reg.end = 1359800, samples = samples[grepl("L23|Mata", samples)])
loc7$transfer <- fct_recode(loc7$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc7$transfer <- factor(loc7$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

ggplot(loc7, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        #annotate("rect", x = NULL, xmin = 966800, xmax = 967870, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        #annotate("rect", x = NULL, xmin = 2495250, xmax = 2495341, y = NULL, ymin = 0, ymax = 20, fill = "orange", alpha = 0.2) +    
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        #ggtitle("Line 5, chr 3: 966500 - 968100") +     
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)


#This is one possibility
loc7 <- load.methylome.region.MA(chr = "Supercontig_12.2", reg.start = 1172401, reg.end = 1175000, samples = samples[grepl("L23|Mata", samples)])
loc7$transfer <- fct_recode(loc7$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc7$transfer <- factor(loc7$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

ggplot(loc7, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        #annotate("rect", x = NULL, xmin = 966800, xmax = 967870, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        #annotate("rect", x = NULL, xmin = 2495250, xmax = 2495341, y = NULL, ymin = 0, ymax = 20, fill = "orange", alpha = 0.2) +    
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        #ggtitle("Line 5, chr 3: 966500 - 968100") +     
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)


#One possibility
loc8 <- load.methylome.region.MA(chr = "Supercontig_12.2", reg.start = 1174001, reg.end = 1179800, samples = samples[grepl("L23|Mata", samples)])
loc8$transfer <- fct_recode(loc8$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc8$transfer <- factor(loc8$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

ggplot(loc8, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        #annotate("rect", x = NULL, xmin = 966800, xmax = 967870, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        #annotate("rect", x = NULL, xmin = 2495250, xmax = 2495341, y = NULL, ymin = 0, ymax = 20, fill = "orange", alpha = 0.2) +    
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        #ggtitle("Line 5, chr 3: 966500 - 968100") +     
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)

filter(DMRs.cent.nano, Chromosome == "Supercontig_12.6", DMR.type == "gain")[60:65,]


#reg.start = 1179801, reg.end = 1180500
loc9 <- load.methylome.region.MA(chr = "Supercontig_12.2", reg.start = 1176001, reg.end = 1185020, samples = samples[grepl("L25|Mata", samples)])
loc9$transfer <- fct_recode(loc9$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc9$transfer <- factor(loc9$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

ggplot(loc9, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        #annotate("rect", x = NULL, xmin = 966800, xmax = 967870, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        #annotate("rect", x = NULL, xmin = 2495250, xmax = 2495341, y = NULL, ymin = 0, ymax = 20, fill = "orange", alpha = 0.2) +    
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        #ggtitle("Line 5, chr 3: 966500 - 968100") +     
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)


loc10 <- load.methylome.region.MA(chr = "Supercontig_12.6", reg.start = 2872000, reg.end = 2873900, samples = samples[grepl("L25|Mata", samples)])
loc10$transfer <- fct_recode(loc10$transfer, "0" = "Transfer 0", "1" = "Transfer 1", "5" = "Transfer 5", "7" = "Transfer 7", "8" = "Transfer 8", "10" = "Transfer 10", "15" = "Transfer 15")
loc10$transfer <- factor(loc10$transfer, levels = c("0", "1", "5", "7", "8", "10", "15"))

p.loc10 <- ggplot(loc10, aes(y = rcmethlvl*100, x = start)) +
        #geom_line(alpha = 0.8) +
        geom_bar(stat = "identity") +
        annotate("rect", x = NULL, xmin = 2872000, xmax = 2873900, y = NULL, ymin = 0, ymax = 20, fill = "blue", alpha = 0.2) +
        #annotate("rect", x = NULL, xmin = 2495250, xmax = 2495341, y = NULL, ymin = 0, ymax = 20, fill = "orange", alpha = 0.2) +    
        ylab("Methylation (%)") +
        xlab("Position (bp)")  +
        ggtitle("Line 23, chr 6: 2872000 - 2873900") +     
        scale_y_continuous( expand = c(0.0,0), breaks = c(0,20)) +
        facet_grid(transfer ~ .)


#DMR figures
MA.DMR.plot <- plot_grid(p.loc6, p.loc10, p.loc4, p.loc5, ncol = 2, labels = c("A", "B", "C", "D"))

save_plot(filename = "./epigenet_ms/fig/MA.DMR.plot.pdf", MA.DMR.plot, base_height = 7, base_width = 12)

#Make on figure with two examples of DMR changes not compatible with any kind of nucleosome drift:
#loc6
#loc10 (methylation changes happen over a longer regions, losses and gains)

#Make a figure with two examples of methylated regions expanding into unmethylated region
#loc4
#loc5
