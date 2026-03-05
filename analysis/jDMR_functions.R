#Functions dealing with jDMR data
#

#This functions uses gene annotation to make promoter annotations
#promdist is the assumed length of the promoters
annotate.promoters <- function(genes, promdist) {

genes$p.type <- rep("unidirectional", nrow(genes)) #Initialize promoter type vector
promoters <- NULL #Initialize promoters

#Need to check for bidirectional promoters
for(i in 1:length(levels(genes$seqid))) {
    current.chr <- filter(genes, seqid == levels(genes$seqid)[i]) #Take genes from current chromosome
    current.chr$gene1 <- rep(NA, nrow(current.chr))
    current.chr$gene1.type <- rep(NA, nrow(current.chr))
    current.chr$gene1.desc <- rep(NA, nrow(current.chr))
    current.chr$gene2 <- rep(NA, nrow(current.chr))
    current.chr$gene2.type <- rep(NA, nrow(current.chr))
    current.chr$gene2.desc <- rep(NA, nrow(current.chr))
    current.chr$p.start <- rep(NA, nrow(current.chr))
    current.chr$p.end <- rep(NA, nrow(current.chr))
                               

    #Check if genes have a potentially bidirectional promoter
    j <- 1 #initialize row counter
    while(j < nrow(current.chr)) {
        if(current.chr[j,7] == "-" & current.chr[j+1,7] == "+" & (current.chr[j+1,4] - current.chr[j,5] <= promdist)   ) {
            current.chr$p.type[j] <- "bidirectional"
            current.chr$p.type[j+1] <- "bidirectional"
            #Then store information about the genes for bidirectional promoters
            gene1 <- unlist(strsplit(current.chr$attributes[j], ";"))
            gene2 <- unlist(strsplit(current.chr$attributes[j+1], ";"))
            current.chr$gene1[j] <- unlist(strsplit(gene1[1],":"))[2]
            current.chr$gene2[j] <- unlist(strsplit(gene2[1],":"))[2]
            current.chr$gene1.type[j] <- unlist(strsplit(gene1[2],"="))[2]
            current.chr$gene2.type[j] <- unlist(strsplit(gene2[2],"="))[2]
            current.chr$gene1.desc[j] <- unlist(strsplit(gene1[3],"="))[2]
            current.chr$gene2.desc[j] <- unlist(strsplit(gene2[3],"="))[2]
            current.chr$p.start[j] <- current.chr$end[j] + 1 #Start coordinate of the promoter
            current.chr$p.end[j] <- current.chr$start[j+1] - 1 #End coordinate of the promoter
            #Need to deal with overlapping genes
            if(current.chr$p.end[j] - current.chr$p.start[j] < 0) { #If genes overlap set promoters as unidirectional
                current.chr$p.type[j] <- "unidirectional"
                current.chr$p.type[j+1] <- "overlapping"
                current.chr$p.start[j] <- current.chr$end[j] + 1
                current.chr$p.end[j] <- current.chr$end[j] + 1 + promdist
            } #Note, this skips the overlapping gene, but, DMRs in the promoter of that gene should be in this case annotated as genic anyway, so this should be fine
                
            j <- j + 2 } else { j <- j + 1 }
    }
    #browser() #Does not work, there is some problems...
    #Then take all bidirectional promoters, filter out every other row which is not used
    pro.bidirectional <- filter(current.chr, p.type == "bidirectional") #Filter for bidirectional entries
    #Need to filter out overlapping genes #TO DO!
    
    #Take every other element
    take <- seq(1, nrow(pro.bidirectional), by = 2)
    pro.bidirectional <- pro.bidirectional[take,]
    

    #Fix start and end coordinates for the unidirectional promoters
    pro.unidirectional <- filter(current.chr, p.type == "unidirectional")

    #Promoters on the plus strand
    promoters.plus <- filter(pro.unidirectional, strand == "+")
        
    #Start and end coordinates of the promoters
    promoters.plus$p.start <- promoters.plus$start - 1 - promdist
    promoters.plus$p.end <- promoters.plus$start - 1

    #Need to take attributes etc. from lists etc. vectorised
    myattr <- strsplit(promoters.plus$attributes, ";")
    k1 <- lapply(myattr, '[[', 1)
    k2 <- strsplit(unlist(k1), ":")
    promoters.plus$gene1 <- unlist(lapply(k2, '[[', 2))
    t1 <- lapply(myattr, '[[', 2)
    t2 <- strsplit(unlist(t1), "=")
    promoters.plus$gene1.type <- unlist(lapply(t2, '[[', 2))
    d1 <- lapply(myattr, '[[', 3)
    d2 <- strsplit(unlist(d1), "=")
    promoters.plus$gene1.desc <- unlist(lapply(d2, '[[', 2))

    #Promoters on the minus strand
    promoters.minus <- filter(pro.unidirectional, strand == "-")
    
    #Note that for promoters in the "-" strand promoter coordinates are given starting from promoter end
    promoters.minus$p.start <- promoters.minus$end + 1
    promoters.minus$p.end <- promoters.minus$end + 1 + promdist

    #Need to take attributes etc. from lists etc. vectorised
    myattr <- strsplit(promoters.minus$attributes, ";")
    k1 <- lapply(myattr, '[[', 1)
    k2 <- strsplit(unlist(k1), ":")
    promoters.minus$gene1 <- unlist(lapply(k2, '[[', 2))
    t1 <- lapply(myattr, '[[', 2)
    t2 <- strsplit(unlist(t1), "=")
    promoters.minus$gene1.type <- unlist(lapply(t2, '[[', 2))
    d1 <- lapply(myattr, '[[', 3)
    d2 <- strsplit(unlist(d1), "=")
    promoters.minus$gene1.desc <- unlist(lapply(d2, '[[', 2)) 
    

    #Combine all promoters
    promoter.current.chr <- rbind(pro.bidirectional, promoters.plus, promoters.minus)
    promoter.current.chr <- arrange(promoter.current.chr, start)

    #browser()
    #Store promoters
    promoters <- rbind(promoters, promoter.current.chr)
   
    } #Done looping over chromosomes

 #Calculate promoter lengths
 promoters$p.length <- promoters$p.end - promoters$p.start #Calculate promoter lengths

return(promoters) #Return results

} #Done

#Function to check overlap of DMRs and annotations
#There is something weird happening in chipseq peak overlap, check this again!!! TO DO!
check.overlap <- function(x, start.dmr, end.dmr) {
            start.feat <- x[1] #start of feature
            end.feat <- x[2] #end of feature
            ifelse( ( start.dmr < start.feat & end.dmr > start.feat ) | (start.dmr >= start.feat & end.dmr <= end.feat ) | ( start.dmr < end.feat & end.dmr > end.feat ) , TRUE, FALSE) }

#Function to annotate DMRs
annotate.DMRs <- function(DMRs, genes, TEs, promoters) {
    #Loop over all DMRs and check for annotations
    #If DMR within 500 (promdist) of transcription start site, then considered a promoter DMR
    #otherwise intergenic
    nDMR <- nrow(DMRs) #Number of DMRs
    anno.type <- rep("intergenic", nDMR)
    TE.type <- rep(NA, nDMR)
    TE.class <- rep(NA, nDMR)
    gene.id <- rep(NA, nDMR)
    gene.type <- rep(NA, nDMR)
    gene.description <- rep(NA, nDMR)
    p.type <- rep(NA, nDMR)
    gene2.id <- rep(NA, nDMR)
    gene2.type <- rep(NA, nDMR)
    gene2.desc <- rep(NA, nDMR)
    p.length <- rep(NA, nDMR)

    

    #Start the loop
    for(i in 1:nDMR) {
        chr <- as.character(DMRs$seqnames[i])
        DMR.start <- DMRs$start[i]
        DMR.end <- DMRs$end[i]
        
        #First check TEs
        TE.test <- check.overlap(filter(TEs, seqid == chr)[,4:5], DMR.start, DMR.end)
        #apply(filter(TEs, seqid == chr)[,4:5], MARGIN = 1, check.overlap, start.dmr = DMR.start, end.dmr = DMR.end)

        #If DMR overlaps a TE
        if(any(TE.test)) {
            anno.type[i] <- "TE"
            myTE <- filter(TEs, seqid == chr)[TE.test,]
            myattr <- unlist(strsplit(myTE$attributes, ";")) #Get info about TE type and class
            TE.type[i] <- unlist(strsplit(myattr[1],"'"))[2] #Extract the information
            TE.class[i] <- unlist(strsplit(myattr[2],"'"))[2]
            #browser() 
        }

        #If neither gene nor TE, check if promoter
        promoter.test <- check.overlap(filter(promoters, seqid == chr)[,17:18], DMR.start, DMR.end)

        #In some cases a DMR can span two promoters, if two genes are not classified as bidirectional promoter, but they are in the same orientation.
        

        if(any(promoter.test)) {
               #browser()
               mypromoter <- filter(promoters, seqid == chr)[promoter.test,]
               anno.type[i] <- "promoter"
               #Cases where two genes are not classified as bidirectional
               #but they are close enough so DMR overlaps two promoters
               if(nrow(mypromoter) == 2) {
               gene.id[i] <- mypromoter$gene1[1]
               gene.type[i] <- mypromoter$gene1.type[1]
               gene.description[i] <- mypromoter$gene1.desc[1]
               p.type[i] <- mypromoter$p.type[1]
               gene2.id[i] <- mypromoter$gene2[2]
               gene2.type[i] <- mypromoter$gene2.type[2]
               gene2.desc[i] <- mypromoter$gene2.desc[2]
               p.length[i] <- mypromoter$p.length[1]
               
               } else {
               gene.id[i] <- mypromoter$gene1
               gene.type[i] <- mypromoter$gene1.type
               gene.description[i] <- mypromoter$gene1.desc
               p.type[i] <- mypromoter$p.type
               gene2.id[i] <- mypromoter$gene2
               gene2.type[i] <- mypromoter$gene2.type
               gene2.desc[i] <- mypromoter$gene2.desc
               p.length[i] <- mypromoter$p.length
               }
           }

        #Check of DMR overlaps a gene
        gene.test <- check.overlap(filter(genes, seqid == chr)[,4:5], DMR.start, DMR.end)

        #If DMR overlaps a gene
        if(any(gene.test)) {
            anno.type[i] <- "gene"
            #Find the overlapping gene
            mygene <- filter(genes, seqid == chr)[gene.test,]
            myattr <- unlist(strsplit(mygene$attributes, ";"))
            gene.id[i] <- unlist(strsplit(myattr[1],":"))[2]
            gene.type[i] <- unlist(strsplit(myattr[2],"="))[2]
            gene.description[i] <- unlist(strsplit(myattr[3],"="))[2]
            #browser()
        }

        #if(anno.type[i] == "intergenic") {
        #    
        #    #browser
        #    ##TO DO!!!
        #}
        
    } #End loop

    #browser()
    DMRs$anno.type <- anno.type #Store annotation types
    DMRs$TE.type <- TE.type
    DMRs$TE.class <- TE.class
    DMRs$gene.id <- gene.id
    DMRs$gene.type <- gene.type
    DMRs$gene.description <- gene.description
    DMRs$p.type <- p.type
    DMRs$gene2.id <- gene2.id
    DMRs$gene2.type <- gene2.type
    DMRs$gene2.desc <- gene2.desc
    DMRs$p.length <- p.length

    return(DMRs)
} #End function



##How many DMRs in different regions?
DMR.by.domains <- function(DMRs, cent, h3k27.plot, h3k9.plot) {
    domain <- rep("", nrow(DMRs)) #initialize results vector
    #colnames(cent) <- c("Chromosome", "Start", "End") #Fix colnames for centromeric regions
    #Loop over all DMRs
    for(i in 1:nrow(DMRs)) {
        #What is the current chromosome?
        cur.chr <- as.character(DMRs[i,1])
        #Building the domain table for the current chromosome
        cur.cent <- filter(cent, Chromosome == cur.chr)
        cur.h3k27 <- filter(h3k27.plot, Chromosome == cur.chr)
        cur.h3k9 <- filter(h3k9.plot, Chromosome == cur.chr)
        cur.domains <- rbind(cur.cent[,1:3], cur.h3k27[,1:3], cur.h3k9[,1:3])
        cur.domains$domain <- c( rep("Centromeric", nrow(cur.cent)), rep("H3K27", nrow(cur.h3k27)), rep("H3K9", nrow(cur.h3k9)) )
        check <- DMRs[i,]$start >= cur.domains$start & DMRs[i,]$end <= cur.domains$end
        domain.belong <- cur.domains$domain[check]
        if(length(domain.belong) == 0) { domain[i] <- "Euchromatin" }
        if(length(domain.belong) >=1 & any(grepl("Centromeric", domain.belong))) { domain[i] <- "Centromeric" }
        if(length(domain.belong) >=1 & any(grepl("H3K9", domain.belong)) == T & any(grepl("Centromeric", domain.belong)) == F ) { domain[i] <- "H3K9" }
        if(length(domain.belong) >=1 & any(grepl("H3K9", domain.belong)) == F & any(grepl("H3K27", domain.belong)) == T) { domain[i] <- "H3K27" }
    }
    return(cbind(DMRs,domain))
}

annotation.by.domains <- function(annotation, cent, h3k27, h3k9, prom = FALSE) {
    domain <- rep("", nrow(annotation)) #initialize results vector
    #Loop over annotations
    for(i in 1:nrow(annotation)) {
            #What is the current chromosome?
            cur.chr <- as.character(annotation[i,1])
            #Building the domain table for the current chromosome
            cur.cent <- filter(cent, Chromosome == cur.chr)
            cur.h3k27 <- filter(h3k27, Chromosome == cur.chr)
            cur.h3k9 <- filter(h3k9, Chromosome == cur.chr)
            cur.domains <- rbind(cur.cent[,1:3], cur.h3k27[,1:3], cur.h3k9[,1:3])
            cur.domains$domain <- c( rep("Centromeric", nrow(cur.cent)), rep("H3K27", nrow(cur.h3k27)), rep("H3K9", nrow(cur.h3k9)))
            if(prom == TRUE) { check <- annotation[i,]$p.start >= cur.domains$start & annotation[i,]$p.end <= cur.domains$end } else {
            check <- annotation[i,]$start >= cur.domains$start & annotation[i,]$end <= cur.domains$end }
            domain.belong <- cur.domains$domain[check]
            if(length(domain.belong) == 0) { domain[i] <- "Euchromatin" }
            if(length(domain.belong) >=1 & any(grepl("Centromeric", domain.belong))) { domain[i] <- "Centromeric" }
            if(length(domain.belong) >=1 & any(grepl("H3K9", domain.belong)) == T & any(grepl("Centromeric", domain.belong)) == F ) { domain[i] <- "H3K9" }
            if(length(domain.belong) >=1 & any(grepl("H3K9", domain.belong)) == F & any(grepl("H3K27", domain.belong)) == T) { domain[i] <- "H3K27" }
        }
    return(cbind(annotation,domain))
}


#This function check whether a DMR is a gain or loss of methylation
#DMR.gain.loss <- function(x) { ifelse(sum(x) > 34, "loss", "gain") }

###This function loops over all DMRs and check methylation levels of 1 and 0 alleles
###It also classifies the DMR into a loss or gain of methylation
classify.DMRs <- function(DMRs, ind.cols, DMR.CG.meth, DMR.CHG.meth, DMR.CHH.meth) {
    ndmr <- nrow(DMRs) #Number of DMRs
    meth.1 <- rep(0, ndmr) #methylation level of 1 alleles
    meth.0 <- rep(0, ndmr) #methylation level of 0 alleles
    DMR.type <- rep("", ndmr) #DMR class
    DMR.common <- rep(0, ndmr) #Which allele is the most common?

    #The individual genotypes / methylation levels are in columns, specified by ind.cols (old 5:73)

    #Loop over all DMRs
    for(i in 1:ndmr) {
        #First find the methylation levels for the current DMR
        #Need to find the context and coordinates of the current DMR and match that with methylation levels
        cur.dmr <- DMRs[i,]
        cur.context <- cur.dmr$context
        cur.chr <- cur.dmr$Chromosome
        cur.start <- cur.dmr$start
        cur.end <- cur.dmr$end
        alleles <- cur.dmr[,ind.cols] #Get the state calls from current DMR
        #Check context and store allele information
        if(cur.context == "CG") {
            #Context is CG
            index <- which(DMR.CG.meth$seqnames == as.character(cur.chr) & DMR.CG.meth$start == cur.start)
            meth.alleles <- DMR.CG.meth[index,ind.cols]
        }
        if(cur.context == "CHG") {
            #Context is CHG
            #Get the index of the DMR from the corresponding methylation 
            index <- which(DMR.CHG.meth$seqnames == as.character(cur.chr) & DMR.CHG.meth$start == cur.start)
            meth.alleles <- DMR.CHG.meth[index,ind.cols]
        }
        if(cur.context == "CHH") {
            #Context is CHH
            #Get the index of the DMR from the corresponding methylation 
            index <- which(DMR.CHH.meth$seqnames == as.character(cur.chr) & DMR.CHH.meth$start == cur.start)
            meth.alleles <- DMR.CHH.meth[index,ind.cols]   
        }

        #Get median methylation of the two different allele classes
        meth.1[i] <- median(meth.alleles[alleles == 1]) #Methylation level of allele 1
        meth.0[i] <- median(meth.alleles[alleles == 0]) #Methylation level of allele 0
        #Check which one is the more common
        allele.counts <- table(as.numeric(alleles))
        DMR.common[i] <- as.numeric(names(allele.counts)[allele.counts == max(allele.counts)])
        #Check has DNA methylation been gained or lost
        #browser()
        if(DMR.common[i] == 1) { DMR.type[i] <- ifelse(meth.1[i] > meth.0[i], "loss", "gain") }
        if(DMR.common[i] == 0) { DMR.type[i] <- ifelse(meth.1[i] > meth.0[i], "gain", "loss") }
    } #Done looping over all DMRs
    
    #Format results and return
    results <- data.frame(meth.1, meth.0, DMR.common, DMR.type)
    return(results)
} #Done


#Need a function that checks the distance of a DMR to the nearest TE
distance.TE <- function(DMRs, TEs) {
    ndmr <- nrow(DMRs) #Number of DMRs
    dist.TE <- rep(0, ndmr)
    closest.TE <- rep("", ndmr)

    #Loop over all DMRs
    for(i in 1:ndmr) {
        cur.chr <- as.character(DMRs$Chromosome[i]) #Current chromosome
        cur.start <- DMRs$start[i] #start
        cur.end <- DMRs$end[i] #end

        #In case a TE overlaps the DMR
        if(is.na(DMRs$TE.class[i]) == FALSE) {
            dist.TE[i] <- 0
            closest.TE[i] <- DMRs$TE.class[i]
        } else { #In case a TE does not overlap the current DMR
        cur.TE <- filter(TEs, seqid == cur.chr) #Only TEs in this particular chromosome
        dist.starts.2.dmr.start <- abs(cur.TE$start - cur.start)
        dist.starts.2.dmr.end <- abs(cur.TE$start - cur.end)
        dist.ends.2.dmr.start <- abs(cur.TE$end - cur.start)
        dist.ends.2.dmr.end <- abs(cur.TE$end - cur.end)

        #Need to check which distance is the closest one
        distances <- cbind(dist.starts.2.dmr.start, dist.starts.2.dmr.end, dist.ends.2.dmr.start, dist.ends.2.dmr.end)
        check.dist <- apply(distances, MARGIN = 2, min)
        index <- which(check.dist == min(check.dist))

        #Index of the TE closest to the DMR
        dist.TE[i] <- min(distances[,index]) #Distance to the nearest TE
        index.TE <- which(distances[,index] == min(distances[,index]))

        #Get the type of the closest TE
        nearest.TE <- cur.TE[index.TE,]
        myattr <- unlist(strsplit(as.character(nearest.TE$attributes), ";")) #Get info about TE type and class
        #    TE.type[i] <- unlist(strsplit(myattr[1],"'"))[2] #Extract the information
        closest.TE[i] <- unlist(strsplit(myattr[2],"'"))[2]
        } #End if-statement    
    } #Done looping over DMRs

    results <- data.frame(dist.TE, closest.TE)
    return(results)
} #Done

#This function checks the distance between TE and sampled fetures
distance.TE.feat <- function(feat, TEs) {
    nfeat <- nrow(feat) #Number of DMRs
    dist.TE <- rep(0, nfeat)
    closest.TE <- rep("", nfeat)

    #Loop over all features
    for(i in 1:nfeat) {
        cur.chr <- as.character(feat$seqid[i]) #Current chromosome
        cur.start <- feat$start[i] #start
        cur.end <- feat$end[i] #end

        #First, check for overlap of TE
        TE.test <- check.overlap(filter(TEs, seqid == cur.chr)[,4:5], cur.start, cur.end)

        #If a feature overlaps a TE
        if(any(TE.test)) {
            dist.TE[i] <- 0
            myTE <- filter(TEs, seqid == cur.chr)[TE.test,]
            myattr <- unlist(strsplit(myTE$attributes, ";")) #Get info about TE type and class
            closest.TE[i] <- unlist(strsplit(myattr[2],"'"))[2]
            #browser() 
        } else { #In case a TE does not overlap the current DMR
        cur.TE <- filter(TEs, seqid == cur.chr) #Only TEs in this particular chromosome
        dist.starts.2.dmr.start <- abs(cur.TE$start - cur.start)
        dist.starts.2.dmr.end <- abs(cur.TE$start - cur.end)
        dist.ends.2.dmr.start <- abs(cur.TE$end - cur.start)
        dist.ends.2.dmr.end <- abs(cur.TE$end - cur.end)

        #Need to check which distance is the closest one
        distances <- cbind(dist.starts.2.dmr.start, dist.starts.2.dmr.end, dist.ends.2.dmr.start, dist.ends.2.dmr.end)
        check.dist <- apply(distances, MARGIN = 2, min)
        index <- which(check.dist == min(check.dist))

        #Index of the TE closest to the DMR
        dist.TE[i] <- min(distances[,index]) #Distance to the nearest TE
        index.TE <- which(distances[,index] == min(distances[,index]))

        #Get the type of the closest TE
        nearest.TE <- cur.TE[index.TE,]
        myattr <- unlist(strsplit(as.character(nearest.TE$attributes), ";")) #Get info about TE type and class
        #    TE.type[i] <- unlist(strsplit(myattr[1],"'"))[2] #Extract the information
        closest.TE[i] <- unlist(strsplit(myattr[2],"'"))[2]
        } #End if-statement    
    } #Done looping over DMRs

    results <- data.frame(dist.TE, closest.TE)
    return(results)
} #Done

#This function loads the chipseq peaks called by Mariana
load.chipseq.peaks <- function(peakdatafolder, peaksamples, fileext) {
    files <- paste0(peakdatafolder, peaksamples, fileext) #File names to be loaded
    nsamples <- length(files) #Number of samples
    foo <- list(0)
    reslist <- rep(foo, nsamples) #List for storing the data
    #Loop over all files, load, and do some processing
    for(i in 1:nsamples) {
        current <- read.table(files[i], header = F, skip = 1)
        current <- current[,1:3] #Drop other columns, they are not needed
        colnames(current)[1:3] <- c("chr", "start", "end") #Adjust col names
        #Take only the seven chromosomes
        current <- filter(current, chr %in% c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6", "Supercontig_12.7"))
        current$sampleID <- peaksamples[i]
        reslist[[i]] <- current
    } #Done looping over all files

    allpeaks <- do.call("rbind", reslist) #Bind the peaks of all samples together
    allpeaks <- arrange(allpeaks, chr, start) #Arrange by chr and position

    return(allpeaks)
} #Done
    
#Need to function that check chipseq peak overlap among all peaks called from all samples
merge.peaks <- function(allpeaks, samples, overl.tre = 0.2) {
    allpeaks$chr <- factor(allpeaks$chr) #Make sure that chr column format is OK
    nchr <- length(levels(allpeaks$chr)) #Number of chrosomes
    chrs <- levels(allpeaks$chr)
    allpeaks$len <- allpeaks$end - allpeaks$start #Calculate peak lengths
    finalres <- NULL #initialize final result
    #First need to loop over the seven chromosomes
    for(j in 1:nchr) {
        #Take only the current chromosome
        cur.chr <- filter(allpeaks, chr == chrs[j])
        cur.chr$peak.group <- "" #Initialize result columns
        cur.chr$exp.cont <- NA
        #Then loop across peaks
        npeaks <- nrow(cur.chr) #Number of peaks in this chromosome
        for(i in 1:npeaks) {
            #Check peak overlap
            o.check <- apply(cur.chr[,2:3], MARGIN = 1, check.overlap, start.dmr = cur.chr[i,2], end.dmr = cur.chr[i,3])
            overlap.index <- (1:npeaks)[o.check] #Get the indices of peaks that overlap this peak
            #Need to deal with the case with no overlap?
            #There now should be overlap at least with itself -> No need to worry about this
            
            #Overlapping peaks
            overlap.peaks <- cur.chr[overlap.index,]
            overlap.peaks$peak.group <- paste("chr", j,":",paste(overlap.index, collapse = "/"),sep = "") #Make a string that stores overlapping peaks
            #Check if length of any peak differs more than 20% (or overl.tre) from the median peak length
            med.len <- median(overlap.peaks$len)
            overlap.peaks$exp.cont <- abs(1 - (overlap.peaks$len / med.len)) > overl.tre #Check if length is larger than overlap threshold
            #Store modified values
            cur.chr[overlap.index, 6:7] <- overlap.peaks[,6:7]
        } #Done loopint over all peaks

        #Need to store the results
        finalres <- rbind(finalres, cur.chr)
    } #Done looping over all chromosomes

    return(finalres)
}

#This function annotates the merged chipseq peaks, and makes a list of peaks
annotate.peaks <- function(merged.peaks, samples) {
    merged.peaks$peak.group <- factor(merged.peaks$peak.group) #Change to a factor
    uniq.peaks <- length(levels(merged.peaks$peak.group)) #Get the number of unique peaks
    #Make results matrix for the unique peaks
    peakresults <- data.frame( chr = rep("", uniq.peaks), start = rep(0, uniq.peaks), end = rep(0, uniq.peaks), len = rep(0, uniq.peaks), group.id = rep("", uniq.peaks), is.polymorphic = rep(F, uniq.peaks), IGV.check = rep(F, uniq.peaks), strains.IGV = rep(NA, uniq.peaks), stringsAsFactors = F )
    #First loop over all unique peaks
    #Store chr, median, start, end, and length
    for(i in 1:uniq.peaks) { #Loop over all unique peaks
        #take current group
        cur.peaks <- filter(merged.peaks, peak.group == levels(merged.peaks$peak.group)[i])
        peakresults[i,1] <- as.character(cur.peaks[1,1]) #Store the chromosome
        peakresults[i,2] <- round(median(cur.peaks[,2])) #Store median of start coordinates
        peakresults[i,3] <- round(median(cur.peaks[,3])) #Store median of end coordinates
        peakresults[i,4] <- peakresults[i,3] - peakresults[i,2] #Store peaks lengths
        peakresults[i,5] <- as.character(cur.peaks[1,6]) #Store peak group ID
        #Check whether a peak is polymorphic among the samples
        #if a peak is present in fewer samples than the total sample then it is polymorphic
        if(nrow(cur.peaks) < length(samples)) { peakresults[i,6] <- TRUE }
        #Then check does the peak need to be cheked manually in IGV
        #Store also the samples that are suspicious (this is based on peak having smaller overlap)
        if(any(cur.peaks[,7]) == T) {
            peakresults[i,7] <- TRUE
            peakresults[i,8] <- paste(cur.peaks[,4][cur.peaks[,7]], collapse = ",")
        }
        #
    } #Done looping over all unique peaks
    #
    #Sort peaks by chr and start coordinate
    peakresults <- arrange(peakresults, chr, start)
    #
    return(peakresults) #return results
} #End


#This functions makes a presence / absence data matrix of the variable peaks
peak.matrix <- function(variable, merged, samples) {
    npeaks <- nrow(variable) #Number of variable peaks
    nsamples <- length(samples) #Number of samples
    #Preparare the results matrix
    res.mat <- variable[,1:5] #Take the chr, start, end, len, and peak ids
    res.mat2 <- data.frame(matrix(rep(0,npeaks*nsamples), ncol = nsamples))
    colnames(res.mat2) <- samples
    res.mat <- cbind(res.mat, res.mat2) #Final results matrix
    ##Go over the variable peaks, then for each peak check in which individuals it occurs using the merged dataframe. Record the samples in the matrix
    for(i in 1:npeaks) {
        current <- variable[i,] #Take the current peak
        #Get all of the individuals that belong to the same peak group
        inds <- filter(merged, peak.group == current$group.id)
        samples.presence <- inds$sampleID
        res.mat[i,6:(nsamples+5)] <- as.numeric(samples %in% samples.presence) #match samples that have peaks and store that information
    } #Done looping over all peaks
    #Return the results
    return(res.mat)
}

#This function calculates the number of mitoses that separate two samples
calc.dt <- function(x1,x2) {
    #get the line numbers and transfer numbers
    #First check if focal is ancestor
    if(grepl("anc", x1) == T) {
        t1 <- 0
        l1 <- "anc"
    } else { #Otherwise do 
    split1 <- unlist(strsplit(x1, "_"))
    t1 <- as.numeric(gsub("G", "", split1[1]))
    l1 <- split1[2]
    }

    #Second sample
    #First check if focal is ancestor
    if(grepl("anc", x2) == T) {
        t2 <- 0
        l2 <- "anc"
    } else { #Otherwise do 
    split2 <- unlist(strsplit(x2, "_"))
    t2 <- as.numeric(gsub("G", "", split2[1]))
    l2 <- split2[2]
    }

    #Then calculate divergence time
    #Check are samples from the same line
    if(l1 == l2) {
        dt <- max(c(t1,t2)) - min(c(t1,t2)) + 2 #2 comes from DNA extraction for both samples
    } else { #If samples are from different lines, divergence calculated as a sum
        dt <- t1 + t2 + 2 #2 comes from DNA extraction for both samples
    }

    return(dt)
}

#This function calculates divergence between samples in H3K9 using presence and absence of peaks 
calc.h3k9.div <- function(peak.matrix, sample.comb) {
    div <- rep(0, nrow(sample.comb)) #initialize divergence
    npeaks <- nrow(peak.matrix) #number of peaks
    ncomb <- nrow(sample.comb) #number of combinations
    for(i in 1:ncomb) { #Loop over all sample combinations
        current <- peak.matrix[,sample.comb[i,]] #Take peaks of current combination
        div[i] <- sum(current[,1] != current[,2]) / npeaks #Calculate divergence: how many are different / total number of peaks
    }
    return(div)
}

#Calculates divergence, needs sample combinations and normalized H3K9 read counts as input
h3k9.count.divergence <- function(samples, counts) {
    ncomb <- nrow(samples)
    div.results <- rep(0, ncomb) #Initialize results vector
    #Loop over all sample combinations
    for(i in 1:ncomb) {
        s1 <- samples[i,1] #Take the samples
        s2 <- samples[i,2]

        #Get the sample counts
        s1.count <- counts[,s1]
        s2.count <- counts[,s2]

        div.results[i] <- sum((abs(s1.count - s2.count)))/length(s1.count) #Calculate divergence in read counts across all windows
    } #Done looping across all combinations
    return(div.results) #Return divergence
} #End


##Need to make a function that check DMR changes in all intervals
intervals.DMR.check <- function(aineisto, DMRs.all.nanopore, ind.cols) {
    #ind.cols <- 5:42
    inds <- DMRs.all.nanopore[,ind.cols]
    ind.names <- colnames(inds)

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
            cur.ind <- grepl(paste0(current.gen,"_"), ind.names) & grepl(current.line, ind.names)
            cur.ind <- which(cur.ind)
            current <- inds[,cur.ind] #CUrrent focal sample

            changes <- current != DMRs.all.nanopore$DMR.common
            res.mat[i,1] <- sum(changes) #Total number of changes in interval
            types <- DMRs.all.nanopore$DMR.type[changes]
            res.mat[i,2] <- sum(types == "gain") #Number of DMR gains
            res.mat[i,3] <- sum(types == "loss") #Number of DMR losses
        } #Done counting DMR changes

        #If not in first interval, need to compare with previous
        if(current.gen != "G1") {
            #Get the current ind
            cur.ind <- grepl(paste0(current.gen,"_"), ind.names) & grepl(current.line, ind.names)
            cur.ind <- which(cur.ind)
            current <- inds[,cur.ind] #Current focal sample
            #Check which was the previous transfer
            if(current.gen == "G5") { previous.gen <- "G1" }
            if(current.gen == "G7") { previous.gen <- "G5" }
            if(current.gen == "G8") { previous.gen <- "G7" }
            if(current.gen == "G10") { previous.gen <- "G8" }
            if(current.gen == "G15") { previous.gen <- "G10" }
            prev.ind <- grepl(paste0(previous.gen,"_"), ind.names) & grepl(current.line, ind.names)
            prev.ind <- which(prev.ind)
            previous <- inds[,prev.ind] #Sample from previous transfer

            changes <- current != previous
            res.mat[i,1] <- sum(changes) #Total number of DMR changes
            cur.cha <- current[changes]
            prev.cha <- previous[changes]
            types <- ifelse(cur.cha == 0 & prev.cha == 1, "loss", "gain")
            res.mat[i,2] <- sum(types == "gain") #Number of DMR gains
            res.mat[i,3] <- sum(types == "loss") #Number of DMR losses
        } #Done counting DMR changes
    } #Done looping over all intervals

    #Return results
    return(cbind(aineisto, res.mat))
} #Done
