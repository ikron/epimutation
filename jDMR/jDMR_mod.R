### Modified version of the jDMR package

### * CfromFASTAv4.R
#' Extract cytosines from FASTA
#'
#' This function slides along the genome and extracts all Cytosines
#' @param fasta
#' @param chr
#' @param out.dir
#' @param write.output
#' @import Biostrings
#' @import R.utils
#' @import data.table
#' @export
#'

### CfromFASTAv1:
### - Option to return tri-nucleotide contexts if needed
### - Bioconstrings for finding positions of contexts

 CfromFASTAv4<-function(fasta, chr, out.dir, write.output)
 {

   # Defining function for cleaning objects in the global environment from within the function
   CleanEnvir <- function(pattern)
   {
     rm(list = ls(envir=globalenv())[
       grep(pattern, ls(envir=globalenv()))], envir = globalenv())
   }

      cat("- Processing chr:", chr, "\n")
      cat("-------------------------", "\n")

      cat("- Converting DNA .....", "\n")
      fasta.minus<-Biostrings::complement(fasta)
      fasta.plus<-Biostrings::DNAString(paste(fasta, collapse=""))
      fasta.minus<-Biostrings::DNAString(paste(fasta.minus, collapse=""))

      # Removing the .fasta file from global environment
      CleanEnvir(pattern="fasta")

      # Listing the cytosine patterns by context
      plus.pattern<-list("CG",  c("CAG", "CTG", "CCG"), c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC"))
      names(plus.pattern)<-c("CG", "CHG", "CHH")

      minus.pattern<-list("GC", c("GTC","GAC", "GCC"), c("AAC", "TAC", "CAC", "ATC", "TTC", "CTC", "ACC", "TCC","CCC"))
      names(minus.pattern)<-c("CG", "CHG", "CHH")

      cytosine.patterns<-list(plus.pattern, minus.pattern)
      names(cytosine.patterns)<-c("+", "-")

        collect.positions<-list()
        counter=0

        for (s in 1:length(cytosine.patterns))
        {
            strand.temp<-cytosine.patterns[[s]]

              for (c in 1:length(strand.temp))
              {
                cat("- Scanning", names(strand.temp)[c], names(cytosine.patterns)[s], "strand", ".....", "\n")
                context.temp<-strand.temp[[c]]

                for (p in 1:length(context.temp))
                {
                  counter<-counter+1

                   if (names(cytosine.patterns)[[s]] == "+")
                   {
                     match.out<-Biostrings::matchPattern(context.temp[p], subject=fasta.plus)
                     match.out<-as.data.frame(slot(match.out, name="ranges"))[,1]
                     match.out<-data.table(rep(chr, length(match.out)), match.out, rep(names(cytosine.patterns)[[s]], length(match.out)),
                                           rep(names(strand.temp)[c], length(match.out)))
                     colnames(match.out)<-c("chr", "pos", "strand", "context")
                     collect.positions[[counter]]<-match.out
                   }
                   if (names(cytosine.patterns)[[s]] == "-")
                   {
                     match.out<-Biostrings::matchPattern(context.temp[p], subject=fasta.minus)
                     match.out<-as.data.frame(slot(match.out, name="ranges"))[,2]
                     match.out<-data.table(rep(chr, length(match.out)), match.out, rep(names(cytosine.patterns)[[s]], length(match.out)),
                                           rep(names(strand.temp)[c], length(match.out)))
                     colnames(match.out)<-c("chr", "pos", "strand", "context")
                     collect.positions[[counter]]<-match.out
                   }
                }

              }

        }

        cat("- Combining contexts .....", "\n")
        C.all<-data.table(do.call("rbind", collect.positions))

        # Reading out the data
        if (write.output == TRUE)
        {
          cat("- Writing out file .....", "\n")
          fwrite(C.all, paste(out.dir, "/cytosine_positions_chr", chr, ".csv", sep=""), row.names = FALSE)
        }

        #return(C.all)

        # Cleaning up
        rm(fasta, fasta.plus, fasta.minus, collect.positions, match.out, C.all)


 }

### * MethimputeReg.R

#'
#' @param distcor
#' @param skip
#' @param plot.parameters
#' @import ggplot2
#' @import stats
#' @import minpack.lm
#' @export
#'
modified.estimateTransDist <- function(distcor, skip=2, plot.parameters=TRUE) {

  ## Context correlation fits and plots
  contexts <- dimnames(distcor$data)[[1]]
  cor.array <- distcor$data
  maxweights <- numeric()
  params.list <- list()
  miny <- min(cor.array, na.rm = TRUE)
  dfs <- list()
  for (c1 in 1:length(contexts)) {
    for (c2 in 1:length(contexts)) {
      context.transition <- paste0(contexts[c1], '-', contexts[c2])
      if (distcor$separate.contexts) {
        if (c1 != c2) {
          next
        }
      }
      if (c1 <= c2) {
        df <- data.frame(distance = as.numeric(dimnames(cor.array)[[3]]),
          correlation = cor.array[c1,c2,,'correlation'],
          weight = cor.array[c1,c2,,'weight'],
          from = contexts[c1], to = contexts[c2])

        ## Fit
        y <- df$correlation[(skip+1):nrow(df)]
        x <- df$distance[(skip+1):nrow(df)]
        weight <- df$weight[(skip+1):nrow(df)]
        startvalues <- list(a0 = stats::na.omit(y)[1], D = 50)
        p <- NULL

        if (is.null(p)) {
          startvalues <- list(a0 = stats::na.omit(y)[1])
          p <- tryCatch({
            fit <- minpack.lm::nlsLM(y ~ a0 * exp(-x/Inf), start=startvalues, weights=weight)
            s <- summary(fit)
            c <- stats::coefficients(s)
            params <- c[1:length(startvalues)]
            names(params) <- names(startvalues)
            params <- as.list(params)
            params$D <- Inf
            params
          }, error = function(e) {
            startvalues$D <- Inf
            startvalues
          })
        }

        ## Check if we have negative D
        if (p$D <= 0) {
          p$D <- Inf
        }
        params.list[[context.transition]] <- p

        ## Plot
        df$correlation.fit <- p$a0 * exp(-df$distance/p$D)
        df$logweight <- log(df$weight+1)
        dfs[[context.transition]] <- df
        maxweights[context.transition] <- max(df$logweight, na.rm = TRUE)
      }
    }
  }
  maxweight <- max(maxweights, na.rm = TRUE)

  ## Plot correlation
  df <- do.call(rbind, dfs)
  df$a0 <- round(sapply(params.list[paste0(df$from, '-', df$to)], '[[', 'a0'), 2)
  df$D <- round(sapply(params.list[paste0(df$from, '-', df$to)], '[[', 'D'), 0)
  df$params <- paste0("a0 = ", df$a0, ", D = ", df$D)
  ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(x='distance', y='correlation', alpha='logweight'))
  ggplt <- ggplt + geom_line(aes_string(x='distance', y='correlation.fit'), col='blue')
  if (plot.parameters) {
    ggplt <- ggplt + geom_text(aes_string(label='params'), x=max(df$distance, na.rm = TRUE), y=max(df$correlation, na.rm = TRUE), vjust=1, hjust=1)
  }
  ggplt <- ggplt + xlab('distance in [bp]')
  ggplt <- ggplt + facet_grid(from ~ to)
  if (miny < 0) {
    ggplt <- ggplt + geom_hline(aes_string('yintercept'=0), linetype=2, alpha=0.5)
  }

  transDist <- sapply(params.list, '[[', 'D')
  return(list(transDist=transDist, plot=ggplt))
}

#--------------------------------------------------------------------------
#' @param model
#' @param out.dir
#' @param context
#' @param name
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom stringr str_replace_all
#' @export
#'
modifiedExportMethylome <- function(model, out.dir, context, name) {
    #data <- model$data
    data <- model
    final_dataset <- as(data, 'data.frame')
    final_dataset <- final_dataset[,c('seqnames','start','end','strand',
      'context','counts.methylated','counts.total',
      'posteriorMax','status','rc.meth.lvl')]

    # dropping columns
    drops <- c('width','strand','clusterlen','counts.methylated',
      'counts.total', 'distance', 'transitionContext', 'posteriorMeth','posteriorUnmeth')
    final_dataset <- final_dataset[ , !(names(final_dataset) %in% drops)]
    #------------------------------------------------------------------
    # convert full string into M/U/I
    final_dataset <- statusStringCheck(final_dataset)
    #------------------------------------------------------------------
    # take 4 digit of decimal value posteriorMax column
    final_dataset$posteriorMax <-floorDec(as.numeric(as.character(final_dataset$posteriorMax)),5)
    final_dataset$rc.meth.lvl <- floorDec(as.numeric(as.character(final_dataset$rc.meth.lvl)),5)
    final_dataset$seqnames <- as.character(final_dataset$seqnames)

    saveFile <- paste0(out.dir, "/", name, "_", context, ".txt")
    fwrite(final_dataset, file = saveFile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    return (final_dataset)
}

#--------------------------------------------------------------------------
#' @param df
#' @param refRegion
#' @param context
#' @param mincov
#' @param nCytosines
#' @import dplyr
#' @importFrom data.table fread
#' @import stats
#' @import GenomicRanges
#' @import S4Vectors
#' @export
#'
makeRegionsImpute <- function(df, context, refRegion, mincov, nCytosines) {

  #regions file
  tmp_reg <- dget(refRegion)
  data <- as.data.frame(tmp_reg$reg.obs)
  #colnames(data)[which(names(data) == "cluster.size")] <- "cluster.length"

  #reference methimpute file
  ref_data <- fread(df, skip = 1, select = c("V1","V2","V3","V4","V5","V6"))

  #remove Mt and chloroplast coordinates. Following is for Arabidopsis only
  ref_data <- ref_data %>% dplyr::filter(ref_data$V1 != "M" & ref_data$V1 != "C")

  #filtering by coverage
  if (mincov>0){
    cat(paste0("Filtering for coverage: ", mincov,"\n"), sep = "")
    ref_data <- ref_data[which(ref_data$V6 >= mincov),]
  }

  ref_data <- ref_data[which(ref_data$V4==context),]

  data_gr <- GRanges(seqnames=data$chr,
                     ranges=IRanges(start=data$start, end=data$end),
                     clusterlen=data$cluster.length,
                     context=as.factor(context))

  ref_gr <- GRanges(seqnames=ref_data$V1,
                    ranges=IRanges(start=ref_data$V2, width=1),
                    context=as.factor(context),
                    methylated=ref_data$V5,
                    total=ref_data$V6)

  data_gr$cytosineCount <- GenomicRanges::countOverlaps(data_gr, ref_gr)

  #filtering by number of cytosines
  if (nCytosines>0){
    cat(paste0("Minimum cytosines covered per region: ", nCytosines,"\n"), sep = "")
    data_gr <- data_gr[which(data_gr$cytosineCount >= nCytosines),]
  }

  counts <- array(NA, dim=c(length(data_gr), 2), dimnames=list(NULL, c("methylated", "total")))

  overlaps <- IRanges::findOverlaps(ref_gr, data_gr)

  overlaps.hits <- ref_gr[S4Vectors::queryHits(overlaps)]
  if (NROW(overlaps.hits) != 0){

    methylated <- stats::aggregate(overlaps.hits$methylated, list(S4Vectors::subjectHits(overlaps)), FUN=sum)
    total <- stats::aggregate(overlaps.hits$total, list(S4Vectors::subjectHits(overlaps)), FUN=sum)

    if (NROW(methylated) != NROW(counts) ){
      missingr <- which(!rownames(data.frame(data_gr)) %in% methylated$Group.1)

      for(item in seq_len(NROW((missingr)))){
        methylated <- rbind (c(missingr[item],0), methylated)
        methylated <- methylated[order(methylated$Group.1),]
        total <- rbind (c(missingr[item],0), total)
        total <- total[order(total$Group.1),]
      }
    }
    counts[,"methylated"] <- methylated$x
    counts[,"total"] <- total$x
    data_gr$counts <- counts
  }
  rm(ref_data, ref_gr, overlaps, overlaps.hits)
  return(data_gr)
}

#--------------------------------------------------------------------------
#' @param df
#' @param context Cytosine context
#' @param fit.plot
#' @param fit.name
#' @param refRegion
#' @param include.intermediate
#' @param probability
#' @param out.dir
#' @param name
#' @param mincov
#' @param nCytosines
#' @importFrom methimpute callMethylation
#' @importFrom methimpute distanceCorrelation
#' @export

makeMethimpute <- function(df, context, fit.plot, fit.name, refRegion,
                         include.intermediate, probability, out.dir, name, mincov, nCytosines){
  methylome.data <- makeRegionsImpute(df, context, refRegion, mincov, nCytosines)
  if (!is.null(methylome.data$counts)) {
    quant.cutoff <- as.numeric(quantile(methylome.data$counts[,"total"], probs = c(0.96), na.rm=TRUE))
    distcor <- distanceCorrelation(methylome.data, distances=0:100)
    fit <- modified.estimateTransDist(distcor)

    if (fit.plot==TRUE){
      print(paste0("Generating fit plot...", name))
      pdf(paste0(out.dir, "/", fit.name, "-fit.pdf", sep = ""))
      print(fit)
      dev.off()
    }

    model <- callMethylation(data = methylome.data,
                             transDist = fit$transDist,
                             count.cutoff = quant.cutoff,
                             max.time = Inf,
                             max.iter = Inf,
                             include.intermediate = include.intermediate,
                             update = probability)

    methFile <- modifiedExportMethylome(model=model$data, out.dir=out.dir, context=context, name=name)
    rm(model)
  }
  rm(methylome.data)
}

### * MethimputeRegTobedGraph.R



MethimputeRegTobedGraph.rcmethlvl <- function(regfile, out.dir) {
  cat(paste0("Reading Methimpute file: ", regfile, "\n"))
  fname <- fread(regfile, skip = 1, select = c("V1","V2","V3","V7"))
  name <- gsub(pattern = "\\.txt$", "", basename(regfile))
  fname$V2 <- format(fname$V2, scientific = FALSE) 
  fname$V3 <- format(fname$V3, scientific = FALSE)
  print(paste0("Writing to bedGraph format....",name))
  fwrite(x= fname, file=paste0(out.dir, "/", name,"-rcmethlvl.bedGraph"), sep="\t", 
         quote=FALSE, row.names=FALSE, col.names=FALSE)
}

MethimputeRegTobedGraph.stateCalls <- function(regfile, out.dir) {
  cat(paste0("Reading Methimpute file: ", regfile, "\n"))
  fname <- fread(regfile, skip = 1, select = c("V1","V2","V3","V6"))
  fname$V6<- ifelse(fname$V6 == "U", yes=0, no=1)
  name <- gsub(pattern = "\\.txt$", "", basename(regfile))
  fname$V2 <- format(fname$V2, scientific = FALSE) 
  fname$V3 <- format(fname$V3, scientific = FALSE)
  print(paste0("Writing to bedGraph format....",name))
  fwrite(x= fname, file=paste0(out.dir, "/", name,"-stateCalls.bedGraph"), sep="\t", 
         quote=FALSE, row.names=FALSE, col.names=FALSE)
}

### * annotateDMRs.R

#'
#' @param gff
#' @importFrom  rtracklayer import.gff3
#' @export
#' @return merge all supplied gff3 annotations into one
#merge the input gff3 files into one
gff3.in <- function(gff){
  input.gff <- lapply(gff, function(x){
    import.gff3(x, colnames=c("type", "ID"))
  })
  merged.gff <- do.call(c, input.gff)
  return(merged.gff)
}

# Check Annotation levels here. Supply annotation terms for e.g genes, TEs
#levels(elementMetadata(merged.gff)[,"type"])
#available annotations
#"chromosome","gene","mRNA","five_prime_UTR","exon","CDS",
#"three_prime_UTR","ncRNA_gene","lnc_RNA","miRNA","tRNA","ncRNA",
#"snoRNA","snRNA","rRNA","TE","promoters"

#' @param gff
#' @param annotation
#' @param grangesObj
#' @param name
#' @param out.dir output directory
#' @import GenomicRanges
#' @importFrom  rtracklayer export.gff
#' @export
#' @return export output files in gff3 format
#output annotated gff3 files
gff3.out <- function(annotation, grangesObj, gff, name, out.dir) {
  getgff3 <- lapply(annotation, function(x){
    idx <- which(elementMetadata(gff)[,"type"] == x)
    gff <- gff[idx,]
    hits <- findOverlaps(grangesObj, gff, ignore.strand=FALSE)
    gr.matched <- grangesObj[queryHits(hits)]
    if (NROW(gr.matched) != 0){
      mcols(gr.matched) <- cbind.data.frame(mcols(gr.matched), mcols(gff[subjectHits(hits)]))
      values(gr.matched) <- cbind(values(gr.matched), region="DMR")
      names(elementMetadata(gr.matched))[names(elementMetadata(gr.matched)) == "type"] <- "annotation"
    }
    return(gr.matched)
  })
  export.gff(do.call(c, getgff3), paste0(out.dir,"/", name, "_annotation.gff3"), version="3")
}

#' @inheritParams gff3.in
#' @param gff
#' @param annotation
#' @param grangesObj
#' @param name
#' @param out.dir output directory
#' @param getAnno
#' @param mygff
#' @param mygr
#' @param annotation
#' @param gff.files
#' @param gff3.out
#' @param input.dir
#' @import  GenomicRanges
#' @export
#' @return annotated list
#extract annotated regions
annotate <- function(getAnno, mygff, mygr){
  lapply(getAnno, function(x){
    idx <- which(elementMetadata(mygff)[,"type"] == x)
    mygff <- mygff[idx,]
    hits <- findOverlaps(mygr, mygff, ignore.strand=FALSE)
    myranges <- subsetByOverlaps(mygr, mygff)

    mcols(myranges)$id <- CharacterList(split(mygff$ID[subjectHits(hits)], queryHits(hits)))
    mcols(myranges)$type <- CharacterList(split(mygff$type[subjectHits(hits)], queryHits(hits)))
    mcols(myranges)$chr <- CharacterList(split(seqnames(mygff)[subjectHits(hits)], queryHits(hits)))
    mcols(myranges)$coord <- CharacterList(split(ranges(mygff)[subjectHits(hits)], queryHits(hits)))
    mcols(myranges)$anno.str <- CharacterList(split(strand(mygff)[subjectHits(hits)], queryHits(hits)))
    df <- as(myranges, "data.frame")
    df <- df %>% mutate(id = strsplit(as.character(id), ","),
                        type = strsplit(as.character(type), ","),
                        chr = strsplit(as.character(chr), ","),
                        coord = strsplit(as.character(coord), ","),
                        anno.str = strsplit(as.character(anno.str), ",")) %>%
      tidyr::unnest(c(id, type, chr, coord, anno.str))
    cleandf <- data.frame(lapply(df, function(k) gsub ("[\\c]|[()]|\"|^ .", "", k)))
    cleandf$Annotation.coord <- apply(cleandf[,c("chr","coord","anno.str")], 1, paste, collapse=":")
    cleandf <- subset(cleandf, select = -c(chr, coord, anno.str))
    return(cleandf)
  })
}

#' Annotate DMRs
#'
#' This function takes gff3 files as input and outputs annotated DMRs in text and gff3 format. Additionally, a DMR count table is generated.
#'
#' @inheritParams gff3.in
#' @inheritParams gff3.out
#' @inheritParams annotate
#' @param out.dir output directory
#' @param annotation annotation terms used to annotate DMRs
#' @param gff.files multiple gff3 files can be supplied as a vector
#' @param gff3.out a logical specifying whether output annotated files in gff3 format
#' @param input.dir input directory containing filtered DMR matrix/matrices. Ideally any file containing 3 columns i.e (chr, start, stop) can be supplied.
#' @import GenomicRanges
#' @importFrom  data.table fread
#' @importFrom  data.table fwrite
#' @importFrom  dplyr group_by
#' @importFrom  dplyr summarize
#' @export
#' @return output files containing annotated DMRs and DMR counts table.

annotateDMRs <- function(annotation, gff.files, gff3.out, input.dir, out.dir) {
  anno.list <- list(); final.df <- list()
  file.list <- list.files(input.dir, pattern="*.txt", full.names = TRUE)
  for (i in 1:length(file.list)){

    cat(paste0("Running file ", file.list[i], "\n"), sep = "")
    tmp.name <- gsub("\\.txt$", "", basename(file.list[i]))
    file <- fread(file.list[i], skip=1, select=c(1,2,3))
    gr <- GRanges(seqnames=file$V1, ranges=IRanges(start=file$V2, end=file$V3))
    gff=gff3.in(gff.files)

    if (gff3.out==TRUE) {
      gff3.out(annotation=annotation,
               gff=gff,
               grangesObj=gr,
               name=tmp.name,
               out.dir=out.dir)
    }

    d <- rbindlist(annotate(getAnno=annotation, mygff=gff, mygr=gr))
    if (nrow(d) != 0){
      out <- d %>%
        dplyr::group_by(seqnames, start, end) %>%
        dplyr::summarize(
          type = paste(type, collapse=","),
          id = paste(id, collapse=","),
          Annotation.coord = paste(Annotation.coord, collapse=",")
          #, .groups = 'drop'
          ) %>% as.data.frame()

      if (NROW(out)!=0){
        for (k1 in 1:NROW(out)){
          out$unique.anno.type[k1] <- paste0(sapply(strsplit(out$type[k1],","),unique), collapse=",")
        }
      }
      # count the DMR overlaps; the output can be used to make a barplot or pie-chart
      # unique annotation overlaps

      out.1 <- out[which(sapply(strsplit(out$type,','), uniqueN)==1),]
      out.2 <- out[which(sapply(strsplit(out$type,','), uniqueN)>=2),]

      #counting unique annotations
      for (k2 in seq_along(annotation)){
        anno.list[[k2]] <- NROW(out.1[grep(annotation[k2], out.1$type),])
        names(anno.list)[[k2]] <- annotation[k2]
      }
      #also counting multiple overlaps
      anno.list[length(annotation)+1] <- list(NROW(out.2))
      names(anno.list)[[length(annotation)+1]] <- "multiple.overlaps"
      df.1 <- do.call(cbind, anno.list)

      df <- cbind(sample=tmp.name, total.DMRs=NROW(gr), df.1)
      final.df <- rbind(final.df, data.frame(df))

      out <- out[,-c(4,6)]
      fwrite(x=out, file=paste0(out.dir, "/", tmp.name, "_annotation.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    }
    fwrite(x=final.df, file=paste0(out.dir, "/DMR-counts.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  }
}

### * filterDMRmatrix.R

#'
#' @param status.collect
#' @param rc.methlevel.collect
#' @param replicate.consensus
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom dplyr semi_join
#' @import R.utils
#' @export
#'

filterReplicateConsensus <- function(status.collect, rc.methlevel.collect, replicate.consensus){

  if (!is.null(status.collect$epiMAF)){
    status.collect <- status.collect[,-c("epiMAF")]
  }
  # deducing replicates info
  mycol <- names(status.collect)[4:NCOL(status.collect)]
  sampleinfo <- data.frame(do.call(rbind, strsplit(as.character(mycol),"_")))
  colnames(sampleinfo) <- c("sample", "replicate")

  dt <- data.frame()
  pb1 <- txtProgressBar(min = 1, max = NROW(status.collect), char = "=", style = 3, file = "")

  q <- lapply(1:NROW(status.collect), function(x){
    mypattern <- unlist(status.collect[x, 4:NCOL(status.collect)])
    df.bind <- cbind(sampleinfo, mypattern)
    for (m in unique(df.bind$sample)){
      rval <- round(replicate.consensus * length(df.bind$mypattern[df.bind$sample==m]))
      pattern.vals <- df.bind$mypattern[df.bind$sample==m]
      tt <- table(pattern.vals)
      if (max(tt) >= rval){
        df.bind$count[df.bind$sample==m] <- 0
      } else {
        df.bind$count[df.bind$sample==m] <- 1
      }
    }
    Sys.sleep(1/NROW(status.collect))
    setTxtProgressBar(pb1, x)
    close(pb1)
    #print(df.bind)
    if (sum(df.bind$count)==0) {
      dt <- rbind(dt, status.collect[x,])
    }
  })
  out <- q[!sapply(q,is.null)]
  #status.collect <- q[-which(sapply(q, is.null))]
  df.status.collect <- rbindlist(out)
  if (NROW(df.status.collect) !=0){
  df.rc.methlevel.collect <- rc.methlevel.collect %>% dplyr::semi_join(df.status.collect, by=c("seqnames","start","end"))
  return(list(df.status.collect, df.rc.methlevel.collect))
  } else {
    message("\nEmpty dataframe. Nothing to write!")
    return(NULL)
  }
}

#' @param mat1
#' @param mat2
#' @param epiMAF
#' @import R.utils
#' @importFrom dplyr semi_join
#' @export
#'

filterEpiMAF <- function(mat1, mat2, epiMAF){

  pb2 <- txtProgressBar(min = 1, max = NROW(mat1), char = "=", style = 3, file = "")
  mat1$epiMAF <- 0

  for (i1 in 1:NROW(mat1)){
    mypattern <- unlist(mat1[i1, 4:(NCOL(mat1)-1)])
    mycount <- table(mypattern)
    epiMAF.out <- min(mycount)/length(mypattern)
    mat1$epiMAF[i1] <- floorDec(as.numeric(as.character(epiMAF.out)),5)

    Sys.sleep(1/NROW(mat1))
    setTxtProgressBar(pb2, i1)
  }
  close(pb2)

  df.status.collect <- mat1[which(mat1$epiMAF < epiMAF),]
  if (NROW(df.status.collect) !=0){
    df.rc.methlevel.collect <- mat2 %>% dplyr::semi_join(df.status.collect, by=c("seqnames","start","end"))
    return(list(df.status.collect, df.rc.methlevel.collect))
  } else {
    message("\nEmpty dataframe. Nothing to write!")
    return(NULL)
  }
}

#'
#' @param rcmethlvl
#' @param statecalls
#' @param gap
#' @import GenomicRanges
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @import R.utils
#' @export
#'
merge.bins <- function(rcmethlvl, statecalls, gap){
  mylist <- list()
  result <- list()

  matrix1 <- as.data.frame(statecalls)
  matrix2 <- as.data.frame(rcmethlvl)

  # extract unique pattern
  extract.pattern <- unique(matrix1[,(4:NCOL(matrix1))])
  extract.pattern$pattern <- apply(extract.pattern[,c(1:NCOL(extract.pattern))], 1, paste, collapse="")

  # the state call matrix
  gr1 <- GRanges(seqnames=matrix1$seqnames, ranges=IRanges(start=matrix1$start, end=matrix1$end))
  values(gr1) <- cbind(values(gr1), pattern=apply(matrix1[,c(4:NCOL(matrix1))], 1, paste, collapse=""))

  # the rcmethlvl matrix, also add the state-call pattern
  gr2 <- GRanges(seqnames=matrix2$seqnames, ranges=IRanges(start=matrix2$start, end=matrix2$end))
  values(gr2) <- cbind(values(gr2), values(gr1), DataFrame(matrix2[,c(4:NCOL(matrix2))]))

  # this is for the state-calls: collapse overlapping bins if pattern is same
  grl_reduce <- unlist(GenomicRanges::reduce(split(gr1, gr1$pattern)))
  result <- sort(grl_reduce)
  result$pattern <- names(result)
  result <- data.frame(result)
  final.status.collect <- result %>% dplyr::left_join(extract.pattern, by=c("pattern"))

  # this is for the rcmethlvl: collapse bins and take average of the bins
  mycols <- colnames(matrix1)[4:NCOL(matrix1)]
  fn = function(u){
    out = GenomicRanges::reduce(u)
    for (x in 1:length(mycols)){
      eval(parse(text=paste0("out$", mycols[x], " = mean(u$", mycols[x], ")")))
    }
    return(out)
  }

  message("\nNow, Merging overlapping and consecutive bins...\n")
  # split rcmthlvl matrix based on pattterns
  grl <- split(gr2, gr2$pattern)

  pb3 <- txtProgressBar(min = 1, max = length(grl), char = "=", style = 3, file = "")

  for (x in 1:length(grl)){
    a <- data.frame(grl[[x]])
    a1 <- a %>% arrange(pattern, start) %>% group_by(pattern) %>% mutate(indx = cumsum(start > lag(end, default = start[1]) + gap))
    a1.gr <- makeGRangesFromDataFrame(a1, keep.extra.columns=TRUE)
    df <- lapply(split(a1.gr, a1.gr$indx), fn)
    mylist[[x]] <- df

    Sys.sleep(0.05)
    setTxtProgressBar(pb3, x)
  }
  close(pb3)

  f.df <- unlist(mylist)
  f.df <- do.call(rbind, lapply(f.df, data.frame))
  f.df <- f.df[order(f.df[,1], f.df[,2]),]
  f.df[,(6:NCOL(f.df))] <- lapply(f.df[,(6:NCOL(f.df))], function(xy){ floorDec(xy,5) })

  final.status.collect <- subset(final.status.collect, select = -c(strand, pattern))
  final.rcmethlvl.collect <- subset(f.df, select = -c(strand))
  return(list(final.status.collect, final.rcmethlvl.collect))
}

export.out <- function(out.rcmethlvl, out.statecalls, context, out.name1, out.name2, data.out){
  fwrite(x=out.statecalls, file=paste0(data.out, "/", context, "_", out.name1, ".txt"),
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  fwrite(x=out.rcmethlvl, file=paste0(data.out, "/", context, "_", out.name2, ".txt"),
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}

#' filter DMR matrix
#'
#' Filters non-polymorphic patterns by default.
#' @param epiMAF.cutoff Filter for Minor Epi-Allele frequency. By default this option is set to NULL. Applicable for population level data.
#' @param replicate.consensus Percentage of concordance in methylation states among samples with multiple replicates. By default this option is set to NULL. Useful for Control/treatment data
#' @param gridDMR Logical specifying if grid DMR approach was used for calling DMRs.
#' @param data.dir Directory containing DMR matrix files. Looks for files with suffix, _StateCalls.txt and _rcMethlvl.txt
#' @importFrom data.table fread
#' @import R.utils
#' @export
#'
filterDMRmatrix <- function(epiMAF.cutoff=NULL, replicate.consensus=NULL, gridDMR, data.dir) {

  list.status <- list.files(data.dir, pattern="_StateCalls.txt", full.names=TRUE)
  if (length(list.status) != 0){
    for (i in seq_along(list.status)){
      context <- gsub("_StateCalls.txt", "", basename(list.status[i]))
      cat("\n")
      cat(paste0("Running DMR matrix for ", context, "\n"), sep = "")
      cat("\n")

      #----------------------------------------------
      # Removing non-polymorphic/unchanged patterns
      #----------------------------------------------
      status.collect <- fread(list.status[i], header=T)
      rc.methlevel.collect <- fread(paste0(data.dir, "/", context, "_rcMethlvl.txt"), header=T)

      cat(paste0("Removing non-polymorphic patterns...\n"))

      index <- which(rowSums(status.collect[,4:NCOL(status.collect)]) != 0 &
                       rowSums(status.collect[,4:NCOL(status.collect)]) != NCOL(status.collect)-3)

      status.collect <- status.collect[index,]
      rc.methlevel.collect <- rc.methlevel.collect[index,]

      if (is.null(epiMAF.cutoff) && is.null(replicate.consensus)) {
        message("Both, epiMAF and replicate consensus set to NULL")
        out1=status.collect
        out2=rc.methlevel.collect
        export.out(out.statecalls=out1,
                   out.rcmethlvl=out2,
                   context=context,
                   out.name1="StateCalls-filtered",
                   out.name2="rcMethlvl-filtered",
                   data.out=data.dir)
      }
      #----------------------------------------------
      # Optional. Filtering out regions with epiMAF < Minor Allele Frequency
      #----------------------------------------------
      if (!is.null(epiMAF.cutoff)) {
        cat(paste0("Filtering for epiMAF: ", epiMAF.cutoff, "\n"))
        cat("\n")
        mydf <- filterEpiMAF(mat1=status.collect, mat2=rc.methlevel.collect, epiMAF=epiMAF.cutoff)
        if (!is.null(mydf)){
          # For Population data remove the epiMAF column
          out1=mydf[[1]][,-c("epiMAF")]
          out2=mydf[[2]]
          export.out(out.statecalls=mydf[[1]],
                     out.rcmethlvl=mydf[[2]],
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        }
      }

      #----------------------------------------------
      # Optional. Retaining samples based on replicate.consensus
      #----------------------------------------------
      if (!is.null(replicate.consensus)) {
        cat(paste0("Filtering for replicate consensus...\n"))
        cat("\n")
        mydf <- filterReplicateConsensus(status.collect, rc.methlevel.collect, replicate.consensus)
        if (!is.null(mydf)){
          out1=mydf[[1]]
          out2=mydf[[2]]
          export.out(out.statecalls=out1,
                     out.rcmethlvl=out2,
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        }
      }

      #----------------------------------------------
      # Merging bins
      #----------------------------------------------
      if (gridDMR==TRUE) {
        if (exists("out1") && exists("out2")){
          out <- merge.bins(statecalls=out1, rcmethlvl=out2, gap=1)
          export.out(out.statecalls=out[[1]],
                     out.rcmethlvl=out[[2]],
                     context=context,
                     out.name1="StateCalls-filtered-merged",
                     out.name2="rcMethlvl-filtered-merged",
                     data.out=data.dir)
        }
      } else {
        message("\ngrid DMR set to FALSE")
        }
    }
  } else {
    message("\nDMR matrix files do not exist!")
  }
}

### * globFun.R


# set number of cores based on cpu, number of jobs
numCore <- function(NumFiles){
  no_cores <- detectCores()
  if (no_cores > NumFiles){
    no_cores <- NumFiles
  }else{
    no_cores <- no_cores - 1
  }
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  getDoParWorkers()
  return(cl)
}

# take 5 digit of decimal value and cut the numbers
floorDec <- function(valParm ,x){
  y <- function(x, level=1) round(x - 5*10^(-level-1), level)
  res <-y(as.numeric(valParm),x)
  return(res)
}


# Name of data files
getNames <- function(fileName){
  tmpName <- gsub(pattern = "\\methylome_", "", basename(fileName))
  fileName <- gsub(pattern = "\\.txt", "",tmpName)
  return(fileName)
}


# Convert all String into U/I/M if there are full text.
statusStringCheck <-  function(file_A){
  list_status <- c("Unmethylated", "Intermediate", "Methylated")
  strTocheckFileA <- utils::head(file_A$status[1])
  if (strTocheckFileA %in% list_status) {
    file_A$status <- str_replace_all(file_A$status,
      pattern = "Unmethylated", replacement = "U")
    file_A$status <- str_replace_all(file_A$status,
      pattern = "Intermediate", replacement = "I")
    file_A$status <- str_replace_all(file_A$status,
      pattern = "Methylated", replacement = "M")
  }
  return(file_A)
}


# saving results- DMR pattern
saveResult<-function(finalDF, name, cytosine, out.dir, current){
  saveFile <- paste0(out.dir, basename(name), cytosine,"_", current, ".csv")
  fwrite(finalDF, file = saveFile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  message("DMR Results saved in: ", saveFile, "\n")

}
# saving results- same pattern
saveSameCyt<-function(samePattern, name, cytosine, out.dir, current){
  saveSamepattern <- paste0(out.dir, name, cytosine, "_", current, ".csv")
  fwrite(samePattern,file = saveSamepattern ,quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
  message("Same Cytosines: ", saveSamepattern, "\n")

}

### * makeDMRmatrix.R


#' @param filepath
#' @param colm
#' @param include.intermediate
#' @param mincov
#' @param nCytosines
#' @importFrom dplyr inner_join
#' @export
#'

# This function will merge (column 6) state calls and (column 7) rc.meth.lvl from all samples into one dataframe
# makes list of 2 dataframes
merge.cols <- function(filepath, colm, include.intermediate) {

  mylist <- list()
  for (l in 1:length(colm)){
    extract <- lapply(filepath, function(k){
      f <- fread(k, header=FALSE, skip=1, select=c(1, 2, 3, colm[l]))
      if (colm[l]==6) {
        if (include.intermediate==TRUE) {
          f[,4] <- ifelse(f[,4] == "U", yes = 0, (ifelse(f[,4] == "I", yes = 0.5, no = 1)))
        } else {
          f[,4] <- ifelse(f[,4] == "U", yes = 0, no = 1)
        }
      }
      colnames(f)[4] <- basename(k)
      return(f)

    })
    df <- Reduce(function(x, y) {
      dplyr::inner_join(x, y, by=c("V1","V2","V3"))
    }, extract)

    mylist[[l]] <- df
  }
  return(mylist)
}

#' Builds a DMR matrix for all samples
#'
#' This function generates a binary matrix, a matrix of recalibrated methylation levels and posterior probabilities for all samples.
#'
#' @param context cytosine context
#' @param samplefiles file containing full path of base level methylation calls, sample names and replicates(optional)
#' @param include.intermediate A logical specifying whether or not the intermediate component should be included in the HMM.
#' @param input.dir input directory containing all region level methylome calls
#' @param out.dir output directory
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table rbindlist
#' @export
#'
makeDMRmatrix <- function(context, samplefiles, input.dir, out.dir, include.intermediate=FALSE) {
  # Read the sample file with filenames
  samplelist <- fread(samplefiles, header=T)
  for (j in  1:length(context)){

    # list all files in the input directory
    extractflist <- list.files(input.dir, pattern=paste0(context[j],".txt"), full.names=TRUE)

    #extract file subsets for construction of DMRmatrix
    if (length(extractflist) != 0){
      mynames <- gsub(paste0("_", context[j], ".txt$"), "", basename(extractflist))
      selectlist <- list()
      message("\nExtracting filenames and matching them....")
      for (a1 in seq_along(mynames)){
        as <- samplelist[grepl(paste0("_",mynames[a1]), samplelist$file),]
        if (NROW(as)==1){
          as$full.path.MethReg <- grep(paste0("/", mynames[a1], "_", context[j], ".txt", sep=""), extractflist, value=TRUE)
          message("\n", basename(as$full.path.MethReg)," found !")
          selectlist[[a1]] <- as
        } else {
          message("\nMultiple files with string match ", mynames[a1]," found !")
        }
      }
      flist <- data.table::rbindlist(selectlist)
      #print(flist)

      # Assign unique names for samples with or without replicate data
      if (!is.null(flist$replicate)) {
        message(paste0("\nRunning context ", context[j], ". Input data with replicates, creating unique sample names...\n"), sep = "")
        flist$name <- paste0(flist$sample,"_", flist$replicate)
      } else {
        flist$name <- flist$sample
      }

      message(paste0("\nNow, constructing DMR matrix for ", context[j]), sep = "")

      # merge samples by Chr coordinates
      #(column 6) state-calls and (column 7) rc.meth.lvl
      mydf <- merge.cols(filepath=flist$full.path.MethReg, include.intermediate=include.intermediate, colm=c(5, 6, 7))

     # list containing state calls
      status.collect <- mydf[[2]]
      # renaming file names with sample names
      for (a in 4:length(colnames(status.collect))) {
        for (n in 1:length(flist$name)) {
          if (colnames(status.collect)[a] == basename(flist$full.path.MethReg)[n]) {
            colnames(status.collect)[a] = flist$name[n]
          }
        }
      }
      # list containing rcmethlvls
      rc.methlevel.collect <- mydf[[3]]
      # renaming file names with sample names
      for (a in 4:length(colnames(rc.methlevel.collect))) {
        for (n in 1:length(flist$name)) {
          if (colnames(rc.methlevel.collect)[a] == basename(flist$full.path.MethReg)[n]) {
            colnames(rc.methlevel.collect)[a] = flist$name[n]
          }
        }
      }

       # list containing postmax
      postMax.collect <- mydf[[1]]
      # renaming file names with sample names
      for (a in 4:length(colnames(postMax.collect))) {
        for (n in 1:length(flist$name)) {
          if (colnames(postMax.collect)[a] == basename(flist$full.path.MethReg)[n]) {
            colnames(postMax.collect)[a] = flist$name[n]
          }
        }
      }

      names(status.collect)[1] <- "seqnames"
      names(status.collect)[2] <- "start"
      names(status.collect)[3] <- "end"

      names(rc.methlevel.collect)[1] <- "seqnames"
      names(rc.methlevel.collect)[2] <- "start"
      names(rc.methlevel.collect)[3] <- "end"

      names(postMax.collect)[1] <- "seqnames"
      names(postMax.collect)[2] <- "start"
      names(postMax.collect)[3] <- "end"

      fwrite(x=status.collect, file=paste0(out.dir,"/", context[j],"_StateCalls.txt"),
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
      fwrite(x=rc.methlevel.collect, file=paste0(out.dir,"/", context[j],"_rcMethlvl.txt"),
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
      fwrite(x=postMax.collect, file=paste0(out.dir,"/", context[j],"_postMax.txt"),
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
      message(paste0("\nDone! \n"), sep = "")
    } else{
      message(paste0("Files for context ",context[j]," do not exist\n"), sep="")
    }
  }
}

### * makeReg.R

#' Regnerates cytosine clusters in genome
#'
#' Uses the extracted cytosines from CfromFASTA and finds clsuters of cytosines along the genome.
#'
#' @param ref.genome
#' @param contexts
#' @param chr
#' @param min.C
#' @param N.boot
#' @param N.sim.C
#' @param fp.rate
#' @param set.tol
#' @param out.dir
#' @param out.name
#' @param makeRegnull
#' @import R.utils
#' @import GenomicRanges
#' @export
#'

makeReg<-function(ref.genome, contexts, chr, min.C, N.boot, N.sim.C, fp.rate, set.tol, out.dir, out.name, makeRegnull)
{

  # Defining function for cleaning objects in the global environment from within the function
  CleanEnvir <- function(pattern)
  {
    rm(list = ls(envir=globalenv())[
      grep(pattern, ls(envir=globalenv()))], envir = globalenv())
  }

  #ref.genome<-out
  #context="CHH"
  #chr=1
  #min.C=5
  #N.boot=10^5
  #fp.rate=0.01
  #set.tol<-0.001

  # This function constructs the regions based on a sequence input
  conReg<-function(seq.in, min.C, win, chr)
  {

    start<-NULL
    stop<-NULL
    index<-rep(0, length(seq.in))

    pb <- txtProgressBar(min = 1, max = length(seq.in), char = "=", style = 3, file = "")
    for (i in seq_len(c(length(seq.in)- c(min.C-1))))
    {
      if (as.numeric(seq.in[i+ c(min.C-1)])-as.numeric(seq.in[i]) <= win){index[i:c(i + c(min.C-1))]<-rep(1, min.C)}

      Sys.sleep(1/length(seq.in)); setTxtProgressBar(pb, i)

    }

    close(pb)

    runs<-rle(index > 0)
    myruns = which(runs$values == TRUE & runs$lengths >= min.C)
    runs.lengths.cumsum = cumsum(runs$lengths)
    ends = runs.lengths.cumsum[myruns]
    newindex = ifelse(myruns>1, myruns-1, 0)
    starts = runs.lengths.cumsum[newindex] + 1
    if (0 %in% newindex) starts = c(1,starts)


    start.pos<-as.numeric(unlist(seq.in[starts]))
    end.pos<-as.numeric(unlist(seq.in[ends]))
    cluster.length<-end.pos - start.pos

    cluster<-cbind(chr, start.pos, end.pos, cluster.length)
    colnames(cluster)<-c("chr","start", "end", "cluster.length")
    cluster<-as.data.frame(cluster, stringsAsFactors = FALSE)
    cluster$chr<-as.character(cluster$chr)
    cluster$start<-as.integer(cluster$start)
    cluster$end<-as.integer(cluster$end)
    cluster$cluster.length<-as.integer(cluster$cluster.length)
    cluster$region<-paste("reg", 1:nrow(cluster), sep="")

    return(cluster)
  }




  cat("- Reading in the data", "\n")
  cat("-------------------------", "\n")
  cat("                         ", "\n")

  # Subselect genome based on chr.id
  #$ref.genome <- ref.genome %>% filter(ref.genome$chr==chr)

  # Making the position vectors
  CHH.pos.plus<-ref.genome[which(ref.genome$strand  == "+" &  ref.genome$context == "CHH"),]$pos
  CHG.pos.plus<-ref.genome[which(ref.genome$strand  == "+" &  ref.genome$context == "CHG"),]$pos
  CG.pos.plus<-ref.genome[which(ref.genome$strand  == "+" &  ref.genome$context == "CG"),]$pos
  CHH.pos.minus<-ref.genome[which(ref.genome$strand  == "-" &  ref.genome$context == "CHH"),]$pos
  CHG.pos.minus<-ref.genome[which(ref.genome$strand  == "-" &  ref.genome$context == "CHG"),]$pos
  CG.pos.minus<-ref.genome[which(ref.genome$strand  == "-" &  ref.genome$context == "CG"),]$pos

  # Determining the number of CHH, CHG and CG sites (+) strand
  N.CHH.obs<-length(CHH.pos.plus)
  N.CHG.obs<-length(CHG.pos.plus)
  N.CG.obs<-length(CG.pos.plus)


  if (N.sim.C == "all")
  {
    chr.length<-max(ref.genome$pos)
    N.CHH<-N.CHH.obs
    N.CHG<-N.CHG.obs
    N.CG<-N.CG.obs
    set.seed(123)
    rp.CHH<-sample(1:chr.length, size=N.CHH, replace = FALSE)
    rp.CHG<-sample(1:chr.length, size=N.CHG, replace = FALSE)
    rp.CG<-sample(1:chr.length, size=N.CG, replace = FALSE)
  }else{
    if (is.numeric(N.sim.C) == FALSE){stop("N.sim.C needs to be numeric")}else{
      chr.length<-max(ref.genome$pos)
      scale.factor<-N.sim.C/c(N.CHH.obs + N.CHG.obs + N.CG.obs)
      chr.length<-round(chr.length*scale.factor,0)
      N.CHH<-round(N.CHH.obs*scale.factor,0)
      N.CHG<-round(N.CHG.obs*scale.factor,0)
      N.CG<-round(N.CG.obs*scale.factor,0)
      set.seed(123)
      rp.CHH<-sample(1:chr.length, size=N.CHH, replace = FALSE)
      rp.CHG<-sample(1:chr.length, size=N.CHG, replace = FALSE)
      rp.CG<-sample(1:chr.length, size=N.CG, replace = FALSE)
    }
  }


  # Remove this large data. It's no longer needed
  CleanEnvir(pattern = "ref.genome")

  cat("- Simulating cytosines", "\n")
  cat("-------------------------", "\n")

  rp.CHH<-sort(rp.CHH, method="quick")
  rp.CHG<-sort(rp.CHG, method="quick")
  rp.CG<-sort(rp.CG, method="quick")

  CHH.lost<-0
  CHG.lost<-0
  CG.lost<-0

  tolerance = set.tol + 0.00001

  while (tolerance > set.tol)
  {
    ### CHH
    # Draw random positions (rp) of CHH sites
    rp.CHH<-c(rp.CHH, sample(1:chr.length, size= CHH.lost, replace =FALSE))
    if  (length(rp.CHH) > N.CHH){rp.CHH<-sample(rp.CHH, size = N.CHH, replace=F)}
    rp.CHH<-sort(rp.CHH)

    ### CHG
    # Draw random positions (rp) of CHG sites
    # Rules 1 to 3
    rp.CHG<-c(rp.CHG, sample(1:chr.length, size= CHG.lost, replace =FALSE))
    if  (length(rp.CHG) > N.CHG){rp.CHG<-sample(rp.CHG, size = N.CHG, replace=F)}
    rp.CHG<-sort(rp.CHG)
    rp.CHG<-rp.CHG[which(diff(rp.CHG, lag=1) >= 3)]

    # Rule 4
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, rp.CHG))]

    # Rule 5
    # no action

    # Rule 6
    # no action

    # Rule 7
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, c(rp.CHG+1)))]

    # Rule 8
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, c(rp.CHG+2)))]


    ### CG
    # Draw random positions (rp) of CG sites
    # Rules 1 to 3
    rp.CG<-c(rp.CG, sample(1:chr.length, size= CG.lost, replace =FALSE))
    if  (length(rp.CG) > N.CG){rp.CG<-sample(rp.CG, size = N.CG, replace=F)}
    rp.CG<-sort(rp.CG)
    rp.CG<-rp.CG[which(diff(rp.CG, lag=1) >= 2)]

    # Rule 4
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, rp.CG))]

    # Rule 5
    rp.CHH.temp<-rp.CHH+1
    rp.CHH.temp<-rp.CHH.temp[which(is.element(rp.CHH.temp, rp.CG))]
    rp.CHH<-rp.CHH[which(!is.element(c(rp.CHH+1), rp.CG))]
    rp.CHG<-c(rp.CHG, c(rp.CHH.temp-1))
    rp.CHG<-unique(rp.CHG)
    if (length(rp.CHG) > N.CHG){rp.CHG<-sample(rp.CHG, size = N.CHG, replace =F)}
    rp.CHG<-sort(rp.CHG)
    rp.CHG<-rp.CHG[which(diff(rp.CHG, lag=1) >= 3)]

    # Rule 6
    # no action

    # Rule 7
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, c(rp.CG+1)))]

    # Rule 8
    rp.CHG<-rp.CHG[which(!is.element(rp.CHG, rp.CG))]

    # Rule 9
    # no action

    # Rule 10
    rp.CHG.temp<-rp.CHG+2
    rp.CHG.temp<-rp.CHG.temp[which(is.element(rp.CHG.temp, rp.CG))]
    rp.CHG<-rp.CHG[which(!is.element(c(rp.CHG+2), rp.CG))]
    rp.CG<-c(rp.CG, c(rp.CHG.temp-2))
    rp.CG<-unique(rp.CG)
    if (length(rp.CG) > N.CG){rp.CG<-sample(rp.CG, size = N.CG, replace =F)}
    rp.CG<-sort(rp.CG)
    rp.CG<-rp.CG[which(diff(rp.CG, lag=1) >= 2)]

    # Rule 11
    rp.CHG<-rp.CHG[which(!is.element(rp.CHG, c(rp.CG+1)))]

    # Counting the lost sites
    CHH.lost<-N.CHH - length(rp.CHH)
    CHG.lost<-N.CHG - length(rp.CHG)
    CG.lost<-N.CG - length(rp.CG)
    N.lost<-CHH.lost + CHG.lost + CG.lost

    tolerance<-N.lost/(N.CG + N.CHG + N.CHH)

    cat(tolerance, "\n")

  }

  cat("                         ", "\n")
  cat("- Converged", "\n")
  cat("-------------------------", "\n")

  sim.out<-data.frame(c(N.CG.obs, N.CHG.obs, N.CHH.obs), c(length(rp.CG), length(rp.CHG), length(rp.CHH)))
  sim.out<-rbind(sim.out, colSums(sim.out))
  rownames(sim.out)<-c("CG", "CHG", "CHH", "Total")
  sim.out$percent<-sim.out[,2]/sim.out[,1]*100
  colnames(sim.out)<-c("N.observed", "N.simulated", "Percent")

  print(sim.out)

  cat("                         ", "\n")

  cat("- Calling regions:", "\n")
  cat("-------------------------", "\n")
  cat("                         ", "\n")


  for (icontext in 1:length(contexts))
  {
    context.temp<-contexts[icontext]


    if (context.temp == "C")
    {
      sim.geno.out<-c(rp.CG, rp.CHG, rp.CHH)
      sim.geno.out<-sort(sim.geno.out, method="quick")

      obs.geno.plus<-sort(c(CG.pos.plus, CHG.pos.plus, CHH.pos.plus))
      obs.geno.minus<-sort(c(CG.pos.minus, CHG.pos.minus, CHH.pos.minus))

      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] -sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))

      cat("Building from: C (+) strand", "\n")
      reg.obs.plus<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      cat("Building from: C (-) strand", "\n")
      reg.obs.minus<-conReg(seq.in=obs.geno.minus, min.C=min.C, win=null.dist, chr=chr)
      df <- data.frame(id=c(rep("+", nrow(reg.obs.plus)), rep("-", nrow(reg.obs.minus))),
                       start=c(reg.obs.plus[,2], reg.obs.minus[,2]), end=c(reg.obs.plus[,3], reg.obs.minus[,3]))
      gr <- GRanges(seqnames = rep(1,nrow(df)),IRanges(start = df$start, end = df$end))
      reg.obs<-as.data.frame(reduce(gr))[,1:3]
      reg.obs[,1]<-rep(chr, nrow(reg.obs))
      reg.obs$cluster.length <- reg.obs$end - reg.obs$start
      reg.obs$region<-paste("reg", 1:nrow(reg.obs), sep="")
      colnames(reg.obs)[1]<-"chr"

      if ("C" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("Building: C null regions", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
      }
      if (!is.element("C", contexts[which(makeRegnull == TRUE)]))
      {
        cat("C NULL omitted", "\n")
      }
    }


    if (context.temp == "CG")
    {
      sim.geno.out<-rp.CG

      obs.geno.plus<-sort(CG.pos.plus)
      obs.geno.minus<-sort(CG.pos.minus)

      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] - sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))

      cat("Building from: from CG (+) strand", "\n")
      reg.obs<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      reg.obs$start<-reg.obs$start - 1
      reg.obs$end<-reg.obs$end + 1
      reg.obs$cluster.length<-reg.obs$end - reg.obs$start

      if ("CG" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("Building: CG null regions", "\n")
        cat("", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
        reg.sim$start<-reg.sim$start - 1
        reg.sim$end<-reg.sim$end + 1
        reg.sim$cluster.length<-reg.sim$end - reg.sim$start
      }
      if (!is.element("CG", contexts[which(makeRegnull == TRUE)]))
      {
        cat("CG NULL omitted", "\n")
      }
    }

    if (context.temp == "CHG")
    {
      sim.geno.out<-rp.CHG
      obs.geno.plus<-sort(CHG.pos.plus)
      obs.geno.minus<-sort(CHG.pos.minus)

      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] - sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))

      cat("Building from: CHG (+) strand", "\n")
      reg.obs<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      reg.obs$start<-reg.obs$start - 2
      reg.obs$end<-reg.obs$end + 2
      reg.obs$cluster.length<-reg.obs$end - reg.obs$start

      if ("CHG" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("Building: CHG null regions", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
        reg.sim$start<-reg.sim$start - 2
        reg.sim$end<-reg.sim$end + 2
        reg.sim$cluster.length<-reg.sim$end - reg.sim$start
      }
      if (!is.element("CHG", contexts[which(makeRegnull == TRUE)]))
      {
        cat("CHG NULL omitted", "\n")
      }
    }

    if (context.temp == "CHH")
    {
      sim.geno.out<-rp.CHH
      obs.geno.plus<-sort(CHH.pos.plus)
      obs.geno.minus<-sort(CHH.pos.minus)

      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] - sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))

      cat("Building from: CHH (+) strand", "\n")
      reg.obs.plus<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      cat("Building from: CHH (-) strand", "\n")
      reg.obs.minus<-conReg(seq.in=obs.geno.minus, min.C=min.C, win=null.dist, chr=chr)
      df <- data.frame(id=c(rep("+", nrow(reg.obs.plus)), rep("-", nrow(reg.obs.minus))),
                       start=c(reg.obs.plus[,2], reg.obs.minus[,2]), end=c(reg.obs.plus[,3], reg.obs.minus[,3]))
      gr <- GRanges(seqnames = rep(1,nrow(df)),IRanges(start = df$start, end = df$end))
      reg.obs<-as.data.frame(reduce(gr))[,1:3]
      reg.obs[,1]<-rep(chr, nrow(reg.obs))
      reg.obs$cluster.length <- reg.obs$end - reg.obs$start
      reg.obs$region<-paste("reg", 1:nrow(reg.obs), sep="")
      colnames(reg.obs)[1]<-"chr"


      if ("CHH" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("CHH null regions", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
      }
      if (!is.element("C", contexts[which(makeRegnull == TRUE)]))
      {
        cat("CHH NULL omitted", "\n")
      }
    }

    if (context.temp %in% contexts[which(makeRegnull == TRUE)])
    {
      output<-list(reg.obs, reg.sim, chr, context.temp)
      names(output)<-c("reg.obs", "reg.sim", "chr", "context")
      dput(output, paste(out.dir,"/",out.name,"_regions_","chr", chr, "_", context.temp, ".Rdata", sep=""))
    }
    if (!is.element(context.temp, contexts[which(makeRegnull == TRUE)]))
    {
      reg.sim<-NA
      output<-list(reg.obs, reg.sim, chr, context.temp)
      names(output)<-c("reg.obs", "reg.sim", "chr", "context")
      dput(output, paste(out.dir,"/",out.name, "_regions_","chr", chr, "_", context.temp, ".Rdata", sep=""))
    }

  } # End of context loop

  # Removing specific objects from global environemnt
  CleanEnvir(pattern = "reg.obs.plus")
  CleanEnvir(pattern = "reg.obs.minus")
  CleanEnvir(pattern = "obs.geno.plus")
  CleanEnvir(pattern = "obs.geno.minus")




}

### * makeRegScaffolds.R




makeReg<-function(ref.genome, contexts, chr, min.C, N.boot, N.sim.C, fp.rate, set.tol, out.dir, out.name, makeRegnull)
{
  
  # Defining function for cleaning objects in the global environment from within the function
  CleanEnvir <- function(pattern)
  {
    rm(list = ls(envir=globalenv())[
      grep(pattern, ls(envir=globalenv()))], envir = globalenv())
  }
  
  #ref.genome<-out
  #context="CHH"
  #chr=1
  #min.C=5
  #N.boot=10^5
  #fp.rate=0.01
  #set.tol<-0.001
  
  # This function constructs the regions based on a sequence input
  conReg<-function(seq.in, min.C, win, chr)
  {
    
    start<-NULL
    stop<-NULL
    index<-rep(0, length(seq.in))
    
    pb <- txtProgressBar(min = 1, max = length(seq.in), char = "=", style = 3, file = "")
    for (i in seq_len(c(length(seq.in)- c(min.C-1))))
    {
      if (as.numeric(seq.in[i+ c(min.C-1)])-as.numeric(seq.in[i]) <= win){index[i:c(i + c(min.C-1))]<-rep(1, min.C)}
      
      Sys.sleep(1/length(seq.in)); setTxtProgressBar(pb, i)
      
    }
    
    close(pb)
    
    runs<-rle(index > 0)
    myruns = which(runs$values == TRUE & runs$lengths >= min.C)
    runs.lengths.cumsum = cumsum(runs$lengths)
    ends = runs.lengths.cumsum[myruns]
    newindex = ifelse(myruns>1, myruns-1, 0)
    starts = runs.lengths.cumsum[newindex] + 1
    if (0 %in% newindex) starts = c(1,starts)
    
    
    start.pos<-as.numeric(unlist(seq.in[starts]))
    end.pos<-as.numeric(unlist(seq.in[ends]))
    cluster.length<-end.pos - start.pos
    
    cluster<-cbind(chr, start.pos, end.pos, cluster.length)
    colnames(cluster)<-c("chr","start", "end", "cluster.length")
    cluster<-as.data.frame(cluster, stringsAsFactors = FALSE)
    cluster$chr<-as.character(cluster$chr)
    cluster$start<-as.integer(cluster$start)
    cluster$end<-as.integer(cluster$end)
    cluster$cluster.length<-as.integer(cluster$cluster.length)
    cluster$region<-paste("reg", 1:nrow(cluster), sep="")
    
    return(cluster)
  }
  
  
  
  
  cat("- Reading in the data", "\n")
  cat("-------------------------", "\n")
  cat("                         ", "\n")
  
  # Subselect genome based on chr.id
  #$ref.genome <- ref.genome %>% filter(ref.genome$chr==chr)
  
  # Making the position vectors
  CHH.pos.plus<-ref.genome[which(ref.genome$strand  == "+" &  ref.genome$context == "CHH"),]$pos
  CHG.pos.plus<-ref.genome[which(ref.genome$strand  == "+" &  ref.genome$context == "CHG"),]$pos
  CG.pos.plus<-ref.genome[which(ref.genome$strand  == "+" &  ref.genome$context == "CG"),]$pos
  CHH.pos.minus<-ref.genome[which(ref.genome$strand  == "-" &  ref.genome$context == "CHH"),]$pos
  CHG.pos.minus<-ref.genome[which(ref.genome$strand  == "-" &  ref.genome$context == "CHG"),]$pos
  CG.pos.minus<-ref.genome[which(ref.genome$strand  == "-" &  ref.genome$context == "CG"),]$pos
  
  # Determining the number of CHH, CHG and CG sites (+) strand
  N.CHH.obs<-length(CHH.pos.plus)
  N.CHG.obs<-length(CHG.pos.plus)
  N.CG.obs<-length(CG.pos.plus)
  
  
  if (N.sim.C == "all")
  {
    chr.length<-max(ref.genome$pos)
    N.CHH<-N.CHH.obs
    N.CHG<-N.CHG.obs
    N.CG<-N.CG.obs
    set.seed(123)
    rp.CHH<-sample(1:chr.length, size=N.CHH, replace = FALSE)
    rp.CHG<-sample(1:chr.length, size=N.CHG, replace = FALSE)
    rp.CG<-sample(1:chr.length, size=N.CG, replace = FALSE)
  }else{
    if (is.numeric(N.sim.C) == FALSE){stop("N.sim.C needs to be numeric")}else{
      chr.length<-max(ref.genome$pos)
      scale.factor<-N.sim.C/c(N.CHH.obs + N.CHG.obs + N.CG.obs) 
      chr.length<-round(chr.length*scale.factor,0)
      N.CHH<-round(N.CHH.obs*scale.factor,0)
      N.CHG<-round(N.CHG.obs*scale.factor,0)
      N.CG<-round(N.CG.obs*scale.factor,0)
      set.seed(123)
      rp.CHH<-sample(1:chr.length, size=N.CHH, replace = FALSE)
      rp.CHG<-sample(1:chr.length, size=N.CHG, replace = FALSE)
      rp.CG<-sample(1:chr.length, size=N.CG, replace = FALSE)
    }
  }
  
  
  # Remove this large data. It's no longer needed
  CleanEnvir(pattern = "ref.genome")
  tryCatch(
    withTimeout ({
      Sys.sleep(2)
      
      cat("- Simulating cytosines", "\n")
      cat("-------------------------", "\n")
      
      rp.CHH<-sort(rp.CHH, method="quick")
      rp.CHG<-sort(rp.CHG, method="quick")
      rp.CG<-sort(rp.CG, method="quick")
      
      CHH.lost<-0
      CHG.lost<-0
      CG.lost<-0
      
      tolerance = set.tol + 0.00001
      
      while (tolerance > set.tol)
      {
        ### CHH
        # Draw random positions (rp) of CHH sites
        rp.CHH<-c(rp.CHH, sample(1:chr.length, size= CHH.lost, replace =FALSE))
        if  (length(rp.CHH) > N.CHH){rp.CHH<-sample(rp.CHH, size = N.CHH, replace=F)}
        rp.CHH<-sort(rp.CHH)
        
        ### CHG 
        # Draw random positions (rp) of CHG sites
        # Rules 1 to 3
        rp.CHG<-c(rp.CHG, sample(1:chr.length, size= CHG.lost, replace =FALSE))
        if  (length(rp.CHG) > N.CHG){rp.CHG<-sample(rp.CHG, size = N.CHG, replace=F)}
        rp.CHG<-sort(rp.CHG)
        rp.CHG<-rp.CHG[which(diff(rp.CHG, lag=1) >= 3)]
        
        # Rule 4
        rp.CHH<-rp.CHH[which(!is.element(rp.CHH, rp.CHG))]
        
        # Rule 5
        # no action
        
        # Rule 6
        # no action
        
        # Rule 7
        rp.CHH<-rp.CHH[which(!is.element(rp.CHH, c(rp.CHG+1)))]
        
        # Rule 8
        rp.CHH<-rp.CHH[which(!is.element(rp.CHH, c(rp.CHG+2)))]
        
        
        ### CG 
        # Draw random positions (rp) of CG sites
        # Rules 1 to 3
        rp.CG<-c(rp.CG, sample(1:chr.length, size= CG.lost, replace =FALSE))
        if  (length(rp.CG) > N.CG){rp.CG<-sample(rp.CG, size = N.CG, replace=F)}
        rp.CG<-sort(rp.CG)
        rp.CG<-rp.CG[which(diff(rp.CG, lag=1) >= 2)]
        
        # Rule 4
        rp.CHH<-rp.CHH[which(!is.element(rp.CHH, rp.CG))]
        
        # Rule 5
        rp.CHH.temp<-rp.CHH+1
        rp.CHH.temp<-rp.CHH.temp[which(is.element(rp.CHH.temp, rp.CG))] 
        rp.CHH<-rp.CHH[which(!is.element(c(rp.CHH+1), rp.CG))]
        rp.CHG<-c(rp.CHG, c(rp.CHH.temp-1))
        rp.CHG<-unique(rp.CHG)
        if (length(rp.CHG) > N.CHG){rp.CHG<-sample(rp.CHG, size = N.CHG, replace =F)}
        rp.CHG<-sort(rp.CHG)
        rp.CHG<-rp.CHG[which(diff(rp.CHG, lag=1) >= 3)]
        
        # Rule 6
        # no action
        
        # Rule 7
        rp.CHH<-rp.CHH[which(!is.element(rp.CHH, c(rp.CG+1)))]
        
        # Rule 8
        rp.CHG<-rp.CHG[which(!is.element(rp.CHG, rp.CG))]
        
        # Rule 9
        # no action
        
        # Rule 10
        rp.CHG.temp<-rp.CHG+2
        rp.CHG.temp<-rp.CHG.temp[which(is.element(rp.CHG.temp, rp.CG))] 
        rp.CHG<-rp.CHG[which(!is.element(c(rp.CHG+2), rp.CG))]
        rp.CG<-c(rp.CG, c(rp.CHG.temp-2))
        rp.CG<-unique(rp.CG)
        if (length(rp.CG) > N.CG){rp.CG<-sample(rp.CG, size = N.CG, replace =F)}
        rp.CG<-sort(rp.CG)
        rp.CG<-rp.CG[which(diff(rp.CG, lag=1) >= 2)]
        
        # Rule 11
        rp.CHG<-rp.CHG[which(!is.element(rp.CHG, c(rp.CG+1)))]
        
        # Counting the lost sites
        CHH.lost<-N.CHH - length(rp.CHH)
        CHG.lost<-N.CHG - length(rp.CHG)
        CG.lost<-N.CG - length(rp.CG)
        N.lost<-CHH.lost + CHG.lost + CG.lost
        
        tolerance<-N.lost/(N.CG + N.CHG + N.CHH)
        
        cat(tolerance, "\n")
        
      }
      cat("[Done]\n")},
      timeout=10),
    TimeoutException=function(ex) {cat("[skipped due to timeout]\n")})
  
  
  cat("                         ", "\n")
  cat("- Converged", "\n")
  cat("-------------------------", "\n")
  
  sim.out<-data.frame(c(N.CG.obs, N.CHG.obs, N.CHH.obs), c(length(rp.CG), length(rp.CHG), length(rp.CHH)))
  sim.out<-rbind(sim.out, colSums(sim.out))
  rownames(sim.out)<-c("CG", "CHG", "CHH", "Total")
  sim.out$percent<-sim.out[,2]/sim.out[,1]*100
  colnames(sim.out)<-c("N.observed", "N.simulated", "Percent")
  
  print(sim.out)
  
  cat("                         ", "\n")
  
  cat("- Calling regions:", "\n")
  cat("-------------------------", "\n")
  cat("                         ", "\n")
  
  
  for (icontext in 1:length(contexts))
  {
    context.temp<-contexts[icontext]
    
    
    if (context.temp == "C")
    {
      sim.geno.out<-c(rp.CG, rp.CHG, rp.CHH)
      sim.geno.out<-sort(sim.geno.out, method="quick")
      
      obs.geno.plus<-sort(c(CG.pos.plus, CHG.pos.plus, CHH.pos.plus))
      obs.geno.minus<-sort(c(CG.pos.minus, CHG.pos.minus, CHH.pos.minus))
      
      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] -sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))
      
      cat("Building from: C (+) strand", "\n")
      reg.obs.plus<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      cat("Building from: C (-) strand", "\n")
      reg.obs.minus<-conReg(seq.in=obs.geno.minus, min.C=min.C, win=null.dist, chr=chr)
      df <- data.frame(id=c(rep("+", nrow(reg.obs.plus)), rep("-", nrow(reg.obs.minus))), 
                       start=c(reg.obs.plus[,2], reg.obs.minus[,2]), end=c(reg.obs.plus[,3], reg.obs.minus[,3]))
      gr <- GRanges(seqnames = rep(1,nrow(df)),IRanges(start = df$start, end = df$end))
      reg.obs<-as.data.frame(reduce(gr))[,1:3]
      reg.obs[,1]<-rep(chr, nrow(reg.obs))
      reg.obs$cluster.length <- reg.obs$end - reg.obs$start
      reg.obs$region<-paste("reg", 1:nrow(reg.obs), sep="")
      colnames(reg.obs)[1]<-"chr"
      
      if ("C" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("Building: C null regions", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
      }
      if (!is.element("C", contexts[which(makeRegnull == TRUE)]))
      {
        cat("C NULL omitted", "\n")
      }
    }
    
    
    if (context.temp == "CG")
    {
      sim.geno.out<-rp.CG
      
      obs.geno.plus<-sort(CG.pos.plus)
      obs.geno.minus<-sort(CG.pos.minus)
      
      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] - sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))
      
      cat("Building from: from CG (+) strand", "\n")
      reg.obs<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      reg.obs$start<-reg.obs$start - 1
      reg.obs$end<-reg.obs$end + 1
      reg.obs$cluster.length<-reg.obs$end - reg.obs$start
      
      if ("CG" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("Building: CG null regions", "\n")
        cat("", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
        reg.sim$start<-reg.sim$start - 1
        reg.sim$end<-reg.sim$end + 1
        reg.sim$cluster.length<-reg.sim$end - reg.sim$start
      }
      if (!is.element("CG", contexts[which(makeRegnull == TRUE)]))
      {
        cat("CG NULL omitted", "\n")
      }
    }
    
    if (context.temp == "CHG")
    {
      sim.geno.out<-rp.CHG
      obs.geno.plus<-sort(CHG.pos.plus)
      obs.geno.minus<-sort(CHG.pos.minus)
      
      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] - sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))
      
      cat("Building from: CHG (+) strand", "\n")
      reg.obs<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      reg.obs$start<-reg.obs$start - 2
      reg.obs$end<-reg.obs$end + 2
      reg.obs$cluster.length<-reg.obs$end - reg.obs$start
      
      if ("CHG" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("Building: CHG null regions", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
        reg.sim$start<-reg.sim$start - 2
        reg.sim$end<-reg.sim$end + 2
        reg.sim$cluster.length<-reg.sim$end - reg.sim$start
      }
      if (!is.element("CHG", contexts[which(makeRegnull == TRUE)]))
      {
        cat("CHG NULL omitted", "\n")
      }
    }
    
    if (context.temp == "CHH")
    {
      sim.geno.out<-rp.CHH
      obs.geno.plus<-sort(CHH.pos.plus)
      obs.geno.minus<-sort(CHH.pos.minus)
      
      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] - sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))
      
      cat("Building from: CHH (+) strand", "\n")
      reg.obs.plus<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      cat("Building from: CHH (-) strand", "\n")
      reg.obs.minus<-conReg(seq.in=obs.geno.minus, min.C=min.C, win=null.dist, chr=chr)
      df <- data.frame(id=c(rep("+", nrow(reg.obs.plus)), rep("-", nrow(reg.obs.minus))), 
                       start=c(reg.obs.plus[,2], reg.obs.minus[,2]), end=c(reg.obs.plus[,3], reg.obs.minus[,3]))
      gr <- GRanges(seqnames = rep(1,nrow(df)),IRanges(start = df$start, end = df$end))
      reg.obs<-as.data.frame(reduce(gr))[,1:3]
      reg.obs[,1]<-rep(chr, nrow(reg.obs))
      reg.obs$cluster.length <- reg.obs$end - reg.obs$start
      reg.obs$region<-paste("reg", 1:nrow(reg.obs), sep="")
      colnames(reg.obs)[1]<-"chr"
      
      
      if ("CHH" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("CHH null regions", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
      }
      if (!is.element("C", contexts[which(makeRegnull == TRUE)]))
      {
        cat("CHH NULL omitted", "\n")
      }
    }
    
    if (context.temp %in% contexts[which(makeRegnull == TRUE)])
    {
      output<-list(reg.obs, reg.sim, chr, context.temp)
      names(output)<-c("reg.obs", "reg.sim", "chr", "context")
      dput(output, paste(out.dir, out.name, "_regions_","chr", chr, "_", context.temp, ".Rdata", sep=""))
    }
    if (!is.element(context.temp, contexts[which(makeRegnull == TRUE)]))
    {
      reg.sim<-NA
      output<-list(reg.obs, reg.sim, chr, context.temp)
      names(output)<-c("reg.obs", "reg.sim", "chr", "context")
      dput(output, paste(out.dir, out.name, "_regions_","chr", chr, "_", context.temp, ".Rdata", sep=""))
    }
    
  } # End of context loop
  
  # Removing specific objects from global environemnt
  CleanEnvir(pattern = "reg.obs.plus")
  CleanEnvir(pattern = "reg.obs.minus")
  CleanEnvir(pattern = "obs.geno.plus")
  CleanEnvir(pattern = "obs.geno.minus")
  
  
  
  
}   

### * runMethimpute.R

#' Run Methimpute for Regions
#'
#' This function runs a HMM model on identified Cytosine clusters
#'
#' @param samplefiles a text file containing path to samples and sample names, replicate info
#' @param Regionfiles output of makeReg, containing coordinates of cytosine clusters/regions
#' @param genome genome label for .e.g Arabidopsis
#' @param context cytosine context
#' @param out.dir output directory
#' @param include.intermediate A logical specifying wheter or not the intermediate component should be included in the HMM.By default it is set as FALSE.
#' @param mincov Minimum read coverage over cytosines
#' @param nCytosines Minimum number of cytsoines
#' @importFrom data.table fread
#' @import GenomicRanges
#' @export
#'
runMethimputeRegions <- function(samplefiles, Regionfiles, genome, context, out.dir, include.intermediate=FALSE, mincov=0, nCytosines=0) {
  df.obs <- list()
  df.sim <- list()
  merge.list <- vector(mode="list")

  # Read the sample file with filenames and file paths
  filelist <- data.table::fread(samplefiles, header=TRUE)

  for (j in 1:length(context)){
    Regfiles <- list.files(Regionfiles, pattern=paste0("_", context[j], ".Rdata"), full.names = TRUE)
    #print(Regfiles)
    cat(paste0("Reading Region files. Merging individual chr data for context ", context[j], " ...\n"), sep="")
    cat("\n")
    for (k1 in 1:length(Regfiles)){
      f.file <- dget(Regfiles[k1])
      if (NROW(f.file$reg.obs)==0) {
        cat(paste0("Empty file ", basename(Regfiles[k1]), " ...\n"), sep="")
      } else {
        df.obs[[k1]] <- as.data.frame(f.file$reg.obs)
        df.sim[[k1]] <- as.data.frame(f.file$reg.sim)
      }
    }
    outlist <- list(reg.obs=do.call(rbind,df.obs),
                    reg.sim=do.call(rbind,df.sim),
                    context=context[j])

    regMerged <- paste(out.dir, "/",genome,"_Merged_Regions_", context[j], ".Rdata", sep="")
    dput(outlist, regMerged)

    for (k2 in 1:length(filelist$file)){
      methfn <- gsub(".*methylome_|\\_All.txt$", "", filelist$file[k2])
      cat(paste0("Now running file: ", methfn, " for context ", context[j], " ...\n"), sep="")
      regions.out <- makeMethimpute(
        df=filelist$file[k2],
        context=context[j],
        refRegion=regMerged,
        fit.plot=FALSE,
        include.intermediate=include.intermediate,
        probability="constrained",
        out.dir=out.dir,
        fit.name=paste0(methfn, "_", context[j]),
        name=methfn,
        nCytosines=nCytosines,
        mincov=mincov)
    }
  }
}

#' @importFrom seqinr read.fasta
#' @importFrom seqinr getLength
#' @import GenomicRanges
#' @import IRanges
#' @export
#'
binGenome <- function(fasta, win, step, genome, out.dir){
  val <- c()
  cat("Extracting chromosomes\n")
  chr.ext <- list.files(fasta, pattern="fa|fasta.gz$|fa.gz$", include.dirs = FALSE, full.names=TRUE)
  if (length(chr.ext)==0) {
    stop ("Empty folder!")
  } else {
    for (x in seq_along(chr.ext)){
      f.name <- gsub("Arabidopsis_thaliana.TAIR10.dna.chromosome.|\\.fa|\\.fa.gz|\\.fasta.gz$", "", basename(chr.ext[x]))
      f <- seqinr::read.fasta(chr.ext[x])
      cat(paste0("Chr: ", f.name, " \n"))
      if (length(f)>1) {
        stop ("Exiting... Multi FASTA file detected! Please, supply individual fasta files\n") }
      else {
        val[x] <- seqinr::getLength(f)
        names(val)[x] <- f.name
      }
    }
  }

  #gr <- GRanges(seqnames=gsub("chr", "", names(chrfile)), ranges=IRanges(start=1, end=chrfile))
  gr <- GRanges(seqnames=names(val), ranges=IRanges(start=1, end=val))

  binned.g <- slidingWindows(gr, width = win, step = step)
  d <- data.frame(unlist(binned.g))
  names(d)[1] <- "chr"
  names(d)[2] <- "start"
  names(d)[3] <- "end"
  names(d)[4] <- "cluster.length"
  cat(paste0("Binning genome with windows of: ", win, " bp and step-size of: ", step, " bp\n"), sep = "")
  new <- list(d)
  names(new) <- "reg.obs"
  out.name <- paste0(out.dir, "/", genome,"_Win", win, "_Step", step, ".Rdata", sep="")
  dput(new, out.name)
}

#' Run Methimpute on binned genome
#'
#' this function runs a HMM model on a genome binned using a sliding/non-sliding window approach
#' @param samplefiles a text file containing path to samples and sample names, replicate info
#' @param context cytosine context
#' @param include.intermediate A logical specifying wheter or not the intermediate component should be included in the HMM.By default it is set as FALSE.
#' @param mincov Minimum read coverage over cytosines
#' @param nCytosines Minimum number of cytsoines
#' @param out.dir output directory
#' @param fasta path to genome fasta files
#' @param win window size
#' @param step window step-size
#' @param genome genome label for .e.g Arabidopsis
#' @importFrom  data.table fread
#' @export
#'
#'
runMethimputeGrid <- function(out.dir, fasta, win, step, genome, samplefiles, context, mincov, include.intermediate=FALSE, nCytosines){
  binGenome(fasta, win, step,genome, out.dir)
  merge.list <- vector(mode="list")
  filelist <- data.table::fread(samplefiles, header=TRUE)
  for (j in 1:length(context)){
    for (i in 1:length(filelist$file)){
      methfn <- gsub(".*methylome_|\\.txt|_All.txt$", "", filelist$file[i])
      cat(paste0("Running file: ",methfn," for context: ",context[j],"\n"), sep = "")
      grid.out <- makeMethimpute(
        df=filelist$file[i],
        context=context[j],
        refRegion=paste0(out.dir,"/",genome,"_Win",win,"_Step",step,".Rdata",sep=""),
        fit.plot=FALSE,
        include.intermediate=include.intermediate,
        probability="constrained",
        out.dir=out.dir,
        fit.name=paste0(methfn, "_", context[j]),
        name=methfn,
        nCytosines=nCytosines,
        mincov=mincov)
    }
  }
}
