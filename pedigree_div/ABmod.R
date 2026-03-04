#This is a modified version of the AlphaBeta package

startTime <- function(...) {
  x <- paste0(..., collapse='')
  message(x, appendLF=FALSE)
  st <- proc.time()
  return(st)
}
stopTime <- function(st) {
  sttime <- proc.time() - st
  message("Total time: ", round(sttime[3],2), "s")
}


# running join
rc.meth.lvl <- function(nodelist, cytosine, posteriorMaxFilter){

      #inputCheck(genTable, cytosine, posteriorMaxFilter)
      genTable <- fread(nodelist)
      #---------------------------Filter based on WGBS
      genTable <- genTable %>% filter(genTable$meth=="Y")
      #---------------------------
      list.rc<-bplapply(genTable$filename,cytosine=cytosine,posteriorMaxFilter= posteriorMaxFilter,
        genTable= genTable, rcRun, BPPARAM = multicoreParam)

      rclvl <- RCsaveResult(list.rc,cytosine,posteriorMaxFilter)

      return(rclvl)

}


rcRun <- function(filename, cytosine, posteriorMaxFilter, genTable){
      file <- RC.dataRead(filename,cytosine, posteriorMaxFilter)
      name <-  getNames(filename,genTable)
      mean.rc <-floorDec(as.numeric(mean(file$rc.meth.lvl)),5)
      res<-list(name,mean.rc)
      return(res)
}


RCsaveResult<-function(list.rc,cytosine,posteriorMaxFilter){
      tmp_dmr <- data.frame(matrix(ncol = 3, nrow =1 ))
      x <- c("Sample_name", "context", "rc.meth.lvls")
      colnames(tmp_dmr) <- x
      mainRC <- tmp_dmr
      for (nlist in seq_len(length(list.rc))){
        tmp_dmr$Sample_name<-list.rc[[nlist]][[1]]$node
        tmp_dmr$context<-cytosine
        tmp_dmr$rc.meth.lvls <-list.rc[[nlist]][[2]]
        mainRC<-rbind(mainRC,tmp_dmr)
        tmp_dmr<-NULL
      }
      mainRC <- mainRC[!(rowSums(x = is.na(x = mainRC)) == ncol(x = mainRC)),]
      rownames(mainRC) <- NULL
      mainRC %>% mutate_if(is.factor, as.character) -> mainRC
      mainRC <- mainRC[mixedorder(mainRC$Sample_name), ]
      return(mainRC)

}



plotPedigree <- function(nodelist, edgelist, sampling.design, out.pdf=NULL, output.dir=NULL,
                         plot.width=11, plot.height=8, vertex.label=NULL, vertex.size=12, aspect.ratio=2.5)
  {
    samples <- fread(nodelist, header=TRUE, skip=0, select = c(2,3,4))
    edges   <- fread(edgelist, header=TRUE, skip=0)

    #graph layout from edges and nodes
    g <- graph_from_data_frame(edges, directed=TRUE, vertices=samples)

    #put root vertex on top with flip.y=TRUE
    lay <- layout_as_tree(g, flip.y=TRUE)

    #-----------------------------------------------------
    if ((sampling.design == "progenitor.intermediate")) {
      df <- lay
    }

    #-----------------------------------------------------
    if ((sampling.design == "progenitor.endpoint"))
    {
      #the generation difference as weights
      weight.scale <- c(1, edges$gendiff)

      #groups from user input
      grp <- c("A", edges$group)
      m <- cbind(lay, weight.scale, grp)
      m <- data.frame(m)

      #split groups as individual dataframes
      gps <- m %>% dplyr::group_by(grp)
      mygps <- dplyr::group_split(gps)

      #set gap-width for complex pedigrees. gap of 2 is good enough.
      gap = 2

      for (j in seq_along(mygps)) {
        if (j==1){
          #relevel the pedigree
          mygps[[j]]$relevel <- as.numeric(mygps[[j]]$V2)
        }
        else {
          #find minimum (last) level of previous group. Use it for re-scaling the pedigree
          mm <- min(as.numeric(mygps[[j-1]]$relevel))
          mygps[[j]]$relevel <- as.numeric(mygps[[j]]$V2)-gap-mm
        }
      }

      m <- bind_rows(mygps)
      m <- as.data.frame(m, stringsAsFactors = FALSE)
      lay[,2] <- m$relevel
      df <- lay
    }

    #-----------------------------------------------------
    if ((sampling.design == "sibling")) {

      #find root of the pedigree
      root <- which(degree(g, v=V(g), mode = "in")==0, useNames=TRUE)
      #find degree of the node
      d <- degree(g, root)

      mg <- merge(samples[2:length(samples$node),], edges, by.x="node",by.y="to", sort=FALSE)

      myg <- c("A", mg$group)
      m <- cbind(lay, samples, myg)
      m <- data.frame(m)

      #split groups as individual dataframes
      gps <- m %>% group_by(myg, meth)
      mygps <-group_split(gps)

      n = length(mygps)
      for (i in 1:n) {
        mygps[[i]]$V1 = i-0.5*(n+1)
      }

      mm <- bind_rows(mygps)
      mm  <- as.data.frame(mm, stringsAsFactors = FALSE)

      if (as.numeric(d) > 1){
        mm[1,1]=0
      }

      df <- merge(m, mm, by="node", sort=FALSE)
      df <- df[,c("V1.y", "V2.y")]
      colnames(df) <- c("V1","V2")
      df <- as.matrix(df)
    }

    #-----------------------------------------------------
    if ((sampling.design == "tree")) {

      #for trees with 1 stem. "Stem" column not needed.
      #---------------------------------------------------
      if (is.null(edges$Stem)) {
        m <- cbind(lay, samples)
        m <- data.frame(m)

        m <- m %>% dplyr::mutate(V1 = ifelse(meth=="N", 0, V1))

        i=1
        j <- length(which(m$meth=="Y"))
        k=0

        while (i <= (length(m$V1))){
          if ("Y" %in% m$meth[i])
          {
            #print(m[i,])
            m$V1[i]=j
            j=j-1
          } else {
            m$V2[i]=k
            k=k-1
          }
          i=i+1
        }
        #5 no. of meth == N
        nl <- length(which(m$meth=="N"))
        m <- m %>% mutate(V2 = ifelse(meth=="Y", -nl, V2))

        mm <- m[,c(1,2)]
        mm[,1] <- as.numeric(as.character(m$V1))
        mm[,2] <- as.numeric(as.character(m$V2))
        df  <- as.matrix(mm)
      }
      else {
        #for trees with 2 stem (rare scenario)
        #---------------------------------------------------
        tot.age <- as.numeric(samples[which(samples$Branchpoint_date==0),]$node)
        mygroup <- c(tot.age, edges$Stem)
        m <- cbind(lay, samples, mygroup)
        m <- data.frame(m)

        #for samples with methylation measurements, level branch to 0
        m <- m %>% mutate(V2 = ifelse(meth=="Y", 0, V2))

        mycount <- m %>% count(mygroup)

        i=0
        while (i < (length(m$V1)-1) ){
          i=i+2
          m$V1[i]=if (i <= mycount$n[1]) 0.5*(-i) else 0.5*(i-mycount$n[1]) # 0.5 for sample pairs
          m$V1[i+1]=if (i <= mycount$n[1]) 0.5*(-i) else 0.5*(i-mycount$n[1])
        }
        mm <- m[,c(1,2)]
        mm[,1] <- as.numeric(as.character(m$V1))
        mm[,2] <- as.numeric(as.character(m$V2))
        df  <- as.matrix(mm)
      }
    }

    #-----------------------------------------------------------

    V(g)$color <- ifelse(samples[V(g), 3] == "Y", "red", "gray")

    if (!is.null(out.pdf)) {
      pdf(paste(output.dir, out.pdf, ".pdf", sep=""), colormodel = 'cmyk', width = plot.width, height = plot.height)
    }

    pl <- plot.igraph(g, layout = df,
                      asp=aspect.ratio,
                      #edge.width = 1,
                      vertex.size = vertex.size,
                      vertex.frame.width = 0,
                      vertex.label = if (vertex.label==FALSE) NA,
                      vertex.label.cex = 0.6,
                      vertex.label.color = "black",
                      vertex.label.dist = 1,
                      vertex.label.degree = -pi/2,
                      edge.arrow.size = 0.1,
                      vertex.color = V(g)$color
    )
    print(pl)

    if (!is.null(out.pdf)) {
      dev.off()
    }
  }





DM.dataRead <- function(dFile, cytosine, posteriorMaxFilter) {
    command <- sprintf(paste0("grep --text -w ", cytosine , " %s"),
      paste(dFile, collapse = ""))

    if (.Platform$OS.type == "unix") {
      tmp <-
        fread(cmd = command,
          select = c("V1", "V2", "V3", "V4", "V7", "V8"))
      cnames <-
        c("seqnames",
          "start",
          "strand",
          "context",
          "posteriorMax",
          "status")
      colnames(tmp) <- cnames

    } else{
      tmp <- fread(
        dFile,
        skip = 0,
        sep = '\t',
        select = c(
          "seqnames",
          "start",
          "strand",
          "context",
          "posteriorMax",
          "status"
        ),
        showProgress = TRUE
      )

      tmp <- tmp %>% filter(tmp$context == cytosine)
    }

    tmp <-
      tmp %>% filter(tmp$posteriorMax >= posteriorMaxFilter &
          tmp$seqnames != "M" &  tmp$seqnames != "C")


    drops <- c("context", "posteriorMax")
    tmp <- tmp %>% select(!all_of(drops))
    #tmp <- tmp[ !(names(tmp) %in% drops)]
    # filter out based on context & posteriorMax

    return(tmp)
}


RC.dataRead <- function(dFile, cytosine, posteriorMaxFilter) {
    command <- sprintf(paste0("grep --text -w ", cytosine , " %s"),
      paste(dFile, collapse = ""))

    if (.Platform$OS.type == "unix") {
      tmp <- fread(cmd = command, select = c("V1", "V4", "V7", "V9"))
      cnames <-
        c("seqnames", "context", "posteriorMax", "rc.meth.lvl")
      colnames(tmp) <- cnames

    } else{
      tmp <- fread(
        dFile,
        skip = 0,
        sep = '\t',
        select = c("seqnames", "context", "posteriorMax", "rc.meth.lvl"),
        showProgress = TRUE
      )

      tmp <- tmp %>% filter(tmp$context == cytosine)
    }

    tmp <-
      tmp %>% filter(tmp$posteriorMax >= posteriorMaxFilter &
          tmp$seqnames != "M" &  tmp$seqnames != "C")

    # filter out based on context & posteriorMax
    drops <- c("context", "posteriorMax")
    #tmp <- tmp[,!(names(tmp) %in% drops)]
    tmp <- tmp %>% select(!all_of(drops))
    return(tmp)
}

# take 5 digit of decimal value posteriorMax column
floorDec<-function(rc.Mean ,x){
  y <- function(x, level=1) round(x - 5*10^(-level-1), level)
  res <-y(as.numeric(rc.Mean),x)
  return(res)
}

convertDMATRIX<-function(samples, edges, dmatrix)
{
    # reading matrix
    g <- graph_from_data_frame(edges, directed=TRUE, vertices=samples)

    # constracting tree
    samples.with.WGBS <- subset(samples, samples$meth == "Y")
    expand.grid.unique <- function(x, y, include.equals=FALSE)
    {
        x <- unique(x)
        y <- unique(y)
        k <- function(i)
        {
            z <- setdiff(y, x[seq_len(i-include.equals)])
            if(length(z)) cbind(x[i], z, deparse.level=0)
        }
        do.call(rbind, lapply(seq_along(x), k))
    }
    # making pedigres
    pairs.for.D <- expand.grid.unique(as.character(samples.with.WGBS$node), as.character(samples.with.WGBS$node))
    t1<-NULL
    t2<-NULL
    t0<-NULL
    s1<-NULL
    s2<-NULL
    for (a in 1:nrow(pairs.for.D))
    {
        # Loop the to and from over a distance matrix that matches D (only measured samples)
        path.this<-all_shortest_paths(g, from=pairs.for.D[a,1], to=pairs.for.D[a,2], mode="all")
        path.this<-names(unlist(path.this$res))
        s1[a]<-pairs.for.D[a,1]
        s2[a]<-pairs.for.D[a,2]
        t1[a]<-samples[which(samples$node == pairs.for.D[a,1]),2]$gen
        t2[a]<-samples[which(samples$node == pairs.for.D[a,2]),2]$gen
        samples.temp<-samples[which(as.character(samples$node) %in% path.this),2]
        t0[a]<-min(samples.temp)
    }
    # making data-frame
    tmpData <- data.frame(s1, s2, t0, t1, t2, stringsAsFactors=FALSE)
    colnames(tmpData)[3]<-"time0"
    colnames(tmpData)[4]<-"time1"
    colnames(tmpData)[5]<-"time2"
    # join Dmatrix with pedigree time
    tm1<-as.data.table(inner_join(dmatrix,tmpData ,  by = c("pair.1"="s1","pair.2"="s2")))
    tm2<-as.data.table(anti_join(dmatrix,tmpData ,  by = c("pair.1"="s1","pair.2"="s2")))
    tm2<-as.data.table(inner_join(tm2,tmpData ,  by = c("pair.1"="s2","pair.2"="s1")))
    dmatrixout<-as.data.table(rbind(tm1,tm2))
    # reordering table
    dmatrixout <- dmatrixout[, c(1, 2, 4, 5, 6, 3)]

    dmatrixout <- as.matrix(dmatrixout[,3:6])

    return(dmatrixout)
}

inputCheck <- function(...) {
  # getting var parameters
  var <- list(...)
  # input file checking

  # check if genTable is exist then check file exist one-by-one
  if (!file.exists(var[[1]])){
    stop("Sample file dosen't exist.")
  }
  gen_tbl <- fread(var[[1]])
  gen_tbl <- gen_tbl %>% filter(gen_tbl$meth=="Y")
  for (i in seq_len(NROW(gen_tbl))){
    if (!file.exists(gen_tbl[[i,1]])){
      stop(paste0("File ",gen_tbl[[i,1]], " Not Exist!"))
    }
  }

  # cytosines part checking
  list_cyt <- c("CHH", "CHG", "CG")
  cytosines <- var[[2]]
  if (!cytosines %in% list_cyt | length(cytosines) > 1) {
    stop("Please enter the valid Cytosines. (CHH/CHG/CG)")
  }

  # posteriorMaxFiltervalue part checking
  pMFilter <- as.numeric(var[[3]])
  if (pMFilter <= 0.80 | pMFilter > 1) {
    stop("posteriorMax value must be < 1 and >= 0.80 ")
  }
  return(0)
}

statusStringCheck <-  function(file_A, file_B){
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

  strTocheckFileB <- utils::head(file_B$status[1])
  if (strTocheckFileB %in% list_status) {
    # replace pattern in Data-set B
    file_B$status <- str_replace_all(file_B$status,
                      pattern = "Unmethylated", replacement = "U")
    file_B$status <- str_replace_all(file_B$status,
                      pattern = "Intermediate", replacement = "I")
    file_B$status <- str_replace_all(file_B$status,
                      pattern = "Methylated", replacement = "M")
  }
  return(list(file_A, file_B))

}

getNames <- function(nameDF,genTable){
  #genTable <- fread(genTable)
  strSearch<-nameDF
  tmp_A <- genTable[genTable$filename == strSearch,][,2]
  #tmp_B <- genTable[genTable$filename == strSearch,][,3]
  #if (tmp_B == "" | is.na(tmp_B) | is.null(tmp_B) ){
  #  tmp_B<-""
  #  gen_name <- paste0(tmp_A,tmp_B)
  #}else{
  #  gen_name <- paste0(tmp_A,"-",tmp_B)
  #}
  #return(gen_name)
  return(tmp_A)
}


buildPedigree <- function(nodelist, edgelist , cytosine="CG", posteriorMaxFilter=0.99) {
    mt <- startTime("constracting pedigree ...\n")
    # constracting D-Matrices from sample file
    dmatrix <- dMatrix(nodelist, cytosine, posteriorMaxFilter)

    # constracting Rc.Meth.lvl from sample file
    message("generating Rc.Meth.lvl from sample file....\n")
    rclvl <- rc.meth.lvl(nodelist, cytosine, posteriorMaxFilter)
    props <- rclvl[which(as.character(rclvl[, 2]) == cytosine), ]
    outliers <- "none"
    props <- rclvl[which(!is.element(rclvl[, 1], outliers) == TRUE), ]
    tmpP0uu<- 1 - mean(as.numeric(as.character(props[, 3])))
    message("finilizing pedegree data...")
    # converting D-matrix into pedegree

    # reading matrix
    edges   <- fread(edgelist, header=TRUE, skip=0)
    samples <- fread(nodelist, header=TRUE, skip=0, select = c(2,3,4))
    # check file 'sample' if there is column 'Branchpoint_date' then it's tree
    if (!is.null(samples$Branchpoint_date)){
        #rename column to 'gen' to run normal covertMatrix fun.
        colnames(samples)[2]<-"gen"
    }

    tmpPedegree <- convertDMATRIX(samples, edges, dmatrix)

    cat(stopTime(mt))
    return(list(Pdata=tmpPedegree,tmpp0=tmpP0uu))
}

FtestRSS<-function(pedigree.select, pedigree.null)
{

    # Reading data
    est<-dget(pedigree.select)
    estN<-dget(pedigree.null)

    if (estN$model == "ABnull.R")
    {
          # Testing
          RSSf<-est$estimates[1,"value"]
          RSSr<-sum((estN$pedigree[,"residual"])^2)
          Npara_r<-1
          Npara_f<-5
          dfF<-length(est$pedigree[,"residual"])-5
          dfR<-length(estN$pedigree[,"residual"])-1
          dfN<-dfR-dfF
          Fvalue<-((RSSr - RSSf)/(Npara_f-Npara_r))/(RSSf/dfF)
          pvalue<-pf(Fvalue, dfN, dfF, lower.tail=FALSE)
          output<-c(RSSf, RSSr, dfF, dfR, Fvalue, pvalue)
          names(output)<-c("RSS_F", "RSS_R", "df_F", "df_R", "Fvalue", "pvalue")
          outfinal<-list(output, est$estimates, estN$estimates)
          names(outfinal)<-c("Ftest", "est.selection", "est.neutral")
    }

    if (estN$model != "ABnull.R")
    {
      # Testing
      RSSf<-est$estimates[1,"value"]
      RSSr<-estN$estimates[1,"value"]
      Npara_r<-4
      Npara_f<-5
      dfF<-length(est$pedigree[,"residual"])-5
      dfR<-length(estN$pedigree[,"residual"])-4
      dfN<-dfR-dfF
      Fvalue<-((RSSr - RSSf)/(Npara_f-Npara_r))/(RSSf/dfF)
      pvalue<-pf(Fvalue, dfN, dfF, lower.tail=FALSE)
      output<-c(RSSf, RSSr, dfF, dfR, Fvalue, pvalue)
      names(output)<-c("RSS_F", "RSS_R", "df_F", "df_R", "Fvalue", "pvalue")
      outfinal<-list(output, est$estimates, estN$estimates)
      names(outfinal)<-c("Ftest", "est.selection", "est.neutral")

    }

    outfinal
}


dMatrix <- function(nodelist, cytosine, posteriorMaxFilter) {
    # checking errors
    inputCheck(nodelist, cytosine, posteriorMaxFilter)
    genTable <- fread(nodelist)
    #---------------------------Filter based on meth
    genTable <- genTable %>% filter(genTable$meth=="Y")
    #---------------------------
    pairs <- combn(genTable$filename, 2)
    final_ds <- runMatrix(pairs, cytosine, posteriorMaxFilter, genTable)

    final_ds<-final_ds[mixedorder(final_ds$X1),]
    colnames(final_ds)<-(c("pair.1", "pair.2", "D.value"))

    #dmatrix <- dMsaveResult(final_ds, cytosine, posteriorMaxFilter)
    message("generating d-matrics done.\n")
    rm(genTable)
    return(final_ds)
}


runMatrix <- function(pairs, cytosine, posteriorMaxFilter,genTable){
     flag=TRUE
     name_ds<-NULL
     pair_len <- length(pairs)/2
       for (i in seq_len(pair_len)){
        df <-pairs[,i]
        name_ds[1] <- getNames(df[1], genTable)
        name_ds[2] <- getNames(df[2], genTable)

        cat(paste0("Reading sample: ",name_ds[[1]], " and ",
         name_ds[[2]], " ( ", i , " out of ",length(pairs)/2," pairs )"),"\n")

        cytosine<-as.character(cytosine)
        file_A <-DM.dataRead(df[1],cytosine, posteriorMaxFilter)
        file_B <-DM.dataRead(df[2],cytosine, posteriorMaxFilter)
        # check ad replace pattern in Data-set A
        returnFile <- statusStringCheck(file_A, file_B)
        file_A <- returnFile[[1]]
        file_B <- returnFile[[2]]
        cat("Computing divergence matrix...\n")
        file_A$seqnames<-as.character(file_A$seqnames)
        file_B$seqnames<-as.character(file_B$seqnames)

        tmp_db <- as.data.table(
                  inner_join(file_A, file_B, by = c("seqnames","start","strand")))
        #set status 0=rows is same, 1=M/U 2=I
        rm(file_A,file_B,returnFile)

        tmp_db$state <- ifelse(tmp_db$status.x==tmp_db$status.y,0,
                        (ifelse((tmp_db$status.x=="I" | tmp_db$status.y=="I" ),2,1)))

        # substract number of Intermediate in data-set
        number_none_inter <- sum(tmp_db$state==1)     #M --> U OR U --> M
        number_intermediate <- sum(tmp_db$state==2)   #Intermediate
        Total <- NROW(tmp_db$state)                   #Total rows
        # number of (M --> U or U --> M  and 1/2 I ) / total of rows
        D <- (number_none_inter+(number_intermediate*0.5))/Total
        # reformat
        D <-floorDec(as.numeric(D),5)
          if (flag == TRUE){
             tmp_big<- data.frame(matrix(ncol = 3, nrow = 1))
             tmp_big$X1<-name_ds[[1]]
             tmp_big$X2<-name_ds[[2]]
             tmp_big$X3<-D
             flag = FALSE
          } else {
          tmp<-NULL
          tmp<-list(name_ds[[1]],name_ds[[2]],D)
          tmp_big<-rbind(tmp_big,tmp)
          rm(tmp_db)
          }
         cat(paste0(name_ds[[1]], " and ",name_ds[[2]], " is done! \n"))
         cat("|--------------------------------------------------|\n")
         name_ds<-NULL
       }
    return(tmp_big)
}






BOOTmodel<-function(pedigree.data, Nboot, out.dir, out.name)
{
  pedigree.data<-dget(pedigree.data)
  model<-pedigree.data$model

#################################
########### ABselectMM ##########
#################################
  if (model == "ABselectMM.R")
  {

      ## Reading the dataset for bootstrapping and extracting the parameter settings
      settings<-pedigree.data$settings
      est<-pedigree.data$estimates
      eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
      eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
      optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
      #p0mm<-round(as.numeric(as.character(settings[which(settings[,1] == "p0mm"),2])),16)
      p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
      #p0um<-round(as.numeric(as.character(settings[which(settings[,1] == "p0um"),2])),15)
      pedigree<-pedigree.data$pedigree
      p0mm = 1-p0uu
      p0um = 0


    	if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

    	## Defining the divergence function
    	divergence <- function(pedigree, p0mm, p0um, p0uu, param)
    	{

    	  ## Initializing parameters
    	  PrMM <- p0mm
    	  PrUM <- p0um
    	  PrUU <- p0uu
    	  alpha <- param[1]
    	  bet   <- param[2]
    	  weight <-param[3]
    	  sel    <-param[4]

    	  ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    	  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

    	  element11<-((1-alpha)^2)*sel
    	  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
    	  element13<-(alpha^2)
    	  rowtotal1<-element11 + element12 + element13

    	  element21<-(1/4*(bet + 1 - alpha)^2)*sel
    	  element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
    	  element23<-(1/4*(alpha + (1 - bet))^2)
    	  rowtotal2<-element21 + element22 + element23

    	  element31<-(bet^2)*sel
    	  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
    	  element33<-((1-bet))^2
    	  rowtotal3<-element31 + element32 + element33

    	  ## Defining the generation (or transition) matrix
    	  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
    	                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
    	                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

    	  ## Calculating the expected divergence for every observed pair in 'pedigree.txt'
    	  Dt1t2<-NULL

    	  for (p in seq_len(NROW(pedigree)))
    	  {

    	    ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
    	    svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
    	    svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    	    svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    	    svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    	    svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    	    svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    	    svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

    	    ## Conditional divergences
    	    dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
    	                        svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

    	    dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
    	                        svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

    	    dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
    	                        svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

    	    ## Total (weighted) divergence
    	    Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


    	  }

    	  ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
    	  puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
    	  puuinf.est<- puuinf.est[1,1]
    	  divout<-list(puuinf.est, Dt1t2)

    	  return(divout)

    	}

    	## Defining the Least Square function to be minimized
    	LSE_intercept<-function(param_int)
    	{
    	  sum((pedigree[,4] - param_int[5] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[2]])^2) +
    	    eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[1]]- eqp)^2)
    	}


      ## Initializing
    	final<-NULL
    	counter<-0


    	## Defining starting values
    	alpha.start  <-est[1,1]
    	beta.start   <-est[1,2]
    	weight.start <-est[1,3]
    	sel.start<-est[1,4]
    	intercept.start <-est[1,5]
    	param_int0 = c(alpha.start, beta.start, weight.start, sel.start, intercept.start)


    	    ## Start of boostrap loops
        	for (booting in seq_len(Nboot))
        	{

        	  opt.out<-NULL
        		pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

        			counter<-counter+1

        			message("Bootstrap interation: ", counter/Nboot, "\n")

        			opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
        			alphafinal<-as.numeric(opt.out[1])
        			betfinal<-as.numeric(opt.out[2])
        			weightfinal<-as.numeric(opt.out[3])
        			selfinal<-as.numeric(opt.out[4])

        			## Calculating equilibrium frequencies based on the model estimates
        			svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

        			element11<-((1-alphafinal)^2)*selfinal
        			element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
        			element13<-(alphafinal^2)
        			rowtotal1<-element11 + element12 + element13

        			element21<-(1/4*(betfinal + 1 - alphafinal)^2)*selfinal
        			element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
        			element23<-(1/4*(alphafinal + (1 - betfinal))^2)
        			rowtotal2<-element21 + element22 + element23

        			element31<-(betfinal^2)*selfinal
        			element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
        			element33<-((1-betfinal))^2
        			rowtotal3<-element31 + element32 + element33

        			## Defining the generation (or transition) matrix
        			Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
        			                      element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
        			                      element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

        			## Note: This is an approximation to the equilibrium values
        			pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
        			PrMMinf<- pinf.vec[1,3]
        			PrUMinf<- pinf.vec[1,2]
        			PrUUinf<- pinf.vec[1,1]
        			opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start,
        			                sel.start, intercept.start)

        			## Collecting the results
        			final<-rbind(final, opt.out)


        	} # End of Bootstrap loops


    	 colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
    	 colnames(final)[14:16]<-c("PrMMinf", "PrUMinf", "PrUUinf")

    	 SE.alpha<-sd(final[,1],na.rm=TRUE)
    	 SE.beta<-sd(final[,2],na.rm=TRUE)
    	 SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
    	 SE.weight<-sd(final[,3],na.rm=TRUE)
    	 SE.sel.coef<-sd(final[,4], na.rm=TRUE)
    	 SE.inter<-sd(final[,5],na.rm=TRUE)
    	 SE.PrMMinf<-sd(final[,14],na.rm=TRUE)
    	 SE.PrUMinf<-sd(final[,15],na.rm=TRUE)
    	 SE.PrUUinf<-sd(final[,16],na.rm=TRUE)

    	 CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
    	 CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
    	 CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
    	 CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
    	 CI.sel.coef<-quantile(final[,4],probs=c(0.025, 0.975))
    	 CI.inter<-quantile(final[,5],probs=c(0.025, 0.975))
    	 CI.PrMMinf<-quantile(final[,14],probs=c(0.025, 0.975))
    	 CI.PrUMinf<-quantile(final[,15],probs=c(0.025, 0.975))
    	 CI.PrUUinf<-quantile(final[,16],probs=c(0.025, 0.975))

    	 SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.sel.coef, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
    	 CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.sel.coef, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

    	 SE.out<-cbind(SE, CI)
    	 colnames(SE.out)[1]<-"SE"
    	 rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "sel.coef", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    	 final<-data.frame(final)
    	 good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

    	 SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
    	 names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

    	 ## Ouputting result datasets
    	 dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

    	 return(SE.out)

  } # End of "ABselectMM" if statement


  #################################
  ########### ABselectUU ##########
  #################################
  if (model == "ABselectUU.R")
  {

    ## Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0


    if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

    ## Defining the divergence function
    divergence <- function(pedigree, p0mm, p0um, p0uu, param)
    {

      ## Initializing parameters
      PrMM <- p0mm
      PrUM <- p0um
      PrUU <- p0uu
      alpha <- param[1]
      bet   <- param[2]
      weight <-param[3]
      sel    <-param[4]

      ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
      svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

      element11<-((1-alpha)^2)
      element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
      element13<-(alpha^2)*sel
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(bet + 1 - alpha)^2)
      element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
      element23<-(1/4*(alpha + (1 - bet))^2)*sel
      rowtotal2<-element21 + element22 + element23

      element31<-(bet^2)
      element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
      element33<-(((1-bet))^2)*sel
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

      ## Calculating the expected divergence for every observed pair in 'pedigree.txt'
      Dt1t2<-NULL

      for (p in seq_len(NROW(pedigree)))
      {

        ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
        svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
        svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
        svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
        svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

        ## Conditional divergences
        dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                            svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

        dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                            svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

        dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                            svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

        ## Total (weighted) divergence
        Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


      }

      ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
      puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      puuinf.est<- puuinf.est[1,1]
      divout<-list(puuinf.est, Dt1t2)

      return(divout)

    }

    ## Defining the Least Square function to be minimized
    LSE_intercept<-function(param_int)
    {
      sum((pedigree[,4] - param_int[5] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[2]])^2) +
        eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[1]]- eqp)^2)
    }



    ## Initializing
    final<-NULL
    counter<-0


    ## Defining starting values
    alpha.start  <-est[1,1]
    beta.start   <-est[1,2]
    weight.start <-est[1,3]
    sel.start<-est[1,4]
    intercept.start <-est[1,5]
    param_int0 = c(alpha.start, beta.start, weight.start, sel.start, intercept.start)


    ## Start of boostrap loops
    for (booting in seq_len(Nboot))
    {

      opt.out<-NULL
      pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

      counter<-counter+1

      message("Bootstrap interation: ", counter/Nboot, "\n")

      opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
      alphafinal<-as.numeric(opt.out[1])
      betfinal<-as.numeric(opt.out[2])
      weightfinal<-as.numeric(opt.out[3])
      selfinal<-as.numeric(opt.out[4])

      ## Calculating equilibrium frequencies based on the model estimates
      svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

      element11<-((1-alphafinal)^2)
      element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
      element13<-(alphafinal^2)*selfinal
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(betfinal + 1 - alphafinal)^2)
      element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
      element23<-(1/4*(alphafinal + (1 - betfinal))^2)*selfinal
      rowtotal2<-element21 + element22 + element23

      element31<-(betfinal^2)
      element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
      element33<-(((1-betfinal))^2)*selfinal
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

      ## Note: This is an approximation to the equilibrium values
      pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      PrMMinf<- pinf.vec[1,3]
      PrUMinf<- pinf.vec[1,2]
      PrUUinf<- pinf.vec[1,1]
      opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start,
                      sel.start, intercept.start)
      final<-rbind(final, opt.out)


    } # End of Bootstrap loops


    colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
    colnames(final)[14:16]<-c("PrMMinf", "PrUMinf", "PrUUinf")

    SE.alpha<-sd(final[,1],na.rm=TRUE)
    SE.beta<-sd(final[,2],na.rm=TRUE)
    SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
    SE.weight<-sd(final[,3],na.rm=TRUE)
    SE.sel.coef<-sd(final[,4], na.rm=TRUE)
    SE.inter<-sd(final[,5],na.rm=TRUE)
    SE.PrMMinf<-sd(final[,14],na.rm=TRUE)
    SE.PrUMinf<-sd(final[,15],na.rm=TRUE)
    SE.PrUUinf<-sd(final[,16],na.rm=TRUE)

    CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
    CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
    CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
    CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
    CI.sel.coef<-quantile(final[,4],probs=c(0.025, 0.975))
    CI.inter<-quantile(final[,5],probs=c(0.025, 0.975))
    CI.PrMMinf<-quantile(final[,14],probs=c(0.025, 0.975))
    CI.PrUMinf<-quantile(final[,15],probs=c(0.025, 0.975))
    CI.PrUUinf<-quantile(final[,16],probs=c(0.025, 0.975))

    SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.sel.coef, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
    CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.sel.coef, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

    SE.out<-cbind(SE, CI)
    colnames(SE.out)[1]<-"SE"
    rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "sel.coef", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    final<-data.frame(final)
    good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

    SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
    names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

    ## Ouputting result datasets
    dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

    return(SE.out)

  } # End of "ABselectUU" if statement





  ###################################
  ########### ABfixselectMM #########
  ###################################
  if (model == "ABfixselectMM.R")
  {

    ## Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0


    if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

    ##### Defining the divergence function
    divergence <- function(pedigree, sel, p0mm, p0um, p0uu, param)
    {

      ## Initializing parameters
      PrMM <- p0mm
      PrUM <- p0um
      PrUU <- p0uu
      sel  <- sel
      alpha <- param[1]
      bet   <- param[2]
      weight <-param[3]


      ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
      svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

      element11<-((1-alpha)^2)*sel
      element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
      element13<-(alpha^2)
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(bet + 1 - alpha)^2)*sel
      element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
      element23<-(1/4*(alpha + (1 - bet))^2)
      rowtotal2<-element21 + element22 + element23

      element31<-(bet^2)*sel
      element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
      element33<-((1-bet))^2
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)


      ## Calculating the expected divergence for every observed pair in 'pedigree.txt'
      Dt1t2<-NULL

      for (p in seq_len(NROW(pedigree)))
      {

        ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
        svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
        svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
        svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
        svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

        ## Conditional divergences
        dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                            svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

        dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                            svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

        dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                            svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

        ## Total (weighted) divergence
        Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


      }

      ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
      puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      puuinf.est<- puuinf.est[1,1]
      divout<-list(puuinf.est, Dt1t2)

      return(divout)

    }

    ###### Defining the Least Square function to be minimized
    ###### Note the equilibrium constraint, which can be made as small as desired.

    LSE_intercept<-function(param_int)
    {
      sum((pedigree[,4] - param_int[4] - divergence(pedigree, sel, p0mm, p0um, p0uu, param_int[1:3])[[2]])^2) +
        eqp.weight*nrow(pedigree)*((divergence(pedigree, sel, p0mm, p0um, p0uu, param_int[1:3])[[1]]- eqp)^2)
    }



    ## Initializing
    final<-NULL
    counter<-0


    ## Defining starting values
    alpha.start  <-est[1,1]
    beta.start   <-est[1,2]
    weight.start <-est[1,3]
    intercept.start <-est[1,4]
    param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)


    ## Start of boostrap loops
    for (booting in NROW(Nboot))
    {

      opt.out<-NULL
      pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

      counter<-counter+1

      message("Bootstrap interation: ", counter/Nboot, "\n")

      opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
      alphafinal<-as.numeric(opt.out[1])
      betfinal<-as.numeric(opt.out[2])
      weightfinal<-as.numeric(opt.out[3])
      selfinal<-sel


      ## Here we want to calculate the equilibrium freqs based on the estimates
      svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

      element11<-((1-alphafinal)^2)*selfinal
      element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
      element13<-(alphafinal^2)
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(betfinal + 1 - alphafinal)^2)*selfinal
      element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
      element23<-(1/4*(alphafinal + (1 - betfinal))^2)
      rowtotal2<-element21 + element22 + element23

      element31<-(betfinal^2)*selfinal
      element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
      element33<-((1-betfinal))^2
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)


      ## Note: This is an approximation to the equilibrium values
      pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      PrMMinf<- pinf.vec[1,3]
      PrUMinf<- pinf.vec[1,2]
      PrUUinf<- pinf.vec[1,1]
      opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
      final<-rbind(final, opt.out)


    } # End of Bootstrap loops


    colnames(final)[1:4]<-c("alpha", "beta", "weight", "intercept")
    colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")

    SE.alpha<-sd(final[,1],na.rm=TRUE)
    SE.beta<-sd(final[,2],na.rm=TRUE)
    SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
    SE.weight<-sd(final[,3],na.rm=TRUE)
    SE.inter<-sd(final[,4],na.rm=TRUE)
    SE.PrMMinf<-sd(final[,13],na.rm=TRUE)
    SE.PrUMinf<-sd(final[,14],na.rm=TRUE)
    SE.PrUUinf<-sd(final[,15],na.rm=TRUE)

    CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
    CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
    CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
    CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
    CI.inter<-quantile(final[,4],probs=c(0.025, 0.975))
    CI.PrMMinf<-quantile(final[,13],probs=c(0.025, 0.975))
    CI.PrUMinf<-quantile(final[,14],probs=c(0.025, 0.975))
    CI.PrUUinf<-quantile(final[,15],probs=c(0.025, 0.975))

    SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
    CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

    SE.out<-cbind(SE, CI)
    colnames(SE.out)[1]<-"SE"
    rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    final<-data.frame(final)
    good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

    SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
    names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

    ## Ouputting result datasets
    dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

    return(SE.out)

  } # End of "ABfixselectMM" if statement






  ###################################
  ########### ABfixselectUU #########
  ###################################
  if (model == "ABfixselectUU.R")
  {

    ## Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0


    if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

    ##### Defining the divergence function
    divergence <- function(pedigree, sel, p0mm, p0um, p0uu, param)
    {

      ## Initializing parameters
      PrMM <- p0mm
      PrUM <- p0um
      PrUU <- p0uu
      sel  <- sel
      alpha <- param[1]
      bet   <- param[2]
      weight <-param[3]


    ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

    element11<-((1-alpha)^2)
    element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
    element13<-(alpha^2)*sel
    rowtotal1<-element11 + element12 + element13

    element21<-(1/4*(bet + 1 - alpha)^2)
    element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
    element23<-(1/4*(alpha + (1 - bet))^2)*sel
    rowtotal2<-element21 + element22 + element23

    element31<-(bet^2)
    element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
    element33<-(((1-bet))^2)*sel
    rowtotal3<-element31 + element32 + element33

    ## Defining the generation (or transition) matrix
    Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                          element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                          element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

    ## Calculating the expected divergence for every observed pair in 'pedigree.txt'
    Dt1t2<-NULL

    for (p in seq_len(NROW(pedigree)))
    {

      ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
      svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
      svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

      ## Conditional divergences
      dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                          svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

      dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                          svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

      dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                          svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

      ## Total (weighted) divergence
      Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


    }

    ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
    puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
    puuinf.est<- puuinf.est[1,1]
    divout<-list(puuinf.est, Dt1t2)

    return(divout)

  }

    ###### Defining the Least Square function to be minimized
    ###### Note the equilibrium constraint, which can be made as small as desired.

    LSE_intercept<-function(param_int)
    {
      sum((pedigree[,4] - param_int[4] - divergence(pedigree, sel, p0mm, p0um, p0uu, param_int[1:3])[[2]])^2) +
        eqp.weight*nrow(pedigree)*((divergence(pedigree, sel, p0mm, p0um, p0uu, param_int[1:3])[[1]]- eqp)^2)
    }



    ## Initializing
    final<-NULL
    counter<-0


    ## Defining starting values
    alpha.start  <-est[1,1]
    beta.start   <-est[1,2]
    weight.start <-est[1,3]
    intercept.start <-est[1,4]
    param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)


    ## Start of boostrap loops
    for (booting in seq_len(Nboot))
    {

      opt.out<-NULL
      pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

      counter<-counter+1

      message("Bootstrap interation: ", counter/Nboot, "\n")

      opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
      alphafinal<-as.numeric(opt.out[1])
      betfinal<-as.numeric(opt.out[2])
      weightfinal<-as.numeric(opt.out[3])
      selfinal<-sel


      ## Here we want to calculate the equilibrium freqs based on the estimates
      svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

      element11<-((1-alphafinal)^2)
      element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
      element13<-(alphafinal^2)*selfinal
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(betfinal + 1 - alphafinal)^2)
      element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
      element23<-(1/4*(alphafinal + (1 - betfinal))^2)*selfinal
      rowtotal2<-element21 + element22 + element23

      element31<-(betfinal^2)
      element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
      element33<-(((1-betfinal))^2)*selfinal
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)


      ## Note: This is an approximation to the equilibrium values
      pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      PrMMinf<- pinf.vec[1,3]
      PrUMinf<- pinf.vec[1,2]
      PrUUinf<- pinf.vec[1,1]
      opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
      final<-rbind(final, opt.out)



    } # End of Bootstrap loops


    colnames(final)[1:4]<-c("alpha", "beta", "weight", "intercept")
    colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")

    SE.alpha<-sd(final[,1],na.rm=TRUE)
    SE.beta<-sd(final[,2],na.rm=TRUE)
    SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
    SE.weight<-sd(final[,3],na.rm=TRUE)
    SE.inter<-sd(final[,4],na.rm=TRUE)
    SE.PrMMinf<-sd(final[,13],na.rm=TRUE)
    SE.PrUMinf<-sd(final[,14],na.rm=TRUE)
    SE.PrUUinf<-sd(final[,15],na.rm=TRUE)

    CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
    CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
    CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
    CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
    CI.inter<-quantile(final[,4],probs=c(0.025, 0.975))
    CI.PrMMinf<-quantile(final[,13],probs=c(0.025, 0.975))
    CI.PrUMinf<-quantile(final[,14],probs=c(0.025, 0.975))
    CI.PrUUinf<-quantile(final[,15],probs=c(0.025, 0.975))

    SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
    CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

    SE.out<-cbind(SE, CI)
    colnames(SE.out)[1]<-"SE"
    rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    final<-data.frame(final)
    good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

    SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
    names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

    ## Ouputting result datasets
    dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

    return(SE.out)

  } # End of "ABfixselectUU" if statement







  ###################################
  ########### ABneutral      ########
  ###################################

  if (model == "ABneutral.R")
  {
  ## Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0


  if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

  ##### Defining the divergence function
  divergence <- function(pedigree, p0mm, p0um, p0uu, param)
  {

    ## Initializing parameters
    PrMM <- p0mm
    PrUM <- p0um
    PrUU <- p0uu
    alpha <- param[1]
    bet <- param[2]
    weight <- param[3]


    ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)



    ## Defining the generation (or transition) matrix
    Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha, alpha^2,
                          1/4*(bet + 1 - alpha)^2, 1/2*(bet + 1 - alpha)*(alpha + 1 - bet), 1/4*(alpha + 1 - bet)^2,
                          bet^2, 2*(1-bet)*bet, (1-bet)^2), nrow=3, byrow=TRUE)


    ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
    Dt1t2<-NULL

    for (p in seq_len(NROW(pedigree)))
    {

      ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
      svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
      svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

      ## Conditional divergences
      dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                          svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

      dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                          svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

      dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                          svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

      ## Total (weighted) divergence
      Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


    }

    # Pr(UU) at equilibrium given alpha and beta
    puuinf.est<- (bet* ((1-bet)^2 - (1-alpha)^2 -1))/((alpha + bet)*((alpha + bet -1)^2 - 2))
    divout<-list(puuinf.est, Dt1t2)

    return(divout)

  }

  ###### Defining the Least Square function to be minimized
  ###### Note the equilibrium constraint, which can be made as small as desired.

  LSE_intercept<-function(param_int)
  {
    sum((pedigree[,4] - param_int[4] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[2]])^2) +
      eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[1]]- eqp)^2)
  }



  ## Initializing
  final<-NULL
  counter<-0


  ## Defining starting values
  alpha.start  <-est[1,1]
  beta.start   <-est[1,2]
  weight.start <-est[1,3]
  intercept.start <-est[1,4]
  param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)


  ## Start of boostrap loops
  for (booting in seq_len(Nboot))
  {

    opt.out<-NULL
    pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

    counter<-counter+1

    message("Bootstrap interation: ", counter/Nboot, "\n")


    opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
    alphafinal<-opt.out[1]
    betfinal<-opt.out[2]
    PrMMinf <- (alphafinal* ((1-alphafinal)^2 - (1-betfinal)^2 -1))/((alphafinal + betfinal)*((alphafinal + betfinal -1)^2 - 2))
    PrUMinf <- (4*alphafinal*betfinal*(alphafinal + betfinal -2))/((alphafinal + betfinal)*((alphafinal + betfinal -1)^2 -2))
    PrUUinf <- (betfinal* ((1-betfinal)^2 - (1-alphafinal)^2 -1))/((alphafinal + betfinal)*((alphafinal + betfinal -1)^2 - 2))
    opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
    final<-rbind(final, opt.out)


  } # End of Bootstrap loops


  colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
  colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")

  SE.alpha<-sd(final[,1],na.rm=TRUE)
  SE.beta<-sd(final[,2],na.rm=TRUE)
  SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
  SE.weight<-sd(final[,3],na.rm=TRUE)
  SE.inter<-sd(final[,4],na.rm=TRUE)
  SE.PrMMinf<-sd(final[,13],na.rm=TRUE)
  SE.PrUMinf<-sd(final[,14],na.rm=TRUE)
  SE.PrUUinf<-sd(final[,15],na.rm=TRUE)

  CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
  CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
  CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
  CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
  CI.inter<-quantile(final[,4],probs=c(0.025, 0.975))
  CI.PrMMinf<-quantile(final[,13],probs=c(0.025, 0.975))
  CI.PrUMinf<-quantile(final[,14],probs=c(0.025, 0.975))
  CI.PrUUinf<-quantile(final[,15],probs=c(0.025, 0.975))

  SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
  CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

  SE.out<-cbind(SE, CI)
  colnames(SE.out)[1]<-"SE"
  rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

  final<-data.frame(final)
  good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

  SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
  names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

  ## Ouputting result datasets
  dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

  return(SE.out)

} # End of "ABneutral if statement







  ###################################
  ########### ABneutralSOMA     #####
  ###################################

 if (model == "ABneutralSOMA.R")
  {

allow.neg.intercept="yes"


  ##### Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0
  #props<-pedigree.data$prop.data.used


  if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}





##### Defining the divergence function
  divergence <- function(pedigree, p0mm, p0um, p0uu, param)
  {

    ## Initializing parameters
    PrMM <- p0mm
    PrUM <- p0um
    PrUU <- p0uu
    alpha <- param[1]
    bet <- param[2]
    weight <- param[3]


## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)



  ## Defining the generation (or transition) matrix for the mitotic case
    Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha,alpha^2,
                          bet*(1-alpha), (1-alpha)*(1-bet)+alpha*bet, alpha*(1-bet),
              bet^2, 2*(1-bet)*bet, (1-bet)^2),nrow=3, byrow=TRUE)


  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
    Dt1t2<-NULL

      for (p in seq_len(NROW(pedigree)) )
      {

      ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
      svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
      svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

      ## Conditional divergences
      dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

      dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

      dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

      ## Total (weighted) divergence
      Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


      }

    # Pr(MM) at equilibrium given alpha and beta
    puuinf.est<-(bet^2)/((alpha+bet)^2)
    divout<-list(puuinf.est, Dt1t2)

    return(divout)

  }


  ###### Defining the Least Square function to be minimized
  ###### Note the equilibrium constraint, which can be made as small as desired.

  LSE_intercept<-function(param_int)
  {
    sum((pedigree[,4] - param_int[4] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[2]])^2) +
      eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[1]]- eqp)^2)
  }






##### Initializing
  final<-NULL
  counter<-0
  opt.out<-NULL

  ## Defining starting values
  alpha.start  <-est[1,1]
  beta.start   <-est[1,2]
  weight.start <-est[1,3]
  intercept.start <-est[1,4]
  param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)


  for (booting in seq_len(Nboot))
  {
    pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

      ## Initializing
      counter<-counter+1

      message("Bootstrap interation: ", counter/Nboot, "\n")

            opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
            alphafinal<-opt.out[1]
            betfinal<-opt.out[2]
            PrMMinf <- (alphafinal^2)/((alphafinal+betfinal)^2)
            PrUMinf <- (2*alphafinal*betfinal)/((alphafinal+betfinal)^2)
            PrUUinf <- (betfinal^2)/((alphafinal+betfinal)^2)
            opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
            final<-rbind(final, opt.out)


  } # End of Bootstrap loops


   colnames(final)[1:4]<-c("alpha", "beta", "weight", "intercept")
   colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")

   if (allow.neg.intercept == "yes")
   { index<-which(final["alpha"] > 0 & final["beta"] > 0 & final["convcode"] == 0)}

   if (allow.neg.intercept == "no")
   {index<-which(final["alpha"] > 0 & final["beta"] > 0 & final["intercept"] > 0 & final["convcode"] == 0)}

   final<-final[index,]

   SE.alpha<-sd(final[,1],na.rm=TRUE)
   SE.beta<-sd(final[,2],na.rm=TRUE)
   SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
   SE.weight<-sd(final[,3],na.rm=TRUE)
   SE.inter<-sd(final[,4],na.rm=TRUE)
   SE.PrMMinf<-sd(final[,13],na.rm=TRUE)
   SE.PrUMinf<-sd(final[,14],na.rm=TRUE)
   SE.PrUUinf<-sd(final[,15],na.rm=TRUE)

   CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
   CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
   CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
   CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
   CI.inter<-quantile(final[,4],probs=c(0.025, 0.975))
   CI.PrMMinf<-quantile(final[,13],probs=c(0.025, 0.975))
   CI.PrUMinf<-quantile(final[,14],probs=c(0.025, 0.975))
   CI.PrUUinf<-quantile(final[,15],probs=c(0.025, 0.975))

   SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
   CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

   SE.out<-cbind(SE, CI)
   colnames(SE.out)[1]<-"SE"
   rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    final<-data.frame(final)
      good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

      SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
      names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

      ## Ouputting result datasets
      dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

  return(SE.out)

  }

} #End of function










ABselectUUSOMA<-function(pedigree.data, p0uu, eqp, eqp.weight, Nstarts, out.dir, out.name)
{


 allow.neg.intercept="no"

##### Defining the divergence function
	divergence <- function(pedigree, p0mm, p0um, p0uu, param)
	{

	  ## Initializing parameters
	  PrMM <- p0mm
	  PrUM <- p0um
	  PrUU <- p0uu
	  alpha <- param[1]
      bet <- param[2]
      weight <- param[3]
      sel    <-param[4]


	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	  element11<-(1-alpha)^2
	  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
	  element13<-(alpha^2)*sel
	  rowtotal1<-element11 + element12 + element13

	  element21<-(bet*(1-alpha))
	  element22<-((1-alpha)*(1-bet)+alpha*bet)*(1/2*(1+sel))
	  element23<-alpha*(1-bet)*sel
	  rowtotal2<-element21 + element22 + element23

	  element31<-(bet^2)
	  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
	  element33<-((1-bet))^2*sel
	  rowtotal3<-element31 + element32 + element33

	## Defining the generation (or transition) matrix
	  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
	                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
	                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)



	## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
	  Dt1t2<-NULL

		  for (p in seq_len(NROW(pedigree)))
		  {

			## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
			svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
			svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
			svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
			svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

			## Conditional divergences
			dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
								svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

			dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
								svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

			dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
								svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

			## Total (weighted) divergence
			Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


		  }

	  # Pr(UU) at equilibrium given alpha and beta
	  puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
	  puuinf.est<- puuinf.est[1,1]
  	  divout<-list(puuinf.est, Dt1t2)

	  return(divout)

	}


###### Defining the Least Square function to be minimized
###### Note the equilibrium constraint, which can be made as small as desired.

		LSE_intercept<-function(param_int)
		{
			sum((pedigree[,4] - param_int[5] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[2]])^2) +
			eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[1]]-eqp)^2)
		}



###### Calculating the initial proportions
###### We always assume that:
		# 1. p0mm is larger than actually observed. This means if p0um is available from measurements,
		#    we will just add it to p0mm.
		# 2. As a consequence of (1.) we also assume that p0um = 0.

		p0uu<-p0uu
		p0mm<-1-p0uu
		p0um<-0


   if(is.null(p0mm ==TRUE | is.null(eqp)==TRUE))
   {stop("Both eqp value AND p0mm have to be supplied")}

   if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1)
  {stop("The initial state probabilities don't sum to 1")}




##### Initializing
	optim.method<-"Nelder-Mead"
	final<-NULL
	counter<-0
	opt.out<-NULL
	pedigree<-pedigree.data


		for (s in seq_len(Nstarts) )
		{

			## Draw random starting values
			alpha.start  <-10^(runif(1, log10(10^-9), log10(10^-2)))
			beta.start   <-10^(runif(1, log10(10^-9), log10(10^-2)))
	   		weight.start <-runif(1,0,0.5)
	     	sel.start <-runif(1,0.1,1)
	    	intercept.start <-runif(1,0,max(pedigree[,4]))
	    	param_int0 = c(alpha.start, beta.start, weight.start, sel.start, intercept.start)

			## Initializing
			counter<-counter+1

			message("Progress: ", counter/Nstarts, "\n")


						opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
						alphafinal<-opt.out[1]
						betfinal<-opt.out[2]
						alphafinal<-as.numeric(opt.out[1])
    					betfinal<-as.numeric(opt.out[2])
    					weightfinal<-as.numeric(opt.out[3])
    					selfinal<-as.numeric(opt.out[4])

    					## Calculating equilibrium frequencies based on the model estimates
    					svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)


							  element11<-(1-alphafinal)^2
							  element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
							  element13<-(alphafinal^2)*selfinal
							  rowtotal1<-element11 + element12 + element13

							  element21<-(betfinal*(1-alphafinal))
							  element22<-((1-alphafinal)*(1-betfinal)+alphafinal*betfinal)*(1/2*(1+selfinal))
							  element23<-alphafinal*(1-betfinal)*selfinal
							  rowtotal2<-element21 + element22 + element23

							  element31<-(betfinal^2)
							  element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
							  element33<-((1-betfinal))^2*selfinal
							  rowtotal3<-element31 + element32 + element33

							## Defining the generation (or transition) matrix
							  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
							                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
							                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)





						## Note: This is an approximation to the equilibrium values
    						pinf.vec <- t(svGzero)  %*% ((Genmatrix)%^% 10000)
    						PrMMinf <- pinf.vec[1,3]
    						PrUMinf <- pinf.vec[1,2]
    						PrUUinf <- pinf.vec[1,1]
    						opt.out <- cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start,
    						                sel.start, intercept.start)
    						final[[s]] <- opt.out

		} # End of Nstarts loop
	  final <- do.call("rbind", final)
    colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
    colnames(final)[14:16]<-c("PrMMinf", "PrUMinf", "PrUUinf")




##### Calculating the least square of the first part of the minimized function
	 lsqpart<-NULL

	 for (l in seq_len(NROW(final)))
	 {
			  PrMM <- p0mm
			  PrUM <- p0um
	          PrUU <- p0uu
			  alpha  <- final[l, "alpha"]
			  bet    <- final[l, "beta"]
			  weight <- final[l, "weight"]
			  sel <-final[l, "sel.coef"]
			  intercept<-final[l,"intercept"]


			## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

			  element11<-(1-alpha)^2
			  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
			  element13<-(alpha^2)*sel
			  rowtotal1<-element11 + element12 + element13

			  element21<-(bet*(1-alpha))
			  element22<-((1-alpha)*(1-bet)+alpha*bet)*(1/2*(1+sel))
			  element23<-alpha*(1-bet)*sel
			  rowtotal2<-element21 + element22 + element23

			  element31<-(bet^2)
			  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
			  element33<-((1-bet))^2*sel
			  rowtotal3<-element31 + element32 + element33

			## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
			                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
			                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

					## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

					## Conditional divergences
					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

					## Total (weighted) divergence
					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


				  }


			 ## Calculating the least square part
			 lsqpart[l]<-sum((pedigree[,4] - intercept - Dt1t2)^2)
		}

	 ## Collecting results and filtering them
	 final<-cbind(final, lsqpart)
	 colnames(final)[ncol(final)]<-c("value.part")
	 final<-final[order(final[,"value"]),]
	 #index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["convcode"] == 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1)
	 index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1 & final["intercept"] > 0)
	 index.2<-setdiff(seq_len(NROW(final)) , index.1)
	 final.1<-final[index.1,]
	 final.2<-final[index.2,]


  ## Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
	 PrMM <- p0mm
	 PrUM <- p0um
	 PrUU <- p0uu
	 alpha  <- final.1[1, "alpha"]
	 bet    <- final.1[1, "beta"]
	 weight <- final.1[1, "weight"]
	 intercept<-final.1[1,"intercept"]
	 sel<-final.1[1,"sel.coef"]

	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	 svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	  		  element11<-(1-alpha)^2
			  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
			  element13<-(alpha^2)*sel
			  rowtotal1<-element11 + element12 + element13

			  element21<-(bet*(1-alpha))
			  element22<-((1-alpha)*(1-bet)+alpha*bet)*(1/2*(1+sel))
			  element23<-alpha*(1-bet)*sel
			  rowtotal2<-element21 + element22 + element23

			  element31<-(bet^2)
			  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
			  element33<-((1-bet))^2*sel
			  rowtotal3<-element31 + element32 + element33

			## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
			                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
			                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL
			  Residual<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

					## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

					## Conditional divergences
					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

					## Total (weighted) divergence
					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

				  }

			 ## Calculating the least square part
			 Residual<-(pedigree[,4] - intercept - Dt1t2)



## Augmenting pedigree
	delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
	pedigree<-cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
	colnames(pedigree)[c(4,5,6,7)]<-c("div.obs", "delta.t","div.pred", "residual")

## Making info about settings
	info<-c("p0mm", "p0um", "p0uu", "eqp", "eqp.weight", "Nstarts", "optim.method")
	info2<-c(p0mm, p0um, p0uu, eqp, eqp.weight, Nstarts, optim.method)
	info.out<-data.frame(info, info2)
	colnames(info.out)<-c("Para", "Setting")


## Generating theoretical fit

		## Reading in pedigree
			obs<-pedigree[,"div.obs"]
			dtime<-pedigree[,"delta.t"]

		## Reading in parameter estimates
			est <-final.1
			alpha <-as.numeric(est[1,1])
	    beta<-as.numeric(est[1,2])
		  weight<-as.numeric(est[1,3])
		  sel<-as.numeric(est[1,4])
	    intercept<-as.numeric(est[1,5])

		## Reading initial state vector
			settings<-info.out
			PrMM<-p0mm<-as.numeric(as.character(settings[1,2]))
			PrUM<-p0um<-as.numeric(as.character(settings[2,2]))
			PrUU<-p0uu<-as.numeric(as.character(settings[3,2]))
			time1<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time2<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time.out<-expand.grid(time1,time2)
			#time0<- rep(min(pedigree[,1]), nrow(time.out))
			time0<- rep(0, nrow(time.out))
			pedigree.new<-as.matrix(cbind(time0,time.out))
			pedigree.new<-cbind(pedigree.new, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			pedigree.new<-pedigree.new[!duplicated(pedigree.new[,4]), ]
			pedigree.new<-pedigree.new[,1:3]

		## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			svGzero   <- c(PrUU, PrMM*weight, (1-weight)*PrMM)

							alphafinal<-alpha
							betfinal<-beta
							selfinal<-sel
							interceptfinal<-intercept

							element11<-(1-alphafinal)^2
							element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
							element13<-(alphafinal^2)*selfinal
							rowtotal1<-element11 + element12 + element13

							element21<-(betfinal*(1-alphafinal))
							element22<-((1-alphafinal)*(1-betfinal)+alphafinal*betfinal)*(1/2*(1+selfinal))
							element23<-alphafinal*(1-betfinal)*selfinal
							rowtotal2<-element21 + element22 + element23

							element31<-(betfinal^2)
							element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
							element33<-((1-betfinal))^2*selfinal
							rowtotal3<-element31 + element32 + element33

							## Defining the generation (or transition) matrix
							  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
							                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
							                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)


							## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
							Dt1t2<-NULL

								for (p in seq_len(NROW(pedigree.new)))
								{

									## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
									svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,1]))
									svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))

									## Conditional divergences
									dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
												 svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

									dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
									             svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

									dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
									             svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

									## Total (weighted) divergence
									Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

								}

			pedigree.new<-cbind(pedigree.new, Dt1t2+interceptfinal, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
			pedigree.new<-pedigree.new[order(pedigree.new[,5]),]

			model<-"ABselectUUSOMA.R"

    	abfreeS.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
    	names(abfreeS.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")

## Ouputting result datasets
	dput(abfreeS.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))
	return(abfreeS.out)



} #End of function



ABselectUU<-function(pedigree.data, p0uu, eqp, eqp.weight, Nstarts, out.dir, out.name)
{


## Defining the divergence function
	divergence <- function(pedigree, p0mm, p0um, p0uu, param)
	{

	 ## Initializing parameters
	  PrMM <- p0mm
	  PrUM <- p0um
	  PrUU <- p0uu
	  alpha <- param[1]
    bet   <- param[2]
    weight <-param[3]
    sel    <-param[4]

	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	  element11<-((1-alpha)^2)
	  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
	  element13<-(alpha^2)*sel
	  rowtotal1<-element11 + element12 + element13

	  element21<-(1/4*(bet + 1 - alpha)^2)
	  element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
	  element23<-(1/4*(alpha + (1 - bet))^2)*sel
	  rowtotal2<-element21 + element22 + element23

	  element31<-(bet^2)
	  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
	  element33<-(((1-bet))^2)*sel
	  rowtotal3<-element31 + element32 + element33

	## Defining the generation (or transition) matrix
	  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
	                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
	                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

	## Calculating the expected divergence for every observed pair in 'pedigree.txt'
	  Dt1t2<-NULL

		  for (p in seq_len(NROW(pedigree)))
		  {

      		## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
      			svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
      			svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      			svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      			svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      			svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      			svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      			svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

      		## Conditional divergences
      			dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
      								svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

      			dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
      								svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

      			dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
      								svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

      		## Total (weighted) divergence
      			Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


		  }

	  ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
  	  puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
  	  puuinf.est<- puuinf.est[1,1]
  	  divout<-list(puuinf.est, Dt1t2)

	  return(divout)

	}

	## Defining the Least Square function to be minimized
		LSE_intercept<-function(param_int)
		{
		  sum((pedigree[,4] - param_int[5] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[2]])^2) +
		    eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[1]]- eqp)^2)
		}



	## Calculating the initial proportions
		p0uu<-p0uu
		p0mm<-1-p0uu
		p0um<-0

		if(is.null(p0uu ==TRUE | is.null(eqp)==TRUE))
		{stop("Both eqp value AND p0uu have to be supplied")}

		if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1)
		{stop("The initial state probabilities don't sum to 1")}

  ## Initializing
  	optim.method<-"Nelder-Mead"
  	final<-NULL
  	counter<-0
  	opt.out<-NULL
  	pedigree<-pedigree.data

		for (s in seq_len(Nstarts))
		{

    		## Draw random starting values
    			alpha.start  <-10^(runif(1, log10(10^-9), log10(10^-2)))
    			beta.start   <-10^(runif(1, log10(10^-9), log10(10^-2)))
    	    weight.start <-runif(1,0,0.1)
    	    sel.start <-runif(1,0.1,1)
    	    intercept.start <-runif(1,0,max(pedigree[,4]))
    			param_int0 = c(alpha.start, beta.start, weight.start, sel.start, intercept.start)

    		## Initializing
    			counter<-counter+1

    			message("Progress: ", counter/Nstarts, "\n")


    						opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
    						alphafinal<-as.numeric(opt.out[1])
    						betfinal<-as.numeric(opt.out[2])
    						weightfinal<-as.numeric(opt.out[3])
    						selfinal<-as.numeric(opt.out[4])

    					## Calculating equilibrium frequencies based on the model estimates
    						svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

    						element11<-((1-alphafinal)^2)
    						element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
    						element13<-(alphafinal^2)*selfinal
    						rowtotal1<-element11 + element12 + element13

    						element21<-(1/4*(betfinal + 1 - alphafinal)^2)
    						element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
    						element23<-(1/4*(alphafinal + (1 - betfinal))^2)*selfinal
    						rowtotal2<-element21 + element22 + element23

    						element31<-(betfinal^2)
    						element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
    						element33<-(((1-betfinal))^2)*selfinal
    						rowtotal3<-element31 + element32 + element33

    					## Defining the generation (or transition) matrix
    						Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
    						                      element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
    						                      element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

    					## Note: This is an approximation to the equilibrium values
    						pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
    						PrMMinf<- pinf.vec[1,3]
    						PrUMinf<- pinf.vec[1,2]
    						PrUUinf<- pinf.vec[1,1]
    						opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start,
    						                sel.start, intercept.start)
    						final[[s]] <- opt.out

		  } # End of Nstarts loop

  	final<-do.call("rbind", final)
    colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
    colnames(final)[14:16]<-c("PrMMinf", "PrUMinf", "PrUUinf")



## Calculating the least square of the first part of the minimized function
	 lsqpart<-NULL

	 for (l in seq_len(NROW(final)))
	 {
			  PrMM <- p0mm
			  PrUM <- p0um
	      PrUU <- p0uu
			  alpha  <- final[l, "alpha"]
			  bet    <- final[l, "beta"]
			  weight <- final[l, "weight"]
			  sel <-final[l, "sel.coef"]
			  intercept<-final[l,"intercept"]

			 ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

			  element11<-((1-alpha)^2)
			  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
			  element13<-(alpha^2)*sel
			  rowtotal1<-element11 + element12 + element13

			  element21<-(1/4*(bet + 1 - alpha)^2)
			  element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
			  element23<-(1/4*(alpha + (1 - bet))^2)*sel
			  rowtotal2<-element21 + element22 + element23

			  element31<-(bet^2)
			  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
			  element33<-(((1-bet))^2)*sel
			  rowtotal3<-element31 + element32 + element33

			 ## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
			                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
			                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

			 ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL

			  for (p in seq_len(NROW(pedigree)))
			  {

    			   ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
    			    svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
    			    svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    			    svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    			    svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    			    svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    			    svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    			    svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))


    			   ## Conditional divergences
    			    dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
    			                        svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

    			    dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
    			                        svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

    			    dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
    			                        svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

    			   ## Total (weighted) divergence
    			    Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

			  } # End of pedigree loop


			## Calculating the least square part
			 lsqpart[l]<-sum((pedigree[,4] - intercept - Dt1t2)^2)

		} # End of final loop

	## Collecting results and filtering them
	 final<-cbind(final, lsqpart)
	 colnames(final)[ncol(final)]<-c("value.part")
	 final<-final[order(final[,"value"]),]
	 #index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["convcode"] == 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1)
	 index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1 & final["intercept"] > 0)
	 index.2<-setdiff(seq_len(NROW(final)), index.1)
	 final.1<-final[index.1,]
	 final.2<-final[index.2,]


  ## Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
	 PrMM <- p0mm
	 PrUM <- p0um
	 PrUU <- p0uu
	 alpha  <- final.1[1, "alpha"]
	 bet    <- final.1[1, "beta"]
	 weight <- final.1[1, "weight"]
	 intercept<-final.1[1,"intercept"]
	 sel<-final.1[1,"sel.coef"]

	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	 svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	 element11<-((1-alpha)^2)
	 element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
	 element13<-(alpha^2)*sel
	 rowtotal1<-element11 + element12 + element13

	 element21<-(1/4*(bet + 1 - alpha)^2)
	 element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
	 element23<-(1/4*(alpha + (1 - bet))^2)*sel
	 rowtotal2<-element21 + element22 + element23

	 element31<-(bet^2)
	 element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
	 element33<-(((1-bet))^2)*sel
	 rowtotal3<-element31 + element32 + element33

	## Defining the generation (or transition) matrix
	 Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
	                       element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
	                       element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

	## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
		Dt1t2<-NULL
		Residual<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

    				## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
    					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
    					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

    				## Conditional divergences
    					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
    										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

    					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
    										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

    					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
    										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

    				## Total (weighted) divergence
    					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

				  }

			## Calculating the least square part
			 Residual<-(pedigree[,4] - intercept - Dt1t2)


## Augmenting pedigree
	delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
	pedigree<-cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
	colnames(pedigree)[c(4,5,6,7)]<-c("div.obs", "delta.t","div.pred", "residual")

## Making info about settings
	info<-c("p0mm", "p0um", "p0uu", "eqp", "eqp.weight", "Nstarts", "optim.method")
	info2<-c(p0mm, p0um, p0uu, eqp, eqp.weight, Nstarts, optim.method)
	info.out<-data.frame(info, info2)
	colnames(info.out)<-c("Para", "Setting")


## Generating theoretical fit

		## Reading in pedigree
			obs<-pedigree[,"div.obs"]
			dtime<-pedigree[,"delta.t"]

		## Reading in parameter estimates
			est <-final.1
			alpha <-as.numeric(est[1,1])
	    beta<-as.numeric(est[1,2])
		  weight<-as.numeric(est[1,3])
		  sel<-as.numeric(est[1,4])
	    intercept<-as.numeric(est[1,5])

		## Reading initial state vector
			settings<-info.out
			PrMM<-p0mm<-as.numeric(as.character(settings[1,2]))
			PrUM<-p0um<-as.numeric(as.character(settings[2,2]))
			PrUU<-p0uu<-as.numeric(as.character(settings[3,2]))
			time1<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time2<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time.out<-expand.grid(time1,time2)
			#time0<- rep(min(pedigree[,1]), nrow(time.out))
			time0<- rep(0, nrow(time.out))
			pedigree.new<-as.matrix(cbind(time0,time.out))
			pedigree.new<-cbind(pedigree.new, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			pedigree.new<-pedigree.new[!duplicated(pedigree.new[,4]), ]
			pedigree.new<-pedigree.new[,1:3]

		## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			svGzero   <- c(PrUU, PrMM*weight, (1-weight)*PrMM)

							alphafinal<-alpha
							betfinal<-beta
							selfinal<-sel
							interceptfinal<-intercept

							element11<-((1-alphafinal)^2)
							element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
							element13<-(alphafinal^2)*selfinal
							rowtotal1<-element11 + element12 + element13

							element21<-(1/4*(betfinal + 1 - alphafinal)^2)
							element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
							element23<-(1/4*(alphafinal + (1 - betfinal))^2)*selfinal
							rowtotal2<-element21 + element22 + element23

							element31<-(betfinal^2)
							element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
							element33<-(((1-betfinal))^2)*selfinal
							rowtotal3<-element31 + element32 + element33

						## Defining the generation (or transition) matrix
							Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
							                      element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
							                      element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

						## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
							Dt1t2<-NULL

								for (p in seq_len(NROW(pedigree.new)))
								{
  								## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
  									svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,1]))
  									svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
  									svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
  									svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
  									svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
  									svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
  									svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))

  								## Conditional divergences
  									dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
  												 svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

  									dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
  									             svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

  									dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
  									             svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

  								## Total (weighted) divergence
  									Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

								}

			pedigree.new<-cbind(pedigree.new, Dt1t2+interceptfinal, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
			pedigree.new<-pedigree.new[order(pedigree.new[,5]),]

			model<-"ABselectUU.R"

    	abfreeS.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
    	names(abfreeS.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")

## Ouputting result datasets
	dput(abfreeS.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))
	return(abfreeS.out)


} #End of function






ABselectMMSOMA<-function(pedigree.data, p0uu, eqp, eqp.weight, Nstarts, out.dir, out.name)
{
 allow.neg.intercept="no"

##### Defining the divergence function
	divergence <- function(pedigree, p0mm, p0um, p0uu, param)
	{

	  ## Initializing parameters
	  PrMM <- p0mm
	  PrUM <- p0um
	  PrUU <- p0uu
	  alpha <- param[1]
      bet <- param[2]
      weight <- param[3]
      sel    <-param[4]


	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	  element11<-(1-alpha)^2*sel
	  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
	  element13<-(alpha^2)
	  rowtotal1<-element11 + element12 + element13

	  element21<-(bet*(1-alpha))*sel
	  element22<-((1-alpha)*(1-bet)+alpha*bet)*(1/2*(1+sel))
	  element23<-alpha*(1-bet)
	  rowtotal2<-element21 + element22 + element23

	  element31<-(bet^2)*sel
	  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
	  element33<-((1-bet))^2
	  rowtotal3<-element31 + element32 + element33

	## Defining the generation (or transition) matrix
	  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
	                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
	                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)



	## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
	  Dt1t2<-NULL

		  for (p in seq_len(NROW(pedigree)))
		  {

			## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
			svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
			svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
			svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
			svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

			## Conditional divergences
			dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
								svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

			dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
								svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

			dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
								svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

			## Total (weighted) divergence
			Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


		  }

	  # Pr(UU) at equilibrium given alpha and beta
	  puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
	  puuinf.est<- puuinf.est[1,1]
  	  divout<-list(puuinf.est, Dt1t2)

	  return(divout)

	}


###### Defining the Least Square function to be minimized
###### Note the equilibrium constraint, which can be made as small as desired.

		LSE_intercept<-function(param_int)
		{
			sum((pedigree[,4] - param_int[5] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[2]])^2) +
			eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[1]]-eqp)^2)
		}



###### Calculating the initial proportions
###### We always assume that:
		# 1. p0mm is larger than actually observed. This means if p0um is available from measurements,
		#    we will just add it to p0mm.
		# 2. As a consequence of (1.) we also assume that p0um = 0.

		p0uu<-p0uu
		p0mm<-1-p0uu
		p0um<-0


   if(is.null(p0mm ==TRUE | is.null(eqp)==TRUE))
   {stop("Both eqp value AND p0mm have to be supplied")}

   if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1)
  {stop("The initial state probabilities don't sum to 1")}




##### Initializing
	optim.method<-"Nelder-Mead"
	final<-NULL
	counter<-0
	opt.out<-NULL
	pedigree<-pedigree.data

		for (s in seq_len(Nstarts) )
		{

			## Draw random starting values
			alpha.start  <-10^(runif(1, log10(10^-9), log10(10^-2)))
			beta.start   <-10^(runif(1, log10(10^-9), log10(10^-2)))
      weight.start <-runif(1,0,0.5)
      sel.start <-runif(1,0.1,1)
      intercept.start <-runif(1,0,max(pedigree[,4]))
      param_int0 = c(alpha.start, beta.start, weight.start, sel.start, intercept.start)
			## Initializing
			counter<-counter+1
			message("Progress: ", counter/Nstarts, "\n")
						opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
						alphafinal<-opt.out[1]
						betfinal<-opt.out[2]
						alphafinal<-as.numeric(opt.out[1])
						betfinal<-as.numeric(opt.out[2])
            weightfinal<-as.numeric(opt.out[3])
            selfinal<-as.numeric(opt.out[4])

    					## Calculating equilibrium frequencies based on the model estimates
    					svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)


							  element11<-(1-alphafinal)^2*selfinal
							  element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
							  element13<-(alphafinal^2)
							  rowtotal1<-element11 + element12 + element13

							  element21<-(betfinal*(1-alphafinal))*selfinal
							  element22<-((1-alphafinal)*(1-betfinal)+alphafinal*betfinal)*(1/2*(1+selfinal))
							  element23<-alphafinal*(1-betfinal)
							  rowtotal2<-element21 + element22 + element23

							  element31<-(betfinal^2)*selfinal
							  element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
							  element33<-((1-betfinal))^2
							  rowtotal3<-element31 + element32 + element33

							## Defining the generation (or transition) matrix
							  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
							                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
							                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)





						## Note: This is an approximation to the equilibrium values
    						pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
    						PrMMinf<- pinf.vec[1,3]
    						PrUMinf<- pinf.vec[1,2]
    						PrUUinf<- pinf.vec[1,1]
    						opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start,
    						                sel.start, intercept.start)
    						final[[s]] <- opt.out

		} # End of Nstarts loop

	  final <- do.call("rbind", final)
    colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
    colnames(final)[14:16]<-c("PrMMinf", "PrUMinf", "PrUUinf")




##### Calculating the least square of the first part of the minimized function
	 lsqpart<-NULL

	 for (l in seq_len(NROW(final)))
	 {
			  PrMM <- p0mm
			  PrUM <- p0um
	          PrUU <- p0uu
			  alpha  <- final[l, "alpha"]
			  bet    <- final[l, "beta"]
			  weight <- final[l, "weight"]
			  sel <-final[l, "sel.coef"]
			  intercept<-final[l,"intercept"]


			## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

			  element11<-(1-alpha)^2*sel
			  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
			  element13<-(alpha^2)
			  rowtotal1<-element11 + element12 + element13

			  element21<-(bet*(1-alpha))*sel
			  element22<-((1-alpha)*(1-bet)+alpha*bet)*(1/2*(1+sel))
			  element23<-alpha*(1-bet)
			  rowtotal2<-element21 + element22 + element23

			  element31<-(bet^2)*sel
			  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
			  element33<-((1-bet))^2
			  rowtotal3<-element31 + element32 + element33

			## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
			                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
			                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

					## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

					## Conditional divergences
					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

					## Total (weighted) divergence
					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


				  }


			 ## Calculating the least square part
			 lsqpart[l]<-sum((pedigree[,4] - intercept - Dt1t2)^2)
		}

	 ## Collecting results and filtering them
	 final<-cbind(final, lsqpart)
	 colnames(final)[ncol(final)]<-c("value.part")
	 final<-final[order(final[,"value"]),]
	 #index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["convcode"] == 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1)
	 index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1 & final["intercept"] > 0)
	 index.2<-setdiff(seq_len(NROW(final)), index.1)
	 final.1<-final[index.1,]
	 final.2<-final[index.2,]


  ## Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
	 PrMM <- p0mm
	 PrUM <- p0um
	 PrUU <- p0uu
	 alpha  <- final.1[1, "alpha"]
	 bet    <- final.1[1, "beta"]
	 weight <- final.1[1, "weight"]
	 intercept<-final.1[1,"intercept"]
	 sel<-final.1[1,"sel.coef"]

	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	 svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	  		element11<-(1-alpha)^2*sel
			  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
			  element13<-(alpha^2)
			  rowtotal1<-element11 + element12 + element13

			  element21<-(bet*(1-alpha))*sel
			  element22<-((1-alpha)*(1-bet)+alpha*bet)*(1/2*(1+sel))
			  element23<-alpha*(1-bet)
			  rowtotal2<-element21 + element22 + element23

			  element31<-(bet^2)*sel
			  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
			  element33<-((1-bet))^2
			  rowtotal3<-element31 + element32 + element33

			## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
			                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
			                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL
			  Residual<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

					## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

					## Conditional divergences
					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

					## Total (weighted) divergence
					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

				  }

			 ## Calculating the least square part
			 Residual<-(pedigree[,4] - intercept - Dt1t2)



## Augmenting pedigree
	delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
	pedigree<-cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
	colnames(pedigree)[c(4,5,6,7)]<-c("div.obs", "delta.t","div.pred", "residual")

## Making info about settings
	info<-c("p0mm", "p0um", "p0uu", "eqp", "eqp.weight", "Nstarts", "optim.method")
	info2<-c(p0mm, p0um, p0uu, eqp, eqp.weight, Nstarts, optim.method)
	info.out<-data.frame(info, info2)
	colnames(info.out)<-c("Para", "Setting")


## Generating theoretical fit

		## Reading in pedigree
			obs<-pedigree[,"div.obs"]
			dtime<-pedigree[,"delta.t"]

		## Reading in parameter estimates
			est <-final.1
			alpha <-as.numeric(est[1,1])
	    beta<-as.numeric(est[1,2])
		  weight<-as.numeric(est[1,3])
		  sel<-as.numeric(est[1,4])
	    intercept<-as.numeric(est[1,5])

		## Reading initial state vector
			settings<-info.out
			PrMM<-p0mm<-as.numeric(as.character(settings[1,2]))
			PrUM<-p0um<-as.numeric(as.character(settings[2,2]))
			PrUU<-p0uu<-as.numeric(as.character(settings[3,2]))
			time1<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time2<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time.out<-expand.grid(time1,time2)
			#time0<- rep(min(pedigree[,1]), nrow(time.out))
			time0<- rep(0, nrow(time.out))
			pedigree.new<-as.matrix(cbind(time0,time.out))
			pedigree.new<-cbind(pedigree.new, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			pedigree.new<-pedigree.new[!duplicated(pedigree.new[,4]), ]
			pedigree.new<-pedigree.new[,1:3]

		## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			svGzero   <- c(PrUU, PrMM*weight, (1-weight)*PrMM)

							alphafinal<-alpha
							betfinal<-beta
							selfinal<-sel
							interceptfinal<-intercept

							element11<-(1-alphafinal)^2*selfinal
							element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
							element13<-(alphafinal^2)
							rowtotal1<-element11 + element12 + element13

							element21<-(betfinal*(1-alphafinal))*selfinal
							element22<-((1-alphafinal)*(1-betfinal)+alphafinal*betfinal)*(1/2*(1+selfinal))
							element23<-alphafinal*(1-betfinal)
							rowtotal2<-element21 + element22 + element23

							element31<-(betfinal^2)*selfinal
							element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
							element33<-((1-betfinal))^2
							rowtotal3<-element31 + element32 + element33

							## Defining the generation (or transition) matrix
							  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
							                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
							                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)


							## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
							Dt1t2<-NULL

								for (p in seq_len(NROW(pedigree.new)))
								{

									## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
									svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,1]))
									svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))

									## Conditional divergences
									dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
												 svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

									dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
									             svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

									dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
									             svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

									## Total (weighted) divergence
									Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

								}

			pedigree.new<-cbind(pedigree.new, Dt1t2+interceptfinal, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
			pedigree.new<-pedigree.new[order(pedigree.new[,5]),]

			model<-"ABselectMMSOMA.R"

    	abfreeS.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
    	names(abfreeS.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")

## Ouputting result datasets
	dput(abfreeS.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))
	return(abfreeS.out)



} #End of function



ABselectMM<-function(pedigree.data, p0uu, eqp, eqp.weight, Nstarts, out.dir, out.name)
{


## Defining the divergence function
	divergence <- function(pedigree, p0mm, p0um, p0uu, param)
	{

	 ## Initializing parameters
	  PrMM <- p0mm
	  PrUM <- p0um
	  PrUU <- p0uu
	  alpha <- param[1]
    bet   <- param[2]
    weight <-param[3]
    sel    <-param[4]

	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	  element11<-((1-alpha)^2)*sel
	  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
	  element13<-(alpha^2)
	  rowtotal1<-element11 + element12 + element13

	  element21<-(1/4*(bet + 1 - alpha)^2)*sel
	  element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
	  element23<-(1/4*(alpha + (1 - bet))^2)
	  rowtotal2<-element21 + element22 + element23

	  element31<-(bet^2)*sel
	  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
	  element33<-((1-bet))^2
	  rowtotal3<-element31 + element32 + element33

	## Defining the generation (or transition) matrix
	  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
	                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
	                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

	## Calculating the expected divergence for every observed pair in 'pedigree.txt'
	  Dt1t2<-NULL

		  for (p in seq_len(NROW(pedigree)))
		  {

      		## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
      			svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
      			svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      			svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      			svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      			svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      			svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      			svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

      		## Conditional divergences
      			dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
      								svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

      			dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
      								svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

      			dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
      								svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

      		## Total (weighted) divergence
      			Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


		  }

	  ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
  	  puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
  	  puuinf.est<- puuinf.est[1,1]
  	  divout<-list(puuinf.est, Dt1t2)

	  return(divout)

	}

	## Defining the Least Square function to be minimized
		LSE_intercept<-function(param_int)
		{
		  sum((pedigree[,4] - param_int[5] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[2]])^2) +
		    eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[1]]- eqp)^2)
		}



	## Calculating the initial proportions
		p0uu<-p0uu
		p0mm<-1-p0uu
		p0um<-0

		if(is.null(p0uu ==TRUE | is.null(eqp)==TRUE))
		{stop("Both eqp value AND p0uu have to be supplied")}

		if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1)
		{stop("The initial state probabilities don't sum to 1")}

  ## Initializing
  	optim.method<-"Nelder-Mead"
  	final<-NULL
  	counter<-0
  	opt.out<-NULL
  	pedigree<-pedigree.data

		for (s in seq_len(Nstarts))
		{

    		## Draw random starting values
    			alpha.start  <-10^(runif(1, log10(10^-9), log10(10^-2)))
    			beta.start   <-10^(runif(1, log10(10^-9), log10(10^-2)))
    	    weight.start <-runif(1,0,0.1)
    	    sel.start <-runif(1,0.1,1)
    	    intercept.start <-runif(1,0,max(pedigree[,4]))
    			param_int0 = c(alpha.start, beta.start, weight.start, sel.start, intercept.start)

    		## Initializing
    			counter<-counter+1

    			message("Progress: ", counter/Nstarts, "\n")


    						opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
    						alphafinal<-as.numeric(opt.out[1])
    						betfinal<-as.numeric(opt.out[2])
    						weightfinal<-as.numeric(opt.out[3])
    						selfinal<-as.numeric(opt.out[4])

    					## Calculating equilibrium frequencies based on the model estimates
    						svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

    						element11<-((1-alphafinal)^2)*selfinal
    						element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
    						element13<-(alphafinal^2)
    						rowtotal1<-element11 + element12 + element13

    						element21<-(1/4*(betfinal + 1 - alphafinal)^2)*selfinal
    						element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
    						element23<-(1/4*(alphafinal + (1 - betfinal))^2)
    						rowtotal2<-element21 + element22 + element23

    						element31<-(betfinal^2)*selfinal
    						element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
    						element33<-((1-betfinal))^2
    						rowtotal3<-element31 + element32 + element33

    					## Defining the generation (or transition) matrix
    						Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
    						                      element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
    						                      element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

    					## Note: This is an approximation to the equilibrium values
    						pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
    						PrMMinf<- pinf.vec[1,3]
    						PrUMinf<- pinf.vec[1,2]
    						PrUUinf<- pinf.vec[1,1]
    						opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start,
    						                sel.start, intercept.start)
    						final[[s]]<-opt.out

		  } # End of Nstarts loop

    final<-do.call("rbind", final)
    colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
    colnames(final)[14:16]<-c("PrMMinf", "PrUMinf", "PrUUinf")



## Calculating the least square of the first part of the minimized function
	 lsqpart<-NULL

	 for (l in seq_len(NROW(final)))
	 {
			  PrMM <- p0mm
			  PrUM <- p0um
	      PrUU <- p0uu
			  alpha  <- final[l, "alpha"]
			  bet    <- final[l, "beta"]
			  weight <- final[l, "weight"]
			  sel <-final[l, "sel.coef"]
			  intercept<-final[l,"intercept"]

			 ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

			  element11<-((1-alpha)^2)*sel
			  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
			  element13<-(alpha^2)
			  rowtotal1<-element11 + element12 + element13

			  element21<-(1/4*(bet + 1 - alpha)^2)*sel
			  element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
			  element23<-(1/4*(alpha + (1 - bet))^2)
			  rowtotal2<-element21 + element22 + element23

			  element31<-(bet^2)*sel
			  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
			  element33<-((1-bet))^2
			  rowtotal3<-element31 + element32 + element33

			 ## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
			                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
			                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

			 ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL

			  for (p in seq_len(NROW(pedigree)))
			  {

    			   ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
    			    svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
    			    svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    			    svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    			    svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    			    svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    			    svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    			    svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))


    			   ## Conditional divergences
    			    dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
    			                        svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

    			    dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
    			                        svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

    			    dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
    			                        svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

    			   ## Total (weighted) divergence
    			    Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

			  } # End of pedigree loop


			## Calculating the least square part
			 lsqpart[l]<-sum((pedigree[,4] - intercept - Dt1t2)^2)

		} # End of final loop

	## Collecting results and filtering them
	 final<-cbind(final, lsqpart)
	 colnames(final)[ncol(final)]<-c("value.part")
	 final<-final[order(final[,"value"]),]
	 #index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["convcode"] == 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1)
	 index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1 & final["intercept"] > 0)
	 index.2<-setdiff(seq_len(NROW(final)), index.1)
	 final.1<-final[index.1,]
	 final.2<-final[index.2,]


  ## Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
	 PrMM <- p0mm
	 PrUM <- p0um
	 PrUU <- p0uu
	 alpha  <- final.1[1, "alpha"]
	 bet    <- final.1[1, "beta"]
	 weight <- final.1[1, "weight"]
	 intercept<-final.1[1,"intercept"]
	 sel<-final.1[1,"sel.coef"]

	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	 svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	 element11<-((1-alpha)^2)*sel
	 element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
	 element13<-(alpha^2)
	 rowtotal1<-element11 + element12 + element13

	 element21<-(1/4*(bet + 1 - alpha)^2)*sel
	 element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
	 element23<-(1/4*(alpha + (1 - bet))^2)
	 rowtotal2<-element21 + element22 + element23

	 element31<-(bet^2)*sel
	 element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
	 element33<-((1-bet))^2
	 rowtotal3<-element31 + element32 + element33

	## Defining the generation (or transition) matrix
	 Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
	                       element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
	                       element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

	## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
		Dt1t2<-NULL
		Residual<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

    				## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
    					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
    					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

    				## Conditional divergences
    					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
    										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

    					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
    										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

    					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
    										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

    				## Total (weighted) divergence
    					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

				  }

			## Calculating the least square part
			 Residual<-(pedigree[,4] - intercept - Dt1t2)


## Augmenting pedigree
	delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
	pedigree<-cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
	colnames(pedigree)[c(4,5,6,7)]<-c("div.obs", "delta.t","div.pred", "residual")

## Making info about settings
	info<-c("p0mm", "p0um", "p0uu", "eqp", "eqp.weight", "Nstarts", "optim.method")
	info2<-c(p0mm, p0um, p0uu, eqp, eqp.weight, Nstarts, optim.method)
	info.out<-data.frame(info, info2)
	colnames(info.out)<-c("Para", "Setting")


## Generating theoretical fit

		## Reading in pedigree
			obs<-pedigree[,"div.obs"]
			dtime<-pedigree[,"delta.t"]

		## Reading in parameter estimates
			est <-final.1
			alpha <-as.numeric(est[1,1])
	    beta<-as.numeric(est[1,2])
		  weight<-as.numeric(est[1,3])
		  sel<-as.numeric(est[1,4])
	    intercept<-as.numeric(est[1,5])

		## Reading initial state vector
			settings<-info.out
			PrMM<-p0mm<-as.numeric(as.character(settings[1,2]))
			PrUM<-p0um<-as.numeric(as.character(settings[2,2]))
			PrUU<-p0uu<-as.numeric(as.character(settings[3,2]))
			time1<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time2<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time.out<-expand.grid(time1,time2)
			#time0<- rep(min(pedigree[,1]), nrow(time.out))
			time0<- rep(0, nrow(time.out))
			pedigree.new<-as.matrix(cbind(time0,time.out))
			pedigree.new<-cbind(pedigree.new, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			pedigree.new<-pedigree.new[!duplicated(pedigree.new[,4]), ]
			pedigree.new<-pedigree.new[,1:3]

		## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			svGzero   <- c(PrUU, PrMM*weight, (1-weight)*PrMM)

							alphafinal<-alpha
							betfinal<-beta
							selfinal<-sel
							interceptfinal<-intercept

							element11<-((1-alphafinal)^2)*selfinal
							element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
							element13<-(alphafinal^2)
							rowtotal1<-element11 + element12 + element13

							element21<-(1/4*(betfinal + 1 - alphafinal)^2)*selfinal
							element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
							element23<-(1/4*(alphafinal + (1 - betfinal))^2)
							rowtotal2<-element21 + element22 + element23

							element31<-(betfinal^2)*selfinal
							element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
							element33<-((1-betfinal))^2
							rowtotal3<-element31 + element32 + element33

						## Defining the generation (or transition) matrix
							Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
							                      element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
							                      element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

						## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
							Dt1t2<-NULL

								for (p in seq_len(NROW(pedigree.new)))
								{
  								## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
  									svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,1]))
  									svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
  									svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
  									svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
  									svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
  									svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
  									svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))

  								## Conditional divergences
  									dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
  												 svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

  									dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
  									             svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

  									dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
  									             svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

  								## Total (weighted) divergence
  									Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

								}

			pedigree.new<-cbind(pedigree.new, Dt1t2+interceptfinal, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
			pedigree.new<-pedigree.new[order(pedigree.new[,5]),]

			model<-"ABselectMM.R"

    	abfreeS.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
    	names(abfreeS.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")

## Ouputting result datasets
	dput(abfreeS.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))
	return(abfreeS.out)


} #End of function



ABplot <- function(pedigree.names, output.dir, out.name, alpha=0.5, geom.point.size=2,
  geom.line.size=0.9, plot.height=8, plot.width=11, plot.type="both", lsq.line="theory", intract=FALSE) {

  #initialize vectors
  dfDiv <- c()
  predDfFits <- c()
  theoryDfFits <- c()
  gg <- ggplot()

  #NCOL and NROW do the same, treating a vector as 1-column matrix.

    #for(j in 1:NCOL(pedigree.names)) {

    datain <- dget(pedigree.names)
    name <- gsub("(_).*", "", basename(as.character(pedigree.names)))
    category <- gsub(pattern=paste0(name, "_|\\.Rdata$"), "", basename(as.character(pedigree.names)))

    div <- as.data.frame(datain$pedigree[,c("delta.t","div.obs")])
    name.div <- cbind(category, name, div)
    dfDiv <- rbind(dfDiv, name.div)
    intercept <- datain$estimates[1, "intercept"]

    #predictive fit
    pred.fit.data <- aggregate(datain$pedigree[,"div.pred"], by=list(datain$pedigree[,"delta.t"]), median)
    pred.fits <- c(intercept, pred.fit.data[,2])
    pred.fit.t <- c(0, pred.fit.data[,1])

    predDfFits <- rbind(predDfFits, data.frame(category, name, pred.fit.t, pred.fits))

    #theoritical fit
    theory.fit.data <- datain$for.fit.plot
    theory.fits <- c(intercept, theory.fit.data[,"div.sim"])
    theory.fit.t <- c(0, theory.fit.data[,"delta.t"])

    theoryDfFits <- rbind(theoryDfFits, data.frame(category, name, theory.fit.t, theory.fits))

  # }

  if ((plot.type == "data.only")) {
    gg <- gg + geom_point(data=dfDiv, aes(x=delta.t, y=div.obs, colour=name), alpha = alpha, size=geom.point.size)
  }

  if ((plot.type == "fit.only") && (lsq.line == "pred")) {
    gg <- gg + geom_line(data=predDfFits, aes(x=pred.fit.t, y=pred.fits, color=name), size=geom.line.size)
  }

  if ((plot.type == "fit.only") && (lsq.line == "theory")) {
    gg <- gg + geom_line(data=theoryDfFits, aes(x=theory.fit.t, y=theory.fits, color=name), size=geom.line.size)
  }

  if ((plot.type == "both") && (lsq.line == "pred")) {
    gg <- gg + geom_point(data=dfDiv, aes(x=delta.t, y=div.obs, colour=name), alpha = alpha, size=geom.point.size) +
      geom_line(data=predDfFits, aes(x=pred.fit.t, y=pred.fits, color=name), size=geom.line.size)
  }

  if ((plot.type == "both") && (lsq.line == "theory")) {
    gg <- gg + geom_point(data=dfDiv, aes(x=delta.t, y=div.obs, colour=name), alpha = alpha, size=geom.point.size) +
      geom_line(data=theoryDfFits, aes(x=theory.fit.t, y=theory.fits, color=name), size=geom.line.size)
  }

  #adjust x, y-axis
  ymax <- max(dfDiv$div.obs) + 0.005
  xmax <- max(dfDiv$delta.t) + 10


  gg <- gg + labs(x="delta t (generations)", y="methylation divergence", color="category") +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(expand = c(0, 0), limits=c(0, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0, ymax)) +
    theme(legend.position="none")
  gg <- gg  + theme_test()

  if(!intract){
  pdf(paste0(output.dir,"/", out.name, ".pdf", sep=""), colormodel = 'cmyk', width = plot.width, height = plot.height)
  print(gg)
  dev.off()
  }else{
    pl <- ggplotly(gg)
    pl
  }

}



ABnull<-function(pedigree.data, out.dir, out.name)
{

    pedigree<-pedigree.data
    delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
    div.pred<-mean(pedigree[,"D.value"])
    residual<-pedigree[,"D.value"] - div.pred
    pedigree<-cbind(pedigree, delta.t, div.pred, residual)
    colnames(pedigree)[4:7]<-c("div.obs", "delta.t", "div.pred", "residual")

    model<-"ABnull.R"
    final.1<-div.pred
    final.2<-NULL
    info.out<-NULL
    pedigree.new<-cbind(pedigree[,1:3], div.pred, delta.t)
    colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
    abfree.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
    names(abfree.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")

    ## Ouputting result datasets
    dput(abfree.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))
    return(abfree.out)

} #End of function



ABneutralSOMA<-function(pedigree.data, p0uu, eqp, eqp.weight, Nstarts, out.dir, out.name)
{

 allow.neg.intercept="no"

##### Defining the divergence function
	divergence <- function(pedigree, p0mm, p0um, p0uu, param)
	{

	  ## Initializing parameters
	  PrMM <- p0mm
	  PrUM <- p0um
	  PrUU <- p0uu
	  alpha <- param[1]
    bet <- param[2]
    weight <- param[3]


	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)



	## Defining the generation (or transition) matrix for the mitotic case
	  Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha,alpha^2,
	                        bet*(1-alpha), (1-alpha)*(1-bet)+alpha*bet, alpha*(1-bet),
							bet^2, 2*(1-bet)*bet, (1-bet)^2),nrow=3, byrow=TRUE)


	## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
	  Dt1t2<-NULL

		  for (p in seq_len(NROW(pedigree)))
		  {

			## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
			svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
			svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
			svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
			svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

			## Conditional divergences
			dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
								svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

			dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
								svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

			dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
								svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

			## Total (weighted) divergence
			Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


		  }

	  # Pr(UU) at equilibrium given alpha and beta
	  puuinf.est<-(bet^2)/((alpha+bet)^2)
	  divout<-list(puuinf.est, Dt1t2)

	  return(divout)

	}


###### Defining the Least Square function to be minimized
###### Note the equilibrium constraint, which can be made as small as desired.

		LSE_intercept<-function(param_int)
		{
			sum((pedigree[,4] - param_int[4] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[2]])^2) +
			eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[1]]-eqp)^2)
		}



###### Calculating the initial proportions
###### We always assume that:
		# 1. p0mm is larger than actually observed. This means if p0um is available from measurements,
		#    we will just add it to p0mm.
		# 2. As a consequence of (1.) we also assume that p0um = 0.

		p0uu<-p0uu
		p0mm<-1-p0uu
		p0um<-0


   if(is.null(p0mm ==TRUE | is.null(eqp)==TRUE))
   {stop("Both eqp value AND p0mm have to be supplied")}

   if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1)
  {stop("The initial state probabilities don't sum to 1")}




##### Initializing
	optim.method<-"Nelder-Mead"
	final<-NULL
	counter<-0
	opt.out<-NULL
	pedigree<-pedigree.data


		for (s in seq_len(Nstarts))
		{

			## Draw random starting values
			alpha.start  <-10^(runif(1, log10(10^-9), log10(10^-2)))
			beta.start   <-10^(runif(1, log10(10^-9), log10(10^-2)))
	    weight.start <-runif(1,0,0.5)
	    intercept.start <-runif(1,0,max(pedigree[,4]))
			param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)

			## Initializing
			counter<-counter+1

			message("Progress: ", counter/Nstarts, "\n")


						opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
						alphafinal<-opt.out[1]
						betfinal<-opt.out[2]
						PrMMinf <- (alphafinal^2)/((alphafinal+betfinal)^2)
						PrUMinf <- (2*alphafinal*betfinal)/((alphafinal+betfinal)^2)
						PrUUinf <- (betfinal^2)/((alphafinal+betfinal)^2)
						opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
						final[[s]] <- opt.out


		} # End of Nstarts loop
    final <- do.call("rbind", final)
    colnames(final)[1:4]<-c("alpha", "beta", "weight", "intercept")
    colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")




##### Calculating the least square of the first part of the minimized function
	 lsqpart<-NULL

	 for (l in seq_len(NROW(final)))
	 {
			  PrMM <- p0mm
			  PrUM <- p0um
	      PrUU <- p0uu
			  alpha  <- final[l, "alpha"]
			  bet    <- final[l, "beta"]
			  weight <- final[l, "weight"]
			  intercept<-final[l,"intercept"]


			## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)


			  ## Defining the generation (or transition) matrix for the mitotic case
				Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha,alpha^2,
	                        bet*(1-alpha), (1-alpha)*(1-bet)+alpha*bet, alpha*(1-bet),
							bet^2, 2*(1-bet)*bet, (1-bet)^2),nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

					## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

					## Conditional divergences
					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

					## Total (weighted) divergence
					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


				  }


			 ## Calculating the least square part
			 lsqpart[l]<-sum((pedigree[,4] - intercept - Dt1t2)^2)
		}

	 final<-cbind(final, lsqpart)
	 colnames(final)[ncol(final)]<-c("value.part")
	 final<-final[order(final[,"value"]),]

	  if (allow.neg.intercept == "yes")
	  { index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["convcode"] == 0)}

	  if (allow.neg.intercept == "no")
	  {index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["intercept"] > 0 & final["convcode"] == 0)}


	 #index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["intercept"] > 0)
	 index.2<-setdiff(seq_len(NROW(final)), index.1)
	 final.1<-final[index.1,]
	 final.2<-final[index.2,]



##### Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
	 PrMM <- p0mm
	 PrUM <- p0um
	 PrUU <- p0uu
	 alpha  <- final.1[1, "alpha"]
	 bet    <- final.1[1, "beta"]
	 weight <- final.1[1, "weight"]
	 intercept<-final.1[1,"intercept"]


			 ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)



			  ## Defining the generation (or transition) matrix for the mitotic case
				Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha,alpha^2,
	                        bet*(1-alpha), (1-alpha)*(1-bet)+alpha*bet, alpha*(1-bet),
							bet^2, 2*(1-bet)*bet, (1-bet)^2),nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL
			  Residual<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

					## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

					## Conditional divergences
					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

					## Total (weighted) divergence
					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

				  }

			 ## Calculating the least square part
			 Residual<-(pedigree[,4] - intercept - Dt1t2)



##### Augmenting pedigree
	delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
	#pedigree<-cbind(pedigree,delta.t)
	pedigree<-cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
	colnames(pedigree)[c(4,5,6,7)]<-c("div.obs", "delta.t","div.pred", "residual")


##### Making info about settings
		info<-c("p0mm", "p0um", "p0uu", "eqp", "eqp.weight", "Nstarts", "optim.method")
		info2<-c(p0mm, p0um, p0uu, eqp, eqp.weight, Nstarts, optim.method)
		info.out<-data.frame(info, info2)
		colnames(info.out)<-c("Para", "Setting")






###### Generating theoretical fit

			## Reading in pedigree
			obs<-pedigree[,"div.obs"]
			dtime<-pedigree[,"delta.t"]

			## Reading in parameter estimates
			est <-final.1
			alpha <-as.numeric(est[1,1])
	    beta<-as.numeric(est[1,2])
		  weight<-as.numeric(est[1,3])
	    intercept<-as.numeric(est[1,4])

			## Reading initial state vector
			settings<-info.out
			PrMM<-p0mm<-as.numeric(as.character(settings[1,2]))
			PrUM<-p0um<-as.numeric(as.character(settings[2,2]))
			PrUU<-p0uu<-as.numeric(as.character(settings[3,2]))
			time1<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time2<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time.out<-expand.grid(time1,time2)
			#time0<- rep(min(pedigree[,1]), nrow(time.out))
			time0<- rep(0, nrow(time.out))
			pedigree.new<-as.matrix(cbind(time0,time.out))
			pedigree.new<-cbind(pedigree.new, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			pedigree.new<-pedigree.new[!duplicated(pedigree.new[,4]), ]
			pedigree.new<-pedigree.new[,1:3]

			## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			  svGzero   <- c(PrUU, weight*PrMM, (1-weight)*PrMM)

							alphafinal<-alpha
							betfinal<-beta
							interceptfinal<-intercept

							## Defining the generation (or transition) matrix for the mitotic case
							Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha,alpha^2, bet*(1-alpha), (1-alpha)*(1-bet)+alpha*bet, alpha*(1-bet),
							bet^2, 2*(1-bet)*bet, (1-bet)^2),nrow=3, byrow=TRUE)

							## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
							Dt1t2<-NULL

								for (p in seq_len(NROW(pedigree.new)))
								{

									## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
									svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,1]))
									svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))

									## Conditional divergences
									dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
												 svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

									dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
									             svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

									dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
									             svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

									## Total (weighted) divergence
									Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

								}

			pedigree.new<-cbind(pedigree.new, Dt1t2+interceptfinal, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
			pedigree.new<-pedigree.new[order(pedigree.new[,5]),]


	model<-"ABneutralSOMA.R"

	abfree.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
	names(abfree.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")



	## Ouputting result datasets
	dput(abfree.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))
	return(abfree.out)


} #End of function


ABneutral<-function(pedigree.data, p0uu, eqp, eqp.weight, Nstarts, out.dir, out.name)
{

##### Defining the divergence function
    divergence <- function(pedigree, p0mm, p0um, p0uu, param)
      {

    ## Initializing parameters
    PrMM <- p0mm
    PrUM <- p0um
    PrUU <- p0uu
    alpha <- param[1]
    bet <- param[2]
    weight <- param[3]


	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)



	## Defining the generation (or transition) matrix
	  Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha, alpha^2,
							1/4*(bet + 1 - alpha)^2, 1/2*(bet + 1 - alpha)*(alpha + 1 - bet), 1/4*(alpha + 1 - bet)^2,
							bet^2, 2*(1-bet)*bet, (1-bet)^2), nrow=3, byrow=TRUE)


	## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
	  Dt1t2<-NULL

		  for (p in seq_len(NROW(pedigree)))
		  {

			## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
			svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
			svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
			svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
			svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
			svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

			## Conditional divergences
			dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
								svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

			dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
								svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

			dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
								svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

			## Total (weighted) divergence
			Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


		  }

	  # Pr(UU) at equilibrium given alpha and beta
	  puuinf.est<- (bet* ((1-bet)^2 - (1-alpha)^2 -1))/((alpha + bet)*((alpha + bet -1)^2 - 2))
	  divout<-list(puuinf.est, Dt1t2)

	  return(divout)

	}

		###### Defining the Least Square function to be minimized
		###### Note the equilibrium constraint, which can be made as small as desired.

		LSE_intercept<-function(param_int)
		{
		  sum((pedigree[,4] - param_int[4] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[2]])^2) +
		    eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[1]]- eqp)^2)
		}



		###### Calculating the initial proportions
		###### We always assume that:
		# 1. p0mm is larger than actually observed. This means if p0um is available from measurements,
		#    we will just add it to p0mm.
		# 2. As a consequence of (1.) we also assume that p0um = 0.

		p0uu<-p0uu
		p0mm<-1-p0uu
		p0um<-0


		if(is.null(p0uu ==TRUE | is.null(eqp)==TRUE))
		{stop("Both eqp value AND p0uu have to be supplied")}

		if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1)
		{stop("The initial state probabilities don't sum to 1")}




##### Initializing
    optim.method<-"Nelder-Mead"
    final<-NULL
    counter<-0
    opt.out<-NULL
    pedigree<-pedigree.data



		for (s in seq_len(Nstarts))
		{
			## Draw random starting values
			alpha.start  <-10^(runif(1, log10(10^-9), log10(10^-2)))
			beta.start   <-10^(runif(1, log10(10^-9), log10(10^-2)))
	    weight.start <-runif(1,0,0.1)
	    intercept.start <-runif(1,0,max(pedigree[,4]))
			param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)

			## Initializing
			counter<-counter+1

			message("Progress: ", counter/Nstarts, "\n")

						opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
						alphafinal<-opt.out[1]
						betfinal<-opt.out[2]
						PrMMinf <- (alphafinal* ((1-alphafinal)^2 - (1-betfinal)^2 -1))/((alphafinal + betfinal)*((alphafinal + betfinal -1)^2 - 2))
						PrUMinf <- (4*alphafinal*betfinal*(alphafinal + betfinal -2))/((alphafinal + betfinal)*((alphafinal + betfinal -1)^2 -2))
						PrUUinf <- (betfinal* ((1-betfinal)^2 - (1-alphafinal)^2 -1))/((alphafinal + betfinal)*((alphafinal + betfinal -1)^2 - 2))
						opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
						final[[s]]<- opt.out
		} # End of Nstarts loop

    final <- do.call("rbind", final)
    colnames(final)[1:4]<-c("alpha", "beta", "weight", "intercept")
    colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")




##### Calculating the least square of the first part of the minimized function
	 lsqpart<-NULL

	 for (l in seq_len(NROW(final)))
	 {
			  PrMM <- p0mm
			  PrUM <- p0um
	      PrUU <- p0uu
			  alpha  <- final[l, "alpha"]
			  bet    <- final[l, "beta"]
			  weight <- final[l, "weight"]
			  intercept<-final[l,"intercept"]


			## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)



			  ## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha, alpha^2,
									1/4*(bet + 1 - alpha)^2, 1/2*(bet + 1 - alpha)*(alpha + 1 - bet), 1/4*(alpha + 1 - bet)^2,
									bet^2, 2*(1-bet)*bet, (1-bet)^2), nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

					## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

					## Conditional divergences
					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

					## Total (weighted) divergence
					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


				  }


			 ## Calculating the least square part
			 lsqpart[l]<-sum((pedigree[,4] - intercept - Dt1t2)^2)
		}

	 final<-cbind(final, lsqpart)
	 colnames(final)[ncol(final)]<-c("value.part")
	 final<-final[order(final[,"value"]),]
	 index.1<- which(final["alpha"] > 0 & final["beta"] > 0 & final["intercept"] > 0 & final[,"weight"] > 0 & final["convcode"] == 0)
	 #index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["intercept"] > 0)
	 index.2<-setdiff(seq_len(NROW(final)) , index.1)
	 final.1<-final[index.1,]
	 final.2<-final[index.2,]



##### Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
	 #cat("Caution: Calculating predicted divergence based on lowest LSQ model: check the biology!", "\n")
	 PrMM <- p0mm
	 PrUM <- p0um
	 PrUU <- p0uu
	 alpha  <- final.1[1, "alpha"]
	 bet    <- final.1[1, "beta"]
	 weight <- final.1[1, "weight"]
	 intercept<-final.1[1,"intercept"]


	 ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	 svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)


			  ## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha, alpha^2,
									1/4*(bet + 1 - alpha)^2, 1/2*(bet + 1 - alpha)*(alpha + 1 - bet), 1/4*(alpha + 1 - bet)^2,
									bet^2, 2*(1-bet)*bet, (1-bet)^2), nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL
			  Residual<-NULL

				  for (p in seq_len(NROW(pedigree)))
				  {

					## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
					svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
					svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
					svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
					svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

					## Conditional divergences
					dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
										svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

					dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
										svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

					dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
										svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

					## Total (weighted) divergence
					Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

				  }

			 ## Calculating the least square part
			 Residual<-(pedigree[,4] - intercept - Dt1t2)



    ##### Augmenting pedigree
    delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
    #pedigree<-cbind(pedigree,delta.t)
    pedigree<-cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
    colnames(pedigree)[c(4,5,6,7)]<-c("div.obs", "delta.t","div.pred", "residual")


    ##### Making info about settings
    info<-c("p0mm", "p0um", "p0uu", "eqp", "eqp.weight", "Nstarts", "optim.method")
    info2<-c(p0mm, p0um, p0uu, eqp, eqp.weight, Nstarts, optim.method)
    info.out<-data.frame(info, info2)
    colnames(info.out)<-c("Para", "Setting")






###### Generating theoretical fit

			## Reading in pedigree
			obs<-pedigree[,"div.obs"]
			dtime<-pedigree[,"delta.t"]

			## Reading in parameter estimates
			est <-final.1
			alpha <-as.numeric(est[1,1])
	    beta<-as.numeric(est[1,2])
		  weight<-as.numeric(est[1,3])
	    intercept<-as.numeric(est[1,4])

			## Reading initial state vector
			settings<-info.out
			PrMM<-p0mm<-as.numeric(as.character(settings[1,2]))
			PrUM<-p0um<-as.numeric(as.character(settings[2,2]))
			PrUU<-p0uu<-as.numeric(as.character(settings[3,2]))
			time1<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time2<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time.out<-expand.grid(time1,time2)
			#time0<- rep(min(pedigree[,1]), nrow(time.out))
			time0<- rep(0, nrow(time.out))
			pedigree.new<-as.matrix(cbind(time0,time.out))
			pedigree.new<-cbind(pedigree.new, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			pedigree.new<-pedigree.new[!duplicated(pedigree.new[,4]), ]
			pedigree.new<-pedigree.new[,1:3]

			## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			svGzero   <- c(PrUU, PrMM*weight, (1-weight)*PrMM)


							alphafinal<-alpha
							betfinal<-beta
							interceptfinal<-intercept

							## Defining the generation (or transition) matrix
							Genmatrix <- matrix(c((1-alphafinal)^2, 2*(1-alphafinal)*alphafinal, alphafinal^2,
												   1/4*(betfinal + 1 - alphafinal)^2, 1/2*(betfinal + 1 - alphafinal)*(alphafinal + 1 - betfinal),
												   1/4*(alphafinal + 1 - betfinal)^2,
												   betfinal^2, 2*(1-betfinal)*betfinal, (1-betfinal)^2), nrow=3, byrow=TRUE)

							## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
							Dt1t2<-NULL

								for (p in seq_len(NROW(pedigree.new)))
								{

									## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
									svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,1]))
									svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))

									## Conditional divergences
									dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
												 svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

									dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
									             svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

									dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
									             svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

									## Total (weighted) divergence
									Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM

								}

			pedigree.new<-cbind(pedigree.new, Dt1t2+interceptfinal, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
			pedigree.new<-pedigree.new[order(pedigree.new[,5]),]

    model<-"ABneutral.R"

    abfree.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
    names(abfree.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")

    ## Ouputting result datasets
    dput(abfree.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))
    return(abfree.out)


} #End of function

