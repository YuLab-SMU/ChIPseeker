##' prepare the promoter regions
##'
##'
##' @title getPromoters
##' @param TxDb TxDb
##' @param upstream upstream from TSS site
##' @param downstream downstream from TSS site
##' @param by one of gene or transcript
##' @return GRanges object
##' @export
getPromoters <- function(TxDb=NULL,
                         upstream=1000,
                         downstream=1000,
                         by = "gene") {
  
  getBioRegion(TxDb = TxDb,
               upstream = upstream,
               downstream = downstream,
               by = by,
               type = "start_site")
}


##' prepare a bioregion of selected feature
##' 
##' this function combined previous functions getPromoters(),getBioRegion(),getGeneBody()
##' https://github.com/GuangchuangYu/ChIPseeker/issues/16
##' https://github.com/GuangchuangYu/ChIPseeker/issues/87
##'
##' @title getBioRegion
##' @param TxDb TxDb
##' @param upstream upstream from start site
##' @param downstream downstream from start site
##' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR'
##' @param type one of "start_site", "end_site", "body"
##' @return GRanges object
##' @import BiocGenerics IRanges GenomicRanges
##' @author Guangchuang Yu, Ming L
##' @export
getBioRegion <- function(TxDb=NULL,
                         upstream=1000,
                         downstream=1000,
                         by="gene",
                         type="start_site"){
  
  by <- match.arg(by, c('gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR'))
  type <- match.arg(type, c("start_site", "end_site", "body"))
  
  TxDb <- loadTxDb(TxDb)
  .ChIPseekerEnv(TxDb)
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  if(type == 'body'){
    if(by %in% c('gene', 'transcript', 'exon', 'intron')){
      label_SS <- paste0("T","SS")
      label_TS <- paste0("T","TS")
      label <- c(label_SS,label_TS)
    }else{
      label_SS <- paste0(by,"_SS")
      label_TS <- paste0(by,"_TS")
      label <- c(label_SS,label_TS)
    }
    
  }else if(type == "start_site"){
    if(by %in% c('gene', 'transcript', 'exon', 'intron')){
      label <- paste0("T","SS")
    }else{
      label <- paste0(by,"_SS") 
    }
    
  }else{
    if(by %in% c('gene', 'transcript', 'exon', 'intron')){
      label <- paste0("T","TS")
    }else{
      label <- paste0(by,"_TS")
    }
  }
  
  
  if(by == 'gene' || by == 'transcript'){
    regions <- getGene(TxDb, by)
  }
  
  if (by == "exon") {
    exonList <- get_exonList(ChIPseekerEnv)
    regions <-  unlist(exonList)
  }
  
  if (by == "intron") {
    intronList <- get_intronList(ChIPseekerEnv)
    regions <- unlist(intronList)
  }
  
  if (by == "3UTR") {
    threeUTRList <- threeUTRsByTranscript(TxDb)
    regions <- unlist(threeUTRList)
  }
  
  if (by == "5UTR") {
    fiveUTRList <- fiveUTRsByTranscript(TxDb)
    regions <- unlist(fiveUTRList)
  }
  
  if(type == "start_site"){
    coordinate<- ifelse(strand(regions) == "+", start(regions), end(regions))
  }else if(type == "end_site"){
    coordinate<- ifelse(strand(regions) == "+", end(regions), start(regions))
  }else{
    ## assign attribute 
    attr(regions, 'type') = type
    attr(regions, 'label') = label
    
    return(regions)
  }
  
  bioRegion <- GRanges(seqnames=seqnames(regions),
                       ranges=IRanges(coordinate-upstream, coordinate+downstream),
                       strand=strand(regions))
  bioRegion <- unique(bioRegion)
  
  ## assign attribute 
  attr(bioRegion, 'type') = type
  attr(bioRegion, 'label') = label
  
  return(bioRegion)
}



##' calculate the tag matrix
##' 
##' 
##' @title getTagMatrix
##'
##' @param peak peak peak file or GRanges object
##' @param upstream the distance of upstream extension
##' @param downstream the distance of downstream extension
##' @param windows a collection of region
##' @param type one of "start_site", "end_site", "body"
##' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR'
##' @param TxDb TxDb
##' @param weightCol column name of weight, default is NULL
##' @param nbin the amount of nbines 
##' @param verbose print message or not
##' @param flip_minor_strand whether flip the orientation of minor strand
##' @return tagMatrix
##' @importFrom ggplot2 rel
##' @export
getTagMatrix <- function(peak, 
                         upstream,
                         downstream, 
                         windows,
                         type,
                         by,
                         TxDb=NULL,
                         weightCol = NULL, 
                         nbin = 800,
                         verbose = TRUE,
                         flip_minor_strand=TRUE){
  
  ## max_region_length decides to use binning method or not
  max_region_length <- getOption("ChIPseeker.max_region_length")
  
  ## check upstream and downstream parameter
  check_upstream_and_downstream(upstream = upstream, downstream = downstream)
  
  if(missingArg(windows)){
    windows <- getBioRegion(TxDb=TxDb,
                            upstream=upstream,
                            downstream=downstream,
                            by=by,
                            type=type)
  }else{
    
    if (! is(windows, "GRanges")) {
      stop("windows should be a GRanges object...")
    }
    
    type <- attr(windows, 'type')
    by <- attr(windows, 'by')
    
  }
  
  if(type != 'body'){
    if(inherits(upstream, 'rel') || is.null(upstream)){
      stop("upstream and downstream for site region should be actual number...")
    }
  }
  
  if (verbose) {
    cat(">> preparing ",type," regions"," by ",by,"... ",
        format(Sys.time(), "%Y-%m-%d %X"), "\n",sep = "")
  }
  
  
  if(type == "body" || width(windows[1]) > max_region_length){
    
    if (verbose) {
      cat(">> preparing tag matrix by binning... ",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    
    tagMatrix <- getTagMatrix.binning.internal(peak = peak, 
                                               weightCol = weightCol, 
                                               windows = windows, 
                                               nbin = nbin,
                                               upstream = upstream,
                                               downstream = downstream)
    
    attr(tagMatrix, "is.binning") <- T
    
  }else{
    
    if (verbose) {
      cat(">> preparing tag matrix... ",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    
    tagMatrix <- getTagMatrix.internal(peak=peak, 
                                       weightCol=weightCol, 
                                       windows=windows, 
                                       flip_minor_strand=flip_minor_strand)
    
    attr(tagMatrix, "is.binning") <- F
  }
  
  ## assign attribute 
  attr(tagMatrix, 'upstream') = upstream
  attr(tagMatrix, 'downstream') = downstream
  attr(tagMatrix, 'type') = attr(windows, 'type')
  attr(tagMatrix, 'label') = attr(windows, 'label')
  
  return(tagMatrix)
}


##' calculate the tag matrix
##'
##'
##' @title getTagMatrix.internal
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight, default is NULL
##' @param windows a collection of region with equal size, eg. promoter region.
##' @param flip_minor_strand whether flip the orientation of minor strand
##' @return tagMatrix
##' @import BiocGenerics S4Vectors IRanges GenomeInfoDb GenomicRanges
##' @author G Yu
getTagMatrix.internal <- function(peak, 
                                  weightCol=NULL, 
                                  windows, 
                                  flip_minor_strand=TRUE) {
  peak.gr <- loadPeak(peak)
  
  if (! is(windows, "GRanges")) {
    stop("windows should be a GRanges object...")
  }
  if (length(unique(width(windows))) != 1) {
    stop("width of windows should be equal...")
  }
  
  ## if (!exists("ChIPseekerEnv", envir = .GlobalEnv)) {
  ##     assign("ChIPseekerEnv", new.env(), .GlobalEnv)
  ## }
  ## ChIPseekerEnv <- get("ChIPseekerEnv", envir = .GlobalEnv)
  
  ## if (exists("peak", envir=ChIPseekerEnv, inherits=FALSE) &&
  ##     exists("promoters", envir=ChIPseekerEnv, inherits=FALSE) &&
  ##     exists("weightCol", envir=ChIPseekerEnv, inherits=FALSE) &&
  ##     exists("tagMatrix", envir=ChIPseekerEnv, inherits=FALSE) ) {
  
  ##     pp <- get("peak", envir=ChIPseekerEnv)
  ##     promoters <- get("promoters", envir=ChIPseekerEnv)
  ##     w <- get("weightCol", envir=ChIPseekerEnv)
  
  ##     if (all(pp == peak)) {
  ##         if (all(windows == promoters)) {
  ##             if ( (is.null(w) && is.null(weightCol)) ||
  ##                 (!is.null(w) && !is.null(weightCol) && w == weightCol)) {
  ##                 tagMatrix <- get("tagMatrix", envir=ChIPseekerEnv)
  ##                 return(tagMatrix)
  ##             } else {
  ##                 assign("weightCol", weightCol, envir=ChIPseekerEnv)
  ##             }
  ##         } else {
  ##             assign("promoters", windows)
  ##             ## make sure it is not conflict with getPromoters
  ##             if ( exists("upstream", envir=ChIPseekerEnv, inherits=FALSE))
  ##                 rm("upstream", envir=ChIPseekerEnv)
  ##         }
  ##     } else {
  ##         assign("peak", peak, envir=ChIPseekerEnv)
  ##     }
  
  ## }
  
  ## if ( !exists("peak", envir=ChIPseekerEnv, inherits=FALSE)) {
  ##     assign("peak", peak, envir=ChIPseekerEnv)
  ## }
  
  ## if ( !exists("promoters", envir=ChIPseekerEnv, inherits=FALSE)) {
  ##     assign("promoters", windows, envir=ChIPseekerEnv)
  ## }
  
  ## if (!exists("weightCol", envir=ChIPseekerEnv, inherits=FALSE)) {
  ##     assign("weightCol", weightCol, envir=ChIPseekerEnv)
  ## }
  if (is.null(weightCol)) {
    peak.cov <- coverage(peak.gr)
  } else {
    weight <- mcols(peak.gr)[[weightCol]]
    peak.cov <- coverage(peak.gr, weight=weight)
  }
  cov.len <- elementNROWS(peak.cov)
  cov.width <- GRanges(seqnames=names(cov.len),
                       IRanges(start=rep(1, length(cov.len)),
                               end=cov.len))
  windows <- subsetByOverlaps(windows, cov.width,
                              type="within", ignore.strand=TRUE)
  
  chr.idx <- intersect(names(peak.cov),
                       unique(as.character(seqnames(windows))))
  
  peakView <- Views(peak.cov[chr.idx], as(windows, "IntegerRangesList")[chr.idx])
  tagMatrixList <- lapply(peakView, function(x) t(viewApply(x, as.vector)))
  tagMatrix <- do.call("rbind", tagMatrixList)
  
  ## get the index of windows, that are reorganized by as(windows, "IntegerRangesList")
  idx.list <- split(1:length(windows),  as.factor(seqnames(windows)))
  idx <- do.call("c", idx.list)
  
  rownames(tagMatrix) <- idx
  tagMatrix <- tagMatrix[order(idx),]
  
  ## minus strand
  if (flip_minor_strand) {
    ## should set to FALSE if upstream is not equal to downstream
    ## can set to TRUE if e.g. 3k-TSS-3k
    ## should set to FALSE if e.g. 3k-TSS-100
    minus.idx <- which(as.character(strand(windows)) == "-")
    tagMatrix[minus.idx,] <- tagMatrix[minus.idx, ncol(tagMatrix):1]
  }
  
  tagMatrix <- tagMatrix[rowSums(tagMatrix)!=0,]
  ## assign("tagMatrix", tagMatrix, envir=ChIPseekerEnv)
  return(tagMatrix)
}


##' calculate the tagMatrix by binning
##' the idea was derived from the function of deeptools
##' https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html 
##' 
##' @title getTagMatrix.binning.internal
##' @param peak peak peak file or GRanges object
##' @param weightCol weightCol column name of weight, default is NULL
##' @param windows windows a collection of region with equal or not equal size, eg. promoter region, gene region.
##' @param nbin the amount of nbines needed to be splited and it should not be more than min_body_length
##' @param upstream rel object, NULL or actual number
##' @param downstream rel object, NULL or actual number
##' @import BiocGenerics S4Vectors IRanges GenomeInfoDb GenomicRanges 
##' @importFrom ggplot2 rel
##' @return tagMatrix 
getTagMatrix.binning.internal <- function(peak, 
                                          weightCol = NULL, 
                                          windows, 
                                          nbin = 800,
                                          upstream = NULL,
                                          downstream = NULL){
  
  min_body_length <- filter_length <- nbin
  peak.gr <- loadPeak(peak)
  
  if (!is(windows, "GRanges")) {
    stop("windows should be a GRanges object...")
  }
  
  if (is.null(weightCol)) {
    peak.cov <- coverage(peak.gr)
  } else {
    weight <- mcols(peak.gr)[[weightCol]]
    peak.cov <- coverage(peak.gr, weight=weight)
  }
  
  
  cov.len <- elementNROWS(peak.cov)
  cov.width <- GRanges(seqnames=names(cov.len),
                       IRanges(start=rep(1, length(cov.len)),
                               end=cov.len))
  
  windows <- subsetByOverlaps(windows, 
                              cov.width,
                              type="within", 
                              ignore.strand=TRUE)
  
  ## extend the windows by rel object
  if(inherits(upstream, 'rel')){
    
    windows1 <- windows
    start(windows1) <- start(windows1) - floor(width(windows)*as.numeric(upstream))
    end(windows1) <- end(windows1) + floor(width(windows)*as.numeric(downstream))
    windows <- windows1
    nbin <- floor(nbin*(1+as.numeric(downstream)+as.numeric(upstream)))
    min_body_length <- min_body_length*(1+as.numeric(upstream)+as.numeric(downstream))
    
    cat(">> preparing matrix with extension from (",attr(windows,'label')[1],"-",
        100*as.numeric(upstream),"%)~(",attr(windows,'label')[2],"+",
        100*as.numeric(downstream),"%)... ",
        format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
  }
  
  ## do not extend
  if(is.null(upstream)){
    if(attr(windows, 'type') == 'body'){
      cat(">> preparing matrix for ",attr(windows, 'type')," region with no flank extension... ",
          format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
    }else{
      cat(">> preparing matrix for ",attr(windows,'type')," region... ",
          format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
    }
  }
  
  ## extend the windows by actual number 
  if(!is.null(upstream) && !inherits(upstream, 'rel') && attr(windows, 'type')== 'body'){
    windows1 <- windows
    start(windows1) <- start(windows1) - upstream
    end(windows1) <- end(windows1) + downstream
    windows <- windows1
    upstreamPer <- floor(upstream/1000)*0.1
    downstreamPer <- floor(downstream/1000)*0.1
    nbin <- floor(nbin*(1+upstreamPer+downstreamPer))
    min_body_length <- min_body_length+upstream+downstream
    
    cat(">> preparing matrix with flank extension from (",attr(windows,'label')[1],"-",
        upstream,"bp)~(",attr(windows,'label')[2],"+",downstream,"bp)... ",
        format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
  }
  
  chr.idx <- intersect(names(peak.cov),
                       unique(as.character(seqnames(windows))))
  
  peakView <- Views(peak.cov[chr.idx], as(windows, "IntegerRangesList")[chr.idx])
  
  ## remove the gene that has no binding proteins
  peakView <- lapply(peakView, function(x) x <- x[viewSums(x)!=0])
  
  tagMatrixList <- lapply(peakView, function(x) viewApply(x, as.vector))
  
  if(!attr(windows, 'type') == 'body'){
    
    tagMatrixList <- lapply(tagMatrixList, function(x) t(x))
    
    # to remove the chromosome that has no binding protein
    tagMatrixList <- tagMatrixList[vapply(tagMatrixList, function(x) length(x)>0, FUN.VALUE = logical(1))]
    suppressWarnings(tagMatrixList <- do.call("rbind", tagMatrixList))
    
    ## create a matrix to receive binning result                 
    tagMatrix  <- matrix(nrow = dim(tagMatrixList)[1],ncol = nbin)
    
    ## this circulation is to deal with different gene             
    for (i in 1:dim(tagMatrixList)[1]) {
      
      ## seq is the distance between different nbines
      seq <- floor(length(tagMatrixList[i,])/nbin)
      
      ## cursor record the position of calculation
      cursor <- 1
      
      ## the third circulation is to calculate the bingding strength
      ## it has two parts
      ## the first part is to for the nbin(1:nbin-1)
      ## because the seq is not derived from exact division
      ## the second part is to compensate the loss of non-exact-division
      
      ## this the first part for 1:(nbin-1)
      for (j in 1:(nbin-1)) {
        
        read <- 0
        
        for (k in cursor:(cursor+seq-1)) {
          read <- read + tagMatrixList[i,k]
        }
        
        tagMatrix [i,j] <- read/seq
        
        cursor <- cursor+seq
      }
      
      ## this the second part to to compensate the loss of non-exact-division
      read <- 0
      for (k in cursor:length(tagMatrixList[i,])) {
        read <- read+tagMatrixList[i,k]
      }
      
      tagMatrix [i,nbin] <- read/(length(tagMatrixList[i,])-cursor)
    }
    
  }else{
    
    ## extend genebody by atual number
    if(!is.null(upstream) & !inherits(upstream, 'rel')){
      
      tagMatrixList <- tagMatrixList[vapply(tagMatrixList, function(x) length(x)>0, FUN.VALUE = logical(1))]
      
      ## count the amount before filtering
      pre_amount <- 0
      for(i in 1:length(tagMatrixList)){
        pre_amount <- pre_amount+length(tagMatrixList[[i]])
      }
      
      tagMatrixList <- lapply(tagMatrixList, function(x) x[vapply(x, function(y) length(y)>min_body_length,FUN.VALUE = logical(1))])
      
      ## count the amount after filtering
      amount <- 0
      for(i in 1:length(tagMatrixList)){
        amount <- amount+length(tagMatrixList[[i]])
      }
      
      cat(">> ",pre_amount-amount," peaks(",100*((pre_amount-amount)/pre_amount),
          "%), having lengths smaller than ",filter_length,"bp, are filtered... ",
          format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
      
      suppressWarnings(tagMatrixList <- do.call("rbind",tagMatrixList))
      
      tagMatrix  <- matrix(nrow = length(tagMatrixList),ncol = nbin)
      
      upstreamnbin <- floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer)))
      bodynbin <- floor(nbin*(1/(1+upstreamPer+downstreamPer)))
      downstreamnbin <- floor(nbin*(downstreamPer/(1+upstreamPer+downstreamPer)))
      
      # count the upstream 
      for (i in 1:length(tagMatrixList)) {
        
        seq <- floor(upstream/upstreamnbin)
        cursor <- 1
        
        for (j in 1:(upstreamnbin-1)) {
          
          read <- 0
          
          for (k in cursor:(cursor+seq-1)) {
            read <- read + tagMatrixList[[i]][k]
          }
          
          tagMatrix [i,j] <- read/seq
          
          cursor <- cursor+seq
        }
        
        
        read <- 0
        for (k in cursor:upstream) {
          read <- read+tagMatrixList[[i]][k]
        }
        
        tagMatrix [i,upstreamnbin] <- read/(upstream-cursor)
      }
      
      
      ## count genebody
      for (i in 1:length(tagMatrixList)) {
        
        seq <- floor((length(tagMatrixList[[i]])-upstream-downstream)/bodynbin)
        cursor <- upstream+1
        
        for (j in (upstreamnbin+1):(upstreamnbin+bodynbin-1)) {
          
          read <- 0
          
          for (k in cursor:(cursor+seq-1)) {
            read <- read + tagMatrixList[[i]][k]
          }
          
          tagMatrix [i,j] <- read/seq
          
          cursor <- cursor+seq
        }
        
        read <- 0
        for (k in cursor:(length(tagMatrixList[[i]])-downstream)) {
          read <- read+tagMatrixList[[i]][k]
        }
        
        tagMatrix [i,bodynbin+upstreamnbin] <- read/(length(tagMatrixList[[i]])-downstream-cursor)
      }
      
      
      ## count downstream
      for (i in 1:length(tagMatrixList)) {
        
        seq <- floor(downstream/downstreamnbin)
        cursor <- length(tagMatrixList[[i]])-downstream+1
        
        for (j in (upstreamnbin+bodynbin+1):(nbin-1)) {
          
          read <- 0
          
          for (k in cursor:(cursor+seq-1)) {
            read <- read + tagMatrixList[[i]][k]
          }
          
          tagMatrix [i,j] <- read/seq
          
          cursor <- cursor+seq
        }
        
        read <- 0
        for (k in cursor:length(tagMatrixList[[i]])) {
          read <- read+tagMatrixList[[i]][k]
        }
        
        tagMatrix [i,nbin] <- read/(length(tagMatrixList[[i]])-cursor)
      }
      
    }else{
      
      tagMatrixList <- tagMatrixList[vapply(tagMatrixList, function(x) length(x)>0, FUN.VALUE = logical(1))]
      
      ## count the amount before filtering
      pre_amount <- 0
      for(i in 1:length(tagMatrixList)){
        pre_amount <- pre_amount+length(tagMatrixList[[i]])
      }
      
      tagMatrixList <- lapply(tagMatrixList, function(x) x[vapply(x, function(y) length(y)>min_body_length, FUN.VALUE = logical(1))])
      
      ## count the amount after filtering
      amount <- 0
      for(i in 1:length(tagMatrixList)){
        amount <- amount+length(tagMatrixList[[i]])
      }
      
      cat(">> ",pre_amount-amount," peaks(",100*((pre_amount-amount)/pre_amount),
          "%), having lengths smaller than ",filter_length,"bp, are filtered... ",
          format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
      
      suppressWarnings(tagMatrixList <- do.call("rbind",tagMatrixList))
      
      tagMatrix  <- matrix(nrow = length(tagMatrixList),ncol = nbin)
      
      for (i in 1:length(tagMatrixList)) {
        
        seq <- floor(length(tagMatrixList[[i]])/nbin)
        cursor <- 1
        
        for (j in 1:(nbin-1)) {
          
          read <- 0
          
          for (k in cursor:(cursor+seq-1)) {
            read <- read + tagMatrixList[[i]][k]
          }
          
          tagMatrix [i,j] <- read/seq
          
          cursor <- cursor+seq
        }
        
        read <- 0
        for (k in cursor:length(tagMatrixList[[i]])) {
          read <- read+tagMatrixList[[i]][k]
        }
        
        tagMatrix [i,nbin] <- read/(length(tagMatrixList[[i]])-cursor)
      }
    }
    
  }
  
  return(tagMatrix)
}