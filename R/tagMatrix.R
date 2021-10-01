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
##' @import BiocGenerics IRanges GenomicRanges
##' @importFrom GenomicFeatures transcriptsBy
getPromoters <- function(TxDb=NULL,
                         upstream=1000,
                         downstream=1000,
                         by = "gene") {
  
  by <- match.arg(by, c("gene", "transcript"))
  
  TxDb <- loadTxDb(TxDb)
  .ChIPseekerEnv(TxDb)
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  if ( exists("upstream", envir=ChIPseekerEnv, inherits=FALSE) &&
       exists("downstream", envir=ChIPseekerEnv, inherits=FALSE) ) {
    us <- get("upstream", envir=ChIPseekerEnv)
    ds <- get("downstream", envir=ChIPseekerEnv)
    if (us == upstream && ds == downstream &&
        exists("promoters", envir=ChIPseekerEnv, inherits=FALSE) ){
      promoters <- get("promoters", envir=ChIPseekerEnv)
      
      ## assign attribute 
      attr(promoters, 'type') = 'TSS'
      
      return(promoters)
    }
  }
  
  Transcripts <- getGene(TxDb, by)
  ## get start position based on strand
  tss <- ifelse(strand(Transcripts) == "+", start(Transcripts), end(Transcripts))
  promoters <- GRanges(seqnames=seqnames(Transcripts),
                       ranges=IRanges(tss-upstream, tss+downstream),
                       strand=strand(Transcripts))
  promoters <- unique(promoters)
  
  assign("promoters", promoters, envir=ChIPseekerEnv)
  assign("upstream", upstream, envir=ChIPseekerEnv)
  assign("downstream", downstream, envir=ChIPseekerEnv)
  
  ## assign attribute 
  attr(promoters, 'type') = 'TSS'
  
  return(promoters)
}


##' get the tts region 
##' 
##' 
##' @title getTTSRegion
##' @param TxDb TxDb
##' @param upstream upstream from TSS site
##' @param downstream downstream from TSS site
##' @param by one of gene or transcript
##' @return GRanges object
##' @import BiocGenerics IRanges GenomicRanges
##' @importFrom GenomicFeatures transcriptsBy
##' @export
##' https://github.com/GuangchuangYu/ChIPseeker/issues/87
getTTSRegion <- function(TxDb=NULL,
                         upstream=1000,
                         downstream=1000,
                         by = "gene"){
  
  by <- match.arg(by, c("gene", "transcript"))
  
  TxDb <- loadTxDb(TxDb)
  .ChIPseekerEnv(TxDb)
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  if ( exists("upstream", envir=ChIPseekerEnv, inherits=FALSE) &&
       exists("downstream", envir=ChIPseekerEnv, inherits=FALSE) ) {
    us <- get("upstream", envir=ChIPseekerEnv)
    ds <- get("downstream", envir=ChIPseekerEnv)
    if (us == upstream && ds == downstream &&
        exists("TTS", envir=ChIPseekerEnv, inherits=FALSE) ){
      TTS_region <- get("TTS", envir=ChIPseekerEnv)
      
      ## assign attribute 
      attr(TTS_region, 'type') = 'TTS'
      
      return(TTS_region)
    }
  }
  
  Transcripts <- getGene(TxDb, by)
  ## get end position based on strand
  tts <- ifelse(strand(Transcripts) == "+", end(Transcripts), start(Transcripts))
  TTS_region <- GRanges(seqnames=seqnames(Transcripts),
                        ranges=IRanges(tts-upstream, tts+downstream),
                        strand=strand(Transcripts))
  TTS_region <- unique(TTS_region)
  
  assign("TTS", TTS_region, envir=ChIPseekerEnv)
  assign("upstream", upstream, envir=ChIPseekerEnv)
  assign("downstream", downstream, envir=ChIPseekerEnv)
  
  ## assign attribute 
  attr(TTS_region, 'type') = 'TTS'
  
  return(TTS_region)
  
}

##' prepare a region center on end site of selected feature
##'
##' 
##' @title getEndRegion
##' @param TxDb TxDb
##' @param upstream upstream from start site
##' @param downstream downstream from start site
##' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR'
##' @return GRanges object
##' @import BiocGenerics IRanges GenomicRanges
##' @export
##' @author Guangchuang Yu
##  https://github.com/GuangchuangYu/ChIPseeker/issues/87
getEndRegion <- function(TxDb=NULL,
                         upstream=1000,
                         downstream=1000,
                         by="gene") {
  
  by <- match.arg(by, c("gene", "transcript", "exon", "intron", "3UTR", "5UTR"))
  
  if (by %in% c("gene", "transcript")) {
    return(getTTSRegion(TxDb, upstream, downstream, by))
  }
  
  TxDb <- loadTxDb(TxDb)
  .ChIPseekerEnv(TxDb)
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  if (by == "exon") {
    exonList <- get_exonList(ChIPseekerEnv)
    regions <-  unlist(exonList)
    
    ## assign attribute 
    attr(regions, 'type') = 'exon_TS'
  }
  
  if (by == "intron") {
    intronList <- get_intronList(ChIPseekerEnv)
    regions <- unlist(intronList)
    
    ## assign attribute 
    attr(regions, 'type') = 'intron_TS'
  }
  
  if (by == "3UTR") {
    threeUTRList <- threeUTRsByTranscript(TxDb)
    regions <- unlist(threeUTRList)
    
    ## assign attribute 
    attr(regions, 'type') = '3UTR_TS'
  }
  
  if (by == "5UTR") {
    fiveUTRList <- fiveUTRsByTranscript(TxDb)
    regions <- unlist(fiveUTRList)
    
    ## assign attribute 
    attr(regions, 'type') = '5UTR_TS'
  }
  
  end_site <- ifelse(strand(regions) == "+", end(regions), start(regions))
  
  endRegion <- GRanges(seqnames=seqnames(regions),
                       ranges=IRanges(end_site-upstream, end_site+downstream),
                       strand=strand(regions))
  endRegion <- unique(endRegion)
  
  ## assign attribute 
  attr(endRegion, 'type') = attr(regions, 'type')
  
  return(endRegion)
}

##' prepare a region center on start site of selected feature
##'
##' 
##' @title getBioRegion
##' @param TxDb TxDb
##' @param upstream upstream from start site
##' @param downstream downstream from start site
##' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR'
##' @return GRanges object
##' @import BiocGenerics IRanges GenomicRanges
##' @export
##' @author Guangchuang Yu
##  https://github.com/GuangchuangYu/ChIPseeker/issues/16
getBioRegion <- function(TxDb=NULL,
                         upstream=1000,
                         downstream=1000,
                         by="gene") {
  
  by <- match.arg(by, c("gene", "transcript", "exon", "intron", "3UTR", "5UTR"))
  
  if (by %in% c("gene", "transcript")) {
    return(getPromoters(TxDb, upstream, downstream, by))
  }
  
  TxDb <- loadTxDb(TxDb)
  .ChIPseekerEnv(TxDb)
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  if (by == "exon") {
    exonList <- get_exonList(ChIPseekerEnv)
    regions <-  unlist(exonList)
    
    ## assign attribute 
    attr(regions, 'type') = 'exon_SS'
  }
  
  if (by == "intron") {
    intronList <- get_intronList(ChIPseekerEnv)
    regions <- unlist(intronList)
    
    ## assign attribute 
    attr(regions, 'type') = 'intron_SS'
  }
  
  if (by == "3UTR") {
    threeUTRList <- threeUTRsByTranscript(TxDb)
    regions <- unlist(threeUTRList)
    
    ## assign attribute 
    attr(regions, 'type') = '3UTR_SS'
  }
  
  if (by == "5UTR") {
    fiveUTRList <- fiveUTRsByTranscript(TxDb)
    regions <- unlist(fiveUTRList)
    
    ## assign attribute 
    attr(regions, 'type') = '5UTR_SS'
  }
  
  start_site <- ifelse(strand(regions) == "+", start(regions), end(regions))
  
  bioRegion <- GRanges(seqnames=seqnames(regions),
                       ranges=IRanges(start_site-upstream, start_site+downstream),
                       strand=strand(regions))
  bioRegion <- unique(bioRegion)
  
  ## assign attribute 
  attr(bioRegion, 'type') = attr(regions, 'type')
  
  return(bioRegion)
}



##' calculate the tag matrix
##'
##'
##' @title getTagMatrix
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight, default is NULL
##' @param windows a collection of region with equal size, eg. promoter region.
##' @param flip_minor_strand whether flip the orientation of minor strand
##' @return tagMatrix
##' @export
##' @import BiocGenerics S4Vectors IRanges GenomeInfoDb GenomicRanges
getTagMatrix <- function(peak, weightCol=NULL, windows, flip_minor_strand=TRUE) {
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


##' prepare the genebody region
##' 
##' 
##' Title getGeneBody
##' @param TxDb TxDb
##' @param type by one of "genes", "exon", "intron"
##' @return GRanges object
##' @import BiocGenerics IRanges GenomicRanges
##' @export
getGeneBody <- function(TxDb=NULL,
                        type="genes"){
  
  type <- match.arg(type, c("genes", "exon", "intron"))
  
  TxDb <- loadTxDb(TxDb)
  .ChIPseekerEnv(TxDb)
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  
  ## select the type of gene to be windows
  if( type == "intron"){
    genebody <- get_intronList(ChIPseekerEnv)
  }else if( type == "exon"){
    genebody <- get_exonList(ChIPseekerEnv)
  }else if(type == "genes"){
    genebody <- getGene(TxDb)
  }else {
    stop("type must be genes, exon or intron..")
  }
  
  ## assign attribute 
  attr(genebody, 'type') = 'genebody'
  
  return(genebody)
}



##' calculate the BioRegionMatrix by binning
##' 
##' 
##' Title getBioRegionMatrix
##'
##' @param peak peak peak file or GRanges object
##' @param weightCol weightCol column name of weight, default is NULL
##' @param windows windows a collection of region with equal or not equal size, eg. promoter region, gene region.
##' @param box the amount of boxes needed to be splited and it should not be more than min_body_length
##' @param min_body_length the minimum length that each gene region should be
##' @param upstream rel object reflects the percentage of flank extension, e.g rel(0.2)
##'                 integer reflects the actual length of flank extension
##'                 default(NULL) reflects the gene body with no extension 
##' @param downstream rel object reflects the percentage of flank extension, e.g rel(0.2)
##'                   integer reflects the actual length of flank extension
##'                   default(NULL) reflects the gene body with no extension
##' @import BiocGenerics S4Vectors IRanges GenomeInfoDb GenomicRanges 
##' @importFrom ggplot2 rel
##' @return bioregionmatrix
##' @export
getBioRegionMatrix <- function(peak, 
                               weightCol = NULL, 
                               windows, 
                               box = 800,
                               min_body_length = 1000,
                               upstream = NULL,
                               downstream = NULL,
                               ...){
  
  ## the idea was derived from the function of deeptools
  ## (https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html)  
  peak.gr <- loadPeak(peak)
  
  if (! is(windows, "GRanges")) {
    stop("windows should be a GRanges object...")
  }
  
  ## check upstream and downstream parameter
  check_upstream_and_downstream(upstream = upstream, downstream = downstream)
  
  ## if getting TSS region, upstream and downstream parameter should only be NULL
  if((attr(windows, 'type') == 'start_region') & (!is.null(upstream) | !is.null(downstream))){
    stop("windows is a GRanges object with equal length(e.g promoter region).\n",
         "If you want to extend it, try to extend it when getting windows.")
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
  
  ## extend the windows by upstream and downstream parameter
  if(inherits(upstream, 'rel') || inherits(downstream, 'rel')){
    
    windows1 <- windows
    start(windows1) <- start(windows1) - floor(width(windows)*as.numeric(upstream))
    end(windows1) <- end(windows1) + floor(width(windows)*as.numeric(downstream))
    windows <- windows1
    box <- floor(box*(1+as.numeric(downstream)+as.numeric(upstream)))
    
    cat(">> preparing matrix with upstream and downstream flank extension from ",
        "(TSS-",100*as.numeric(upstream),"%)~(TES+",100*as.numeric(downstream),"%)...\t",
        format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
  }
  
  if(is.null(upstream) || is.null(downstream)){
    if(attr(windows, 'type') == 'genebody'){
      cat(">> preparing matrix for bioregion with no flank extension...\t",
          format(Sys.time(), "%Y-%m-%d %X"),"\n")
    }else{
      cat(">> preparing matrix for start_site_region...\t",
          format(Sys.time(), "%Y-%m-%d %X"),"\n")
    }
  }
  
  if(!is.null(upstream) & !inherits(upstream, 'rel')){
    windows1 <- windows
    start(windows1) <- start(windows1) - upstream
    end(windows1) <- end(windows1) + downstream
    windows <- windows1
    upstreamPer <- floor(upstream/1000)*0.1
    downstreamPer <- floor(downstream/1000)*0.1
    box <- floor(box*(1+upstreamPer+downstreamPer))
    
    cat(">> preparing matrix with upstream and downstream flank extension from ",
        "(TSS-",upstream,"bp)~(TES+",downstream,"bp)...\t",
        format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
  }
  
  chr.idx <- intersect(names(peak.cov),
                       unique(as.character(seqnames(windows))))
  
  peakView <- Views(peak.cov[chr.idx], as(windows, "IntegerRangesList")[chr.idx])
  
  ## remove the gene that has no binding proteins
  peakView <- lapply(peakView, function(x) x <- x[viewSums(x)!=0])
  
  bioregionList <- lapply(peakView, function(x) viewApply(x, as.vector))
  
  ## the "if" judge statement is to be compatible with 
  ## windows that has the equal size, like promoters(-3000,3000)
  if(class(bioregionList[[1]])=="matrix"){
    
    bioregionList <- lapply(bioregionList, function(x) t(x))
    
    # to remove the chromosome that has no binding protein
    bioregionList <- bioregionList[vapply(bioregionList, function(x) length(x)>0, FUN.VALUE = logical(1))]
    suppressWarnings(bioregionList <- do.call("rbind", bioregionList))
    
    ## create a matrix to receive binning result                 
    bioregionmatrix <- matrix(nrow = dim(bioregionList)[1],ncol = box)
    
    ## this circulation is to deal with different gene             
    for (i in 1:dim(bioregionList)[1]) {
      
      ## seq is the distance between different boxes
      seq <- floor(length(bioregionList[i,])/box)
      
      ## cursor record the position of calculation
      cursor <- 1
      
      ## the third circulation is to calculate the bingding strength
      ## it has two parts
      ## the first part is to for the box(1:box-1)
      ## because the seq is not derived from exact division
      ## the second part is to compensate the loss of non-exact-division
      
      ## this the first part for 1:(box-1)
      for (j in 1:(box-1)) {
        
        read <- 0
        
        for (k in cursor:(cursor+seq-1)) {
          read <- read + bioregionList[i,k]
        }
        
        bioregionmatrix[i,j] <- read/seq
        
        cursor <- cursor+seq
      }
      
      ## this the second part to to compensate the loss of non-exact-division
      read <- 0
      for (k in cursor:length(bioregionList[i,])) {
        read <- read+bioregionList[i,k]
      }
      
      bioregionmatrix[i,box] <- read/(length(bioregionList[i,])-cursor)
    }
    
  }else{
    
    ## extend genebody by atual number
    if(!is.null(upstream) & !inherits(upstream, 'rel')){
      
      bioregionList <- bioregionList[vapply(bioregionList, function(x) length(x)>0, FUN.VALUE = logical(1))]
      
      bioregionList <- lapply(bioregionList, function(x) x[vapply(x, function(y) length(y)>min_body_length+downstream+upstream,FUN.VALUE = logical(1))])
      suppressWarnings(bioregionList <- do.call("rbind",bioregionList))
      
      bioregionmatrix <- matrix(nrow = length(bioregionList),ncol = box)
      
      upstreambox <- floor(box*(upstreamPer/(1+upstreamPer+downstreamPer)))
      bodybox <- floor(box*(1/(1+upstreamPer+downstreamPer)))
      downstreambox <- floor(box*(downstreamPer/(1+upstreamPer+downstreamPer)))
      
      # count the upstream 
      for (i in 1:length(bioregionList)) {
        
        seq <- floor(upstream/upstreambox)
        cursor <- 1
        
        for (j in 1:(upstreambox-1)) {
          
          read <- 0
          
          for (k in cursor:(cursor+seq-1)) {
            read <- read + bioregionList[[i]][k]
          }
          
          bioregionmatrix[i,j] <- read/seq
          
          cursor <- cursor+seq
        }
        
        
        read <- 0
        for (k in cursor:upstream) {
          read <- read+bioregionList[[i]][k]
        }
        
        bioregionmatrix[i,upstreambox] <- read/(upstream-cursor)
      }
      
      
      ## count genebody
      for (i in 1:length(bioregionList)) {
        
        seq <- floor((length(bioregionList[[i]])-upstream-downstream)/bodybox)
        cursor <- upstream+1
        
        for (j in (upstreambox+1):(upstreambox+bodybox-1)) {
          
          read <- 0
          
          for (k in cursor:(cursor+seq-1)) {
            read <- read + bioregionList[[i]][k]
          }
          
          bioregionmatrix[i,j] <- read/seq
          
          cursor <- cursor+seq
        }
        
        read <- 0
        for (k in cursor:(length(bioregionList[[i]])-downstream)) {
          read <- read+bioregionList[[i]][k]
        }
        
        bioregionmatrix[i,bodybox+upstreambox] <- read/(length(bioregionList[[i]])-downstream-cursor)
      }
      
      
      ## count downstream
      for (i in 1:length(bioregionList)) {
        
        seq <- floor(downstream/downstreambox)
        cursor <- length(bioregionList[[i]])-downstream+1
        
        for (j in (upstreambox+bodybox+1):(box-1)) {
          
          read <- 0
          
          for (k in cursor:(cursor+seq-1)) {
            read <- read + bioregionList[[i]][k]
          }
          
          bioregionmatrix[i,j] <- read/seq
          
          cursor <- cursor+seq
        }
        
        read <- 0
        for (k in cursor:length(bioregionList[[i]])) {
          read <- read+bioregionList[[i]][k]
        }
        
        bioregionmatrix[i,box] <- read/(length(bioregionList[[i]])-cursor)
      }
      
    }else{
      
      bioregionList <- bioregionList[vapply(bioregionList, function(x) length(x)>0, FUN.VALUE = logical(1))]
      
      bioregionList <- lapply(bioregionList, function(x) x[vapply(x, function(y) length(y)>min_body_length,FUN.VALUE = logical(1))])
      suppressWarnings(bioregionList <- do.call("rbind",bioregionList))
      
      bioregionmatrix <- matrix(nrow = length(bioregionList),ncol = box)
      
      for (i in 1:length(bioregionList)) {
        
        seq <- floor(length(bioregionList[[i]])/box)
        cursor <- 1
        
        for (j in 1:(box-1)) {
          
          read <- 0
          
          for (k in cursor:(cursor+seq-1)) {
            read <- read + bioregionList[[i]][k]
          }
          
          bioregionmatrix[i,j] <- read/seq
          
          cursor <- cursor+seq
        }
        
        read <- 0
        for (k in cursor:length(bioregionList[[i]])) {
          read <- read+bioregionList[[i]][k]
        }
        
        bioregionmatrix[i,box] <- read/(length(bioregionList[[i]])-cursor)
      }
    }
    
  }
  
  
  ## assign attribute 
  attr(bioregionmatrix, 'type') = attr(windows, 'type')
  
  return(bioregionmatrix)
}

