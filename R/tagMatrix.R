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
    
    return(promoters)
}

##' prepare a region center on start site of selected feature
##'
##' 
##' @title getBioRegion
##' @param TxDb TxDb
##' @param upstream upstream from start site
##' @param downstream downstream from start site
##' @param by one of 'gene', 'transcript', 'exon', 'intron'
##' @return GRanges object
##' @import BiocGenerics IRanges GenomicRanges
##' @export
##' @author Guangchuang Yu
##  https://github.com/GuangchuangYu/ChIPseeker/issues/16
getBioRegion <- function(TxDb=NULL,
                         upstream=1000,
                         downstream=1000,
                         by="gene") {
    
    by <- match.arg(by, c("gene", "transcript", "exon", "intron"))

    if (by %in% c("gene", "transcript")) {
        return(getPromoters(TxDb, upstream, downstream, by))
    }
    
    TxDb <- loadTxDb(TxDb)
    .ChIPseekerEnv(TxDb)
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

    if (by == "exon") {
        exonList <- get_exonList(ChIPseekerEnv)
        regions <-  unlist(exonList)
    }

    if (by == "intron") {
        intronList <- get_intronList(ChIPseekerEnv)
        regions <- unlist(intronList)
    }

    start_site <- ifelse(strand(regions) == "+", start(regions), end(regions))

    bioRegion <- GRanges(seqnames=seqnames(regions),
                         ranges=IRanges(start_site-upstream, start_site+downstream),
                         strand=strand(regions))
    bioRegion <- unique(bioRegion)

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


