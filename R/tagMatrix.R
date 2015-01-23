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
##' @importFrom BiocGenerics unique
##' @importFrom GenomicRanges GRanges
##' @importFrom GenomicRanges unlist
##' @importFrom GenomicRanges strand
##' @importFrom GenomicFeatures transcriptsBy
##' @importFrom IRanges start
##' @importFrom IRanges end
##' @importFrom IRanges IRanges
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

    assign("upstream", upstream, envir=ChIPseekerEnv)
    assign("downstream", downstream, envir=ChIPseekerEnv)

    Transcripts <- getGene(TxDb, by)
    ## get start position based on strand
    tss <- ifelse(strand(Transcripts) == "+", start(Transcripts), end(Transcripts))
    promoters <- GRanges(seqnames=seqnames(Transcripts),
                         ranges=IRanges(tss-upstream, tss+downstream),
                         strand=strand(Transcripts))
    promoters <- unique(promoters)
    assign("promoters", promoters, envir=ChIPseekerEnv)
    
    return(promoters)
}



##' calculate the tag matrix
##'
##'
##' @title getTagMatrix
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight, default is NULL
##' @param windows a collection of region with equal size, eg. promoter region.
##' @return tagMatrix
##' @export
##' @importFrom IRanges subsetByOverlaps
##' @importFrom IRanges elementLengths
##' @importFrom IRanges width
##' @importFrom IRanges coverage
##' @importFrom IRanges IRanges
##' @importFrom IRanges Views
##' @importFrom IRanges viewApply
##' @importFrom IRanges as.vector
##' @importFrom S4Vectors as.factor
##' @importFrom S4Vectors mcols
##' @importFrom GenomicRanges GRanges
##' @importFrom GenomeInfoDb seqnames
##' @importFrom BiocGenerics intersect
##' @importFrom BiocGenerics unique
getTagMatrix <- function(peak, weightCol=NULL, windows) {
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
    cov.len <- elementLengths(peak.cov)
    cov.width <- GRanges(seqnames=names(cov.len),
                         IRanges(start=rep(1, length(cov.len)),
                                 end=cov.len))
    windows <- subsetByOverlaps(windows, cov.width,
                                type="within", ignore.strand=TRUE)

    chr.idx <- intersect(names(peak.cov),
                         unique(as.character(seqnames(windows))))
    
    peakView <- Views(peak.cov[chr.idx], as(windows, "RangesList")[chr.idx])
    tagMatrixList <- lapply(peakView, function(x) t(viewApply(x, as.vector)))
    tagMatrix <- do.call("rbind", tagMatrixList)

    ## get the index of windows, that are reorganized by as(windows, "RangesList")
    idx.list <- split(1:length(windows),  as.factor(seqnames(windows)))
    idx <- do.call("c", idx.list)
    
    rownames(tagMatrix) <- idx
    tagMatrix <- tagMatrix[order(idx),]
    
    ## minus strand
    minus.idx <- which(as.character(strand(windows)) == "-")
    tagMatrix[minus.idx,] <- tagMatrix[minus.idx, ncol(tagMatrix):1]
    tagMatrix <- tagMatrix[rowSums(tagMatrix)!=0,]
    ## assign("tagMatrix", tagMatrix, envir=ChIPseekerEnv)
    return(tagMatrix)
}


