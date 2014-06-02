##' calcuate overlap significant of ChIP experiments based on their nearest gene annotation
##'
##' 
##' @title enrichAnnoOverlap
##' @param queryPeak query bed file
##' @param targetPeak target bed file(s) or folder containing bed files
##' @param TranscriptDb TranscriptDb
##' @param pAdjustMethod pvalue adjustment method
##' @param chainFile chain file for liftOver
##' @return data.frame
##' @export
##' @importFrom rtracklayer import.chain
##' @importFrom rtracklayer liftOver
##' @author G Yu
enrichAnnoOverlap <- function(queryPeak, targetPeak, TranscriptDb=NULL, pAdjustMethod="BH", chainFile=NULL) {
    targetFiles <- parse_targetPeak_Param(targetPeak)
    TranscriptDb <- loadTxDb(TranscriptDb)
 
    query.anno <- annotatePeak(queryPeak, TranscriptDb=TranscriptDb,
                               assignGenomicAnnotation=FALSE, annoDb=NULL, verbose=FALSE)

    target.gr <- lapply(targetFiles, loadPeak)
    if (!is.null(chainFile)) {
        chain <- import.chain(chainFile)
        target.gr <- lapply(target.gr, liftOver, chain=chain)
    }
     
    target.anno <- lapply(target.gr, annotatePeak, TranscriptDb=TranscriptDb,
                          assignGenomicAnnotation=FALSE, annoDb=NULL, verbose=FALSE)
    

    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    features <- get("features", envir=ChIPseekerEnv)
    ol <- lapply(target.anno, function(i) unique(intersect(query.anno$geneId, i$geneId)))
    oln <- unlist(lapply(ol, length))
    N <- length(features)
    ## white ball
    m <- length(unique(query.anno$geneId))
    ## black ball
    n <- N - m
    ## drawn
    k <- unlist(lapply(target.anno, function(i) length(unique(i$geneId))))
    p <- phyper(oln, m, n, k, lower.tail=FALSE)
    qSample <- sub(".+/", "", queryPeak)
    tSample <- sub(".+/", "", targetFiles)
    padj <- p.adjust(p, method=pAdjustMethod)
    res <- data.frame(qSample=qSample,
                      tSample=tSample,
                      qLen=length(unique(query.anno$geneId)),
                      tLen=unlist(lapply(target.anno, function(i) length(unique(i$geneId)))),
                      N_OL=oln,
                      pvalue=p,
                      p.adjust=padj)
    return(res)
}

##' calculate overlap significant of ChIP experiments based on the genome coordinations
##'
##' 
##' @title enrichPeakOverlap
##' @param queryPeak query bed file
##' @param targetPeak target bed file(s) or folder that containing bed files
##' @param TranscriptDb TranscriptDb
##' @param pAdjustMethod pvalue adjustment method
##' @param nShuffle shuffle numbers
##' @param chainFile chain file for liftOver
##' @return data.frame
##' @export
##' @importFrom rtracklayer import.chain
##' @importFrom rtracklayer liftOver
##' @author G Yu
enrichPeakOverlap <- function(queryPeak, targetPeak, TranscriptDb=NULL, pAdjustMethod="BH", nShuffle=1000, chainFile=NULL) {
    targetFiles <- parse_targetPeak_Param(targetPeak)
    TranscriptDb <- loadTxDb(TranscriptDb)
    query.gr <- loadPeak(queryPeak)
    target.gr <- lapply(targetFiles, loadPeak)
    if (!is.null(chainFile)) {
        chain <- import.chain(chainFile)
        target.gr <- lapply(target.gr, liftOver, chain=chain)
    }
    
    p.ol <- enrichOverlap.peak.internal(query.gr, target.gr, TranscriptDb, nShuffle)
    p <- p.ol$pvalue
    ol <- p.ol$overlap
    qSample <- sub(".+/", "", queryPeak)  ## remove path, only keep file name
    tSample <- sub(".+/", "", targetFiles) 
    padj <- p.adjust(p, method=pAdjustMethod)
    res <- data.frame(qSample=qSample,
                      tSample=tSample,
                      qLen=length(query.gr),
                      tLen=unlist(lapply(target.gr, length)),
                      N_OL=ol,
                      pvalue=p,
                      p.adjust=padj)
    
    return(res)
}



##' shuffle the position of peak
##'
##' 
##' @title shuffle
##' @param peak.gr GRanges object
##' @param TranscriptDb TranscriptDb
##' @return GRanges object
##' @export
##' @author G Yu
shuffle <- function(peak.gr, TranscriptDb) {
    chrLens <- seqlengths(TranscriptDb)[names(seqlengths(peak.gr))]
    nn <- as.vector(seqnames(peak.gr))
    ii <- order(nn)
    w <- width(peak.gr)
    nnt <- table(nn)
    jj <- order(names(nnt))
    nnt <- nnt[jj]
    chrLens <- chrLens[jj]
    ss <- unlist(sapply(1:length(nnt), function(i) sample(chrLens[i],nnt[i])))

    res <- GRanges(seqnames=nn[ii], ranges=IRanges(ss, ss+w[ii]), strand="*")
    return(res)   
}




##' @importFrom GenomeInfoDb intersect
##' @importFrom GenomeInfoDb seqlengths
enrichOverlap.peak.internal <- function(query.gr, target.gr, TranscriptDb, nShuffle=1000) {
    idx <- sample(1:length(target.gr), nShuffle, replace=TRUE)

    len <- unlist(lapply(target.gr, length))
    
    rr <- mclapply(idx, function(i) {
        tarShuffle <- shuffle(target.gr[[i]], TranscriptDb)
        length(intersect(query.gr, tarShuffle))/len[i]
    }, mc.cores=detectCores()-1
                   )
    rr <- unlist(rr) ## random ratio

    qLen <- mclapply(target.gr, function(tt) {
        length(intersect(query.gr, tt))
    }, mc.cores=detectCores()-1
                  )
    qLen <- unlist(qLen)
    ## query ratio
    qr <- qLen/len
    p <- lapply(qr, function(q) mean(rr>q))
    res <- list(pvalue=unlist(p), overlap=qLen)
    return(res)
}

