enrichOverlap <- function(peak1, peak2, by="peak", TranscriptDb=NULL, pAdjustMethod="BH", ...) {
    by <- match.arg(by, c("gene", "peak"))
    if (by == "gene") {
        res <- enrichOverlap.gene(peak1, peak2, TranscriptDb=TranscriptDb, pAdjustMethod=pAdjustMethod, ...)
    } else {
        ## not implemented yet
        res <- enrichOverlap.peak(peak1, peak2, TranscriptDb=TranscriptDb, pAdjustMethod=pAdjustMethod, ...)
    } 
    return(res)
}

enrichOverlap.peak <- function(peak1, peak2, TranscriptDb, nShuffle, pAdjustMethod="BH") {
    if (!file.exists(peak1)) {
        stop("peak1 should be peak file in bed format...")
    }
    if (! any(file.exists(unlist(peak2)))) {
        stop("peak2 should be a list of files i bed format...")
    }
    peak1.gr <- loadPeak(peak1)
    peak2 <- peak2[file.exists(unlist(peak2))] 
    peak2.gr <- lapply(peak2, loadPeak)
    p <- sapply(peak2.gr, enrichOverlap.peak.internal,
                peak1=peak1, TranscriptDb=TranscriptDb, nShuffle=nShuffle)

    qSample <- sub(".+/", "", peak1)
    tSample <- sub(".+/", "", peak2)
    ol <- unlist(lapply(peak2.gr, function(i) length(intersect(peak1.gr, i))))
    padj <- p.adjust(p, method=pAdjustMethod)
    res <- data.frame(qSample=qSample,
                      tSample=tSample,
                      qLen=length(peak1.gr),
                      tLen=unlist(lapply(peak2.gr, length)),
                      N_OL=ol,
                      pvalue=p,
                      p.adjust=padj)

    return(res)
}

enrichOverlap.peak.internal <- function(peak1, peak2, TranscriptDb, nShuffle=1000) {
    peak1.gr <- loadPeak(peak1)
    peak2.gr <- loadPeak(peak2)
    chrLens <- seqlengths(TranscriptDb)[names(seqlengths(peak2.gr))]

    shuffle.OL <- mclapply(1:nShuffle, function(i) {
        length(intersect(peak1.gr, shuffle(peak=peak2.gr, chrLens)))
    },
                             mc.cores=detectCores()
                             )
    sol <- unlist(shuffle.OL)
    ol <- length(intersect(peak1.gr, peak2.gr))

    p <- sum(sol > ol)/nShuffle
    return(p)
}

enrichOverlap.gene <- function(peak1, peak2, TranscriptDb, pAdjustMethod="BH") {
    if (!file.exists(peak1)) {
        stop("peak1 should be peak file in bed format...")
    }
    if (! any(file.exists(unlist(peak2)))) {
        stop("peak2 should be a list of files i bed format...")
    }
    peak2 <- peak2[file.exists(unlist(peak2))]
    
    peakAnno1 <- annotatePeak(peak1, TranscriptDb=TranscriptDb, assignGenomicAnnotation=FALSE, annoDb=NULL, verbose=FALSE)
    peakAnno2 <- lapply(peak2, function(i) annotatePeak(i, TranscriptDb=TranscriptDb, assignGenomicAnnotation=FALSE, annoDb=NULL, verbose=FALSE))
    
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    features <- get("features", envir=ChIPseekerEnv)
    ol <- lapply(peakAnno2, function(i) intersect(peakAnno1$geneId, i$geneId))
    oln <- unlist(lapply(ol, length))
    N <- length(features)
    ## white ball
    m <- length(unique(peakAnno1$geneId))
    ## black ball
    n <- N - m
    ## drawn
    k <- unlist(lapply(peakAnno2, function(i) length(unique(i$geneId))))
    p <- phyper(oln, m, n, k, lower.tail=FALSE)
    qSample <- sub(".+/", "", peak1)
    tSample <- sub(".+/", "", peak2)
    padj <- p.adjust(p, method=pAdjustMethod)
    res <- data.frame(qSample=qSample,
                      tSample=tSample,
                      qLen=length(peakAnno1$geneId),
                      tLen=unlist(lapply(peakAnno2, length)),
                      N_OL=oln,
                      pvalue=p,
                      p.adjust=padj)
    return(res)
}
