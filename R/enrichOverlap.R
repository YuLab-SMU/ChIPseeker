##' calcuate overlap significant of ChIP experiments based on their nearest gene annotation
##'
##' 
##' @title enrichAnnoOverlap
##' @param queryPeak query bed file
##' @param targetPeak target bed file(s) or folder containing bed files
##' @param TxDb TxDb
##' @param pAdjustMethod pvalue adjustment method
##' @param chainFile chain file for liftOver
##' @param distanceToTSS_cutoff restrict nearest gene annotation by distance cutoff
##' @return data.frame
##' @export
##' @importFrom rtracklayer import.chain
##' @importFrom rtracklayer liftOver
##' @author G Yu
enrichAnnoOverlap <- function(queryPeak, targetPeak, TxDb=NULL, pAdjustMethod="BH", chainFile=NULL, distanceToTSS_cutoff=NULL) {
    targetFiles <- parse_targetPeak_Param(targetPeak)
    TxDb <- loadTxDb(TxDb)
 
    query.anno <- annotatePeak(queryPeak, TxDb=TxDb,
                               assignGenomicAnnotation=FALSE, annoDb=NULL, verbose=FALSE)

    target.gr <- lapply(targetFiles, loadPeak)
    if (!is.null(chainFile)) {
        chain <- import.chain(chainFile)
        target.gr <- lapply(target.gr, liftOver, chain=chain)
    }
     
    target.anno <- lapply(target.gr, annotatePeak, TxDb=TxDb,
                          assignGenomicAnnotation=FALSE, annoDb=NULL, verbose=FALSE)
    

    if (!is.null(distanceToTSS_cutoff)) {
        query.anno <- dropAnno(query.anno, distanceToTSS_cutoff)
        target.anno <- lapply(target.anno, dropAnno, distanceToTSS_cutoff = distanceToTSS_cutoff)
    }
    
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    if ( exists("Transcripts", envir=ChIPseekerEnv, inherits=FALSE) ) {
        features <- get("Transcripts", envir=ChIPseekerEnv)
    } else {
        features <- transcriptsBy(TxDb)
        features <- unlist(features)
        assign("Transcripts", features, envir=ChIPseekerEnv)
    }
    
    ol <- lapply(target.anno, function(i) unique(intersect(as.GRanges(query.anno)$geneId, as.GRanges(i)$geneId)))
    oln <- unlist(lapply(ol, length))
    N <- length(features)
    ## white ball
    m <- length(unique(as.GRanges(query.anno)$geneId))
    ## black ball
    n <- N - m
    ## drawn
    k <- unlist(lapply(target.anno, function(i) length(unique(as.GRanges(i)$geneId))))
    p <- phyper(oln, m, n, k, lower.tail=FALSE)
    qSample <- sub(".+/", "", queryPeak)
    tSample <- sub(".+/", "", targetFiles)
    padj <- p.adjust(p, method=pAdjustMethod)
    res <- data.frame(qSample=qSample,
                      tSample=tSample,
                      qLen=length(unique(as.GRanges(query.anno)$geneId)),
                      tLen=unlist(lapply(target.anno, function(i) length(unique(as.GRanges(i)$geneId)))),
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
##' @param TxDb TxDb
##' @param pAdjustMethod pvalue adjustment method
##' @param nShuffle shuffle numbers
##' @param chainFile chain file for liftOver
##' @param pool logical, whether pool target peaks
##' @param ... additional parameter
##' @return data.frame
##' @export
##' @importFrom rtracklayer import.chain
##' @importFrom rtracklayer liftOver
##' @author G Yu
enrichPeakOverlap <- function(queryPeak, targetPeak, TxDb=NULL, pAdjustMethod="BH", nShuffle=1000, chainFile=NULL, pool=TRUE, ...) {
    targetFiles <- parse_targetPeak_Param(targetPeak)
    TxDb <- loadTxDb(TxDb)
    query.gr <- loadPeak(queryPeak)
    target.gr <- lapply(targetFiles, loadPeak)
     if (!is.null(chainFile)) {
        chain <- import.chain(chainFile)
        target.gr <- lapply(target.gr, liftOver, chain=chain)
    }

    if (pool) {
        p.ol <- enrichOverlap.peak.internal(query.gr, target.gr, TxDb, nShuffle, ...)
    } else {
        res_list <- lapply(1:length(targetFiles), function(i) {
            enrichPeakOverlap(queryPeak = queryPeak,
                              targetPeak = targetFiles[i],
                              TxDb = TxDb,
                              pAdjustMethod = pAdjustMethod,
                              nShuffle = nShuffle,
                              chainFile = chainFile,
                              ...)
        })
        res <- do.call("rbind", res_list)
        return(res)
    }
    
    if (is.null(p.ol$pvalue)) {
        p <- padj <- NA
    } else {
        p <- p.ol$pvalue
        padj <- p.adjust(p, method=pAdjustMethod)
    }
        
    ol <- p.ol$overlap
    qSample <- sub(".+/", "", queryPeak)  ## remove path, only keep file name
    tSample <- sub(".+/", "", targetFiles) 
    
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
##' @param TxDb TxDb
##' @return GRanges object
##' @export
##' @author G Yu
shuffle <- function(peak.gr, TxDb) {
    chrLens <- seqlengths(TxDb)[names(seqlengths(peak.gr))]
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
##' @importFrom parallel mclapply
##' @importFrom parallel detectCores
enrichOverlap.peak.internal <- function(query.gr, target.gr, TxDb, nShuffle=1000, mc.cores=detectCores()-1, verbose=TRUE) {    
    if (verbose) {
        cat(">> permutation test of peak overlap...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }

    idx <- sample(1:length(target.gr), nShuffle, replace=TRUE)
    len <- unlist(lapply(target.gr, length))

    if(Sys.info()[1] == "Windows") {
        qLen <- lapply(target.gr, function(tt) {
            length(intersect(query.gr, tt))
        })
    } else {
        qLen <- mclapply(target.gr, function(tt) {
            length(intersect(query.gr, tt))
        }, mc.cores=mc.cores
                         )
    }
    qLen <- unlist(qLen)
    ## query ratio
    qr <- qLen/len
    
    if (nShuffle < 1) {
        res <- list(pvalue=NULL, overlap=qLen)
        return(res)
    }

    if (verbose) {
        pb <- txtProgressBar(min=0, max=nShuffle, style=3)
    }
    if(Sys.info()[1] == "Windows") {
        rr <- lapply(seq_along(idx), function(j) {
            if (verbose) {
                setTxtProgressBar(pb, j)
            }
            i <- idx[j]
            tarShuffle <- shuffle(target.gr[[i]], TxDb)
            length(intersect(query.gr, tarShuffle))/len[i]
        })
    } else {
        rr <- mclapply(seq_along(idx), function(j) {
            if (verbose) {
                setTxtProgressBar(pb, j)
            }
            i <- idx[j]
            tarShuffle <- shuffle(target.gr[[i]], TxDb)
            length(intersect(query.gr, tarShuffle))/len[i]
        }, mc.cores=mc.cores
                       )
    }

    if (verbose) {
        close(pb)
    }
    
    rr <- unlist(rr) ## random ratio

    ## p <- lapply(qr, function(q) mean(rr>q))
    p <- lapply(qr, function(q) (sum(rr>q)+1)/(length(rr)+1))
    res <- list(pvalue=unlist(p), overlap=qLen)
    return(res)
}

