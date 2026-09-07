##' Calculate overlap significance of ChIP experiments based on nearest gene annotation
##'
##' This function tests whether two ChIP-seq experiments share significantly
##' more genes than expected by chance, based on the genes nearest to their peaks.
##'
##' @description
##' The function compares the overlap of genes associated with query and target
##' ChIP-seq peaks. It uses a hypergeometric test to determine if the number of
##' overlapping genes is significantly higher than expected by chance, suggesting
##' biological similarity between the experiments.
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Annotates both query and target peaks to find their nearest genes
##'   \item Optionally filters peaks by distance to TSS using
##'         \code{distanceToTSS_cutoff}
##'   \item Counts the number of overlapping genes between query and target
##'   \item Performs a hypergeometric test to calculate p-values
##'   \item Adjusts p-values for multiple testing
##' }
##'
##' The hypergeometric test parameters:
##' \itemize{
##'   \item White balls (m): number of unique genes in query peaks
##'   \item Black balls (n): total genes minus query genes
##'   \item Drawn (k): number of unique genes in target peaks
##'   \item Overlap (q): number of overlapping genes
##' }
##'
##' @param queryPeak character, path to query peak file (BED format), or GRanges
##'   object containing query peaks
##' @param targetPeak character vector of target peak file paths, a folder
##'   containing bed files, or a list of GRanges objects. Multiple target peaks
##'   can be provided for comparison
##' @param TxDb TxDb or EnsDb annotation object (already loaded in R). If NULL,
##'   uses \code{TxDb.Hsapiens.UCSC.hg19.knownGene} as default. Must be an
##'   already-loaded R object, not a file path. To load from a file, use
##'   \code{GenomicFeatures::makeTxDbFromGFF()}, \code{AnnotationDbi::loadDb()},
##'   or similar functions first
##' @param pAdjustMethod character, method for p-value adjustment. Default is
##'   "BH" (Benjamini-Hochberg). See \code{\link[stats]{p.adjust}} for options
##' @param chainFile character, path to chain file for liftOver conversion if
##'   target peaks are in a different genome assembly. Default is NULL
##' @param distanceToTSS_cutoff numeric, distance cutoff in base pairs. Peaks
##'   with absolute distance to TSS greater than this value will be excluded from
##'   the analysis. If NULL, all peaks are included. Default is NULL
##' @return A data.frame with one row per target peak comparison, containing:
##'   \itemize{
##'     \item \code{qSample}: name of the query sample
##'     \item \code{tSample}: name of the target sample(s)
##'     \item \code{qLen}: number of unique genes in query peaks
##'     \item \code{tLen}: number of unique genes in target peaks
##'     \item \code{N_OL}: number of overlapping genes
##'     \item \code{pvalue}: p-value from hypergeometric test
##'     \item \code{p.adjust}: adjusted p-value
##'   }
##' @importFrom stats p.adjust
##' @importFrom stats phyper
##' @export
##' @importFrom rtracklayer import.chain
##' @importFrom rtracklayer liftOver
##' @importFrom yulab.utils get_cache_element
##' @importFrom yulab.utils update_cache_item
##' @author G Yu
enrichAnnoOverlap <- function(queryPeak, targetPeak, TxDb=NULL, pAdjustMethod="BH", chainFile=NULL, distanceToTSS_cutoff=NULL) {

    TxDb <- loadTxDb(TxDb)

    query.anno <- annotatePeak(queryPeak, TxDb=TxDb,
                               assignGenomicAnnotation=FALSE, annoDb=NULL, verbose=FALSE)


    if (is(targetPeak[1], "GRanges") || is(targetPeak[[1]], "GRanges")) {
        target.gr <- targetPeak
        targetFiles <- NULL
    } else {
        targetFiles <- parse_targetPeak_Param(targetPeak)
        target.gr <- lapply(targetFiles, loadPeak)
    }

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

    # ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    features <- get_cache_element(item = ChIPseekerCache, elements = "Transcripts")

    if(is.null(features)){
        features <- transcriptsBy(TxDb)
        features <- unlist(features)
        update_cache_item(item = ChIPseekerCache, list("Transcripts" = features))
    }

    # if ( exists("Transcripts", envir=ChIPseekerEnv, inherits=FALSE) ) {
    #     features <- get("Transcripts", envir=ChIPseekerEnv)
    # } else {
    #     features <- transcriptsBy(TxDb)
    #     features <- unlist(features)
    #     assign("Transcripts", features, envir=ChIPseekerEnv)
    # }

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


    if (is(queryPeak, "GRanges")) {
        qSample <- "queryPeak"
    } else {
        qSample <- basename(queryPeak)
    }

    if (is.null(targetFiles)) {
        tSample <- names(target.gr)
        if(is.null(tSample)) {
            tSample <- paste0("targetPeak", seq_along(target.gr))
        }
    } else {
        tSample <- basename(targetFiles)
    }

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

##' Calculate overlap significance of ChIP experiments based on genomic coordinates
##'
##' This function tests whether two ChIP-seq experiments have significantly more
##' overlapping peaks than expected by chance, using permutation testing with
##' shuffled peak positions.
##'
##' @description
##' The function compares the genomic overlap between query and target ChIP-seq
##' peaks. It uses a permutation test approach: randomly shuffling target peaks
##' across the genome and comparing the observed overlap to the distribution of
##' overlaps from shuffled data. This provides a p-value indicating whether the
##' observed overlap is significantly higher than expected by chance.
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Loads query and target peaks (supports BED files or GRanges objects)
##'   \item Optionally converts target peaks to query genome using liftOver
##'         (if \code{chainFile} is provided)
##'   \item If \code{pool=TRUE}, pools all target peaks and tests overlap with
##'         the pooled set
##'   \item If \code{pool=FALSE}, tests overlap with each target peak set separately
##'   \item For permutation testing:
##'     \itemize{
##'       \item Calculates the observed overlap ratio (overlapping peaks / total
##'             target peaks)
##'       \item Randomly shuffles target peaks \code{nShuffle} times across the
##'             genome (preserving chromosome and peak width)
##'       \item Calculates overlap ratio for each shuffled set
##'       \item Computes p-value as (number of shuffled ratios >= observed ratio + 1) /
##'             (nShuffle + 1)
##'     }
##'   \item Adjusts p-values for multiple testing
##' }
##'
##' @param queryPeak character, path to query peak file (BED format), or GRanges
##'   object containing query peaks
##' @param targetPeak character vector of target peak file paths, a folder
##'   containing bed files, or a list of GRanges objects. Multiple target peaks
##'   can be provided
##' @param TxDb TxDb or EnsDb annotation object (already loaded in R). Required
##'   for shuffling peaks (to get chromosome lengths). If NULL, uses
##'   \code{TxDb.Hsapiens.UCSC.hg19.knownGene} as default. Must be an
##'   already-loaded R object, not a file path
##' @param pAdjustMethod character, method for p-value adjustment. Default is
##'   "BH" (Benjamini-Hochberg). See \code{\link[stats]{p.adjust}} for options
##' @param nShuffle integer, number of permutations for the statistical test.
##'   More shuffles provide more accurate p-values but take longer. Default is 1000.
##'   Set to 0 to skip permutation testing (p-value will be NA)
##' @param chainFile character, path to chain file for liftOver conversion if
##'   target peaks are in a different genome assembly. Default is NULL
##' @param pool logical, whether to pool all target peaks together for a single
##'   test. If TRUE, tests overlap with the combined set of all target peaks.
##'   If FALSE, tests overlap with each target peak set separately. Default is TRUE
##' @param mc.cores integer, number of CPU cores to use for parallel processing.
##'   Default is \code{detectCores()-1}. See \code{\link[parallel]{mclapply}}
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @return A data.frame with one row per target peak comparison, containing:
##'   \itemize{
##'     \item \code{qSample}: name of the query sample
##'     \item \code{tSample}: name of the target sample(s)
##'     \item \code{qLen}: number of peaks in query
##'     \item \code{tLen}: number of peaks in target (or pooled target if pool=TRUE)
##'     \item \code{N_OL}: number of overlapping peaks
##'     \item \code{pvalue}: p-value from permutation test (NA if nShuffle=0)
##'     \item \code{p.adjust}: adjusted p-value
##'   }
##' @export
##' @importFrom rtracklayer import.chain
##' @importFrom rtracklayer liftOver
##' @author G Yu
enrichPeakOverlap <- function(queryPeak, targetPeak, TxDb=NULL, pAdjustMethod="BH", nShuffle=1000,
                              chainFile=NULL, pool=TRUE, mc.cores=detectCores()-1, verbose=TRUE) {
    TxDb <- loadTxDb(TxDb)
    query.gr <- loadPeak(queryPeak)
    if (is(targetPeak[1], "GRanges") || is(targetPeak[[1]], "GRanges")) {
        target.gr <- targetPeak
        targetFiles <- NULL
    } else {
        targetFiles <- parse_targetPeak_Param(targetPeak)
        target.gr <- lapply(targetFiles, loadPeak)
    }

    if (!is.null(chainFile)) {
        chain <- import.chain(chainFile)
        target.gr <- lapply(target.gr, liftOver, chain=chain)
    }

    if (pool) {
        p.ol <- enrichOverlap.peak.internal(query.gr, target.gr, TxDb, nShuffle,
                                            mc.cores=mc.cores,verbose=verbose)
    } else {
        res_list <- lapply(1:length(target.gr), function(i) {
            enrichPeakOverlap(queryPeak = queryPeak,
                              targetPeak = target.gr[i],
                              TxDb = TxDb,
                              pAdjustMethod = pAdjustMethod,
                              nShuffle = nShuffle,
                              chainFile = chainFile,
                              mc.cores = mc.cores,
                              verbose = verbose)
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


    if (is(queryPeak, "GRanges")) {
        qSample <- "queryPeak"
    } else {
        ## remove path, only keep file name
        qSample <- basename(queryPeak)
    }

    if (is.null(targetFiles)) {
        tSample <- names(target.gr)
        if(is.null(tSample)) {
            tSample <- paste0("targetPeak", seq_along(target.gr))
        }
    } else {
        tSample <- basename(targetFiles)
    }

    res <- data.frame(qSample=qSample,
                      tSample=tSample,
                      qLen=length(query.gr),
                      tLen=unlist(lapply(target.gr, length)),
                      N_OL=ol,
                      pvalue=p,
                      p.adjust=padj)

    return(res)
}



##' Shuffle peak positions across the genome
##'
##' This function randomly shuffles peak positions across chromosomes while
##' preserving the original peak widths and chromosome distribution.
##'
##' @description
##' The function generates a null distribution for permutation testing by randomly
##' repositioning peaks within their respective chromosomes. It maintains:
##' \itemize{
##'   \item The number of peaks per chromosome
##'   \item The width of each peak
##'   \item The chromosome assignment
##' }
##'
##' @details
##' The shuffling process:
##' \enumerate{
##'   \item Groups peaks by chromosome
##'   \item For each chromosome, randomly samples new start positions from the
##'         valid range (1 to chromosome length - peak width)
##'   \item Creates new GRanges with shuffled positions but original widths
##'   \item Sets strand to "*" (unstranded)
##' }
##'
##' This is used in permutation testing to generate null distributions for
##' statistical significance testing.
##'
##' @param peak.gr GRanges object containing peaks to shuffle
##' @param TxDb TxDb or EnsDb annotation object providing chromosome lengths
##'   for valid position ranges
##' @return GRanges object with shuffled peak positions. Peaks are randomly
##'   repositioned within their chromosomes, maintaining original widths and
##'   chromosome distribution
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

    res <- GRanges(seqnames=nn[ii], ranges=IRanges(ss, width=w[ii]), strand="*")
    return(res)
}




##' @import GenomeInfoDb
##' @importFrom utils txtProgressBar
##' @importFrom utils setTxtProgressBar
##' @importFrom parallel mclapply
##' @importFrom parallel detectCores
##' @noRd
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

