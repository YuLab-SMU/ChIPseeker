##' Initialize and manage ChIPseeker cache environment
##'
##' This function manages the caching system for TxDb objects and related
##' genomic features to improve performance by avoiding redundant computations.
##'
##' @description
##' The function initializes and maintains a cache for TxDb objects and their
##' associated genomic features. It checks if the cached TxDb matches the
##' provided one, and updates the cache if they differ or if forced to update.
##'
##' @details
##' The function performs the following:
##' \enumerate{
##'   \item Retrieves or creates a cache item
##'   \item If no TxDb is cached, stores the provided TxDb
##'   \item If a TxDb exists in cache, compares metadata to check if it matches
##'   \item Updates cache if TxDb differs or if \code{force=TRUE}
##'   \item Prints the genome version being used
##' }
##'
##' This caching mechanism allows the package to reuse computed features
##' (exons, introns, transcripts, etc.) across multiple function calls when
##' using the same TxDb object.
##'
##' @param TxDb TxDb or EnsDb annotation object to cache
##' @param item character, name of the cache item. Default is "ChIPseekerEnv"
##' @param force logical, whether to force update of the TxDb in cache,
##'   clearing all cached features. Default is FALSE
##' @return Invisibly returns NULL. Prints messages about genome version usage
#' @importFrom yulab.utils get_cache_item
#' @importFrom yulab.utils update_cache_item
#' @importFrom yulab.utils rm_cache_item
#' @importFrom yulab.utils initial_cache_item
#' @importFrom S4Vectors metadata
.ChIPseekerEnv <- function(TxDb, item = "ChIPseekerEnv", force = FALSE) {

    # get cache item
    # it will create a list if there is no a cache item
    cache_item <- get_cache_item(item)

    # if there is no TXDB cached, write in cache
    if (is.null(cache_item$TXDB)) {
        update_cache_item(item = item, list(TXDB = TxDb))
        cat(">> Using Genome:", get_env_genome(),"...\n")
        return(invisible(NULL))
    }

    # force to update item
    if(force){
        cat(">> Force to update txdb in cache...\n")
        rm_cache_item(item)
        initial_cache_item(item)
        update_cache_item(item, list(TXDB = TxDb))
        cat(">> Using Genome:", get_env_genome(),"...\n")
    }

    # if exist TXDB
    TXDB <- cache_item$TXDB
    m1 <- tryCatch(unlist(metadata(TXDB)), error = function(e) NULL)
    m2 <- tryCatch(unlist(metadata(TxDb)),  error = function(e) NULL)
    if (!is.null(m1)) m1 <- m1[!is.na(m1)]
    if (!is.null(m2)) m2 <- m2[!is.na(m2)]

    txdb_flag <- is.character(all.equal(TXDB, TxDb))

    if (is.null(m1) || is.null(m2) || length(m1) != length(m2) || any(m1 != m2) || txdb_flag) {
        cat(">> Update txdb in cache...\n")
        rm_cache_item(item)
        initial_cache_item(item)
        update_cache_item(item, list(TXDB = TxDb))
    }

    cat(">> Using Genome:", get_env_genome(),"...\n")

    invisible(NULL)

    # pos <- 1
    # envir <- as.environment(pos)
    # if (!exists("ChIPseekerEnv", envir=.GlobalEnv)) {
    #     assign("ChIPseekerEnv", new.env(), envir = envir)
    # }

    # ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    # if (!exists("TXDB", envir=ChIPseekerEnv, inherits=FALSE)) {
    #     ## first run
    #     assign("TXDB", TxDb, envir=ChIPseekerEnv)
    # } else {
    #     TXDB <- get("TXDB", envir=ChIPseekerEnv)
    #     m1 <- tryCatch(unlist(metadata(TXDB)), error=function(e) NULL)

    #     m2 <- unlist(metadata(TxDb))

    #     if (!is.null(m1)) {
    #         m1 <- m1[!is.na(m1)]
    #     }
    #     m2 <- m2[!is.na(m2)]

    #     if ( is.null(m1) || length(m1) != length(m2) || any(m1 != m2) ) {
    #         rm(ChIPseekerEnv)
    #         assign("ChIPseekerEnv", new.env(), envir = envir)
    #         ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    #         assign("TXDB", TxDb, envir=ChIPseekerEnv)
    #     }
    # }
}


##' Get exon list from cache or compute if needed
##'
##' This function retrieves the exon list from cache or computes it from the
##' cached TxDb object if not available.
##'
##' @description
##' The function checks the cache for a pre-computed exon list. If not found,
##' it extracts exons from the cached TxDb using \code{exonsBy()} and stores
##' the result in cache for future use.
##'
##' @param item character, name of the cache item. Default is "ChIPseekerEnv"
##' @return A GRangesList object containing exons grouped by transcript. Each
##'   element is a GRanges object with exons for one transcript
##' @importFrom GenomicFeatures exonsBy
##' @importFrom yulab.utils get_cache_element
##' @importFrom yulab.utils update_cache_item
##' @noRd
get_exonList <- function(item = "ChIPseekerEnv") {
    # TxDb <- get("TXDB", envir=ChIPseekerEnv)
    TxDb <- get_cache_element(item = item, elements = "TXDB")

    exonList <- get_cache_element(item = item, elements = "exonList")

    if(is.null(exonList)){
        exonList <- exonsBy(TxDb)
        update_cache_item(item = item, list("exonList" = exonList))
    }

    # if ( exists("exonList", envir=ChIPseekerEnv, inherits=FALSE) ) {
    #     exonList <- get("exonList", envir=ChIPseekerEnv)
    # } else {
    #     exonList <- exonsBy(TxDb)
    #     assign("exonList", exonList, envir=ChIPseekerEnv)
    # }
    return(exonList)
}

##' Get intron list from cache or compute if needed
##'
##' This function retrieves the intron list from cache or computes it from the
##' cached TxDb object if not available.
##'
##' @description
##' The function checks the cache for a pre-computed intron list. If not found,
##' it extracts introns from the cached TxDb using \code{intronsByTranscript()}
##' and stores the result in cache for future use.
##'
##' @param item character, name of the cache item. Default is "ChIPseekerEnv"
##' @return A GRangesList object containing introns grouped by transcript. Each
##'   element is a GRanges object with introns for one transcript
##' @importFrom GenomicFeatures intronsByTranscript
##' @importFrom yulab.utils get_cache_element
##' @importFrom yulab.utils update_cache_item
##' @noRd
get_intronList <- function(item = "ChIPseekerEnv") {

    # TxDb <- get("TXDB", envir=ChIPseekerEnv)
    TxDb <- get_cache_element(item = item, elements = "TXDB")

    intronList <- get_cache_element(item = item, elements = "intronList")

    if(is.null(intronList)){
        intronList <- intronsByTranscript(TxDb)
        update_cache_item(item = item, list("intronList" = intronList))
    }

    # if ( exists("intronList", envir=ChIPseekerEnv, inherits=FALSE) ) {
    #     intronList <- get("intronList", envir=ChIPseekerEnv)
    # } else {
    #     intronList <- intronsByTranscript(TxDb)
    #     assign("intronList", intronList, envir=ChIPseekerEnv)
    # }
    return(intronList)
}


##' Get color palette for plotting
##'
##' Returns a vector of colors from predefined palettes for use in plots.
##'
##' @param n integer, number of colors needed
##' @return Character vector of color codes (hex format) of length \code{n}
##' @noRd
getCols <- function(n) {
    col <- c("#8dd3c7", "#ffffb3", "#bebada",
             "#fb8072", "#80b1d3", "#fdb462",
             "#b3de69", "#fccde5", "#d9d9d9",
             "#bc80bd", "#ccebc5", "#ffed6f")

    col2 <- c("#1f78b4", "#ffff33", "#c2a5cf",
             "#ff7f00", "#810f7c", "#a6cee3",
             "#006d2c", "#4d4d4d", "#8c510a",
             "#d73027", "#78c679", "#7f0000",
             "#41b6c4", "#e7298a", "#54278f")

    col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
              "#33a02c", "#fb9a99", "#e31a1c",
              "#fdbf6f", "#ff7f00", "#cab2d6",
              "#6a3d9a", "#ffff99", "#b15928")

    ## colorRampPalette(brewer.pal(12, "Set3"))(n)
    col3[1:n]
}

##' Get color palette names for RColorBrewer
##'
##' Returns a vector of diverging color palette names from RColorBrewer.
##'
##' @param n integer, number of palette names needed
##' @return Character vector of palette names of length \code{n}
##' @noRd
getPalette <- function(n){

  palette <- c("RdBu", "RdYlGn", "Spectral",
               "RdYlBu", "PiYG", "PRGn",
               "PuOr", "BrBG", "RdGy")

  palette[1:n]

}

getSgn <- function(data, idx){
    d <- data[idx, ]
    ss <- colSums(d)
    ss <- ss / sum(ss)
    return(ss)
}
parseBootCiPerc <- function(bootCiPerc){
    bootCiPerc <- bootCiPerc$percent
    tmp <- length(bootCiPerc)
    ciLo <- bootCiPerc[tmp - 1]
    ciUp <- bootCiPerc[tmp]
    return(c(ciLo, ciUp))
}

##' Estimate confidence intervals for tag matrix using bootstrapping
##'
##' This function calculates confidence intervals for each position in a tag
##' matrix using bootstrap resampling.
##'
##' @description
##' The function performs bootstrap resampling on the tag matrix to estimate
##' confidence intervals for the signal at each genomic position. This is useful
##' for visualizing uncertainty in tag density profiles.
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Performs bootstrap resampling (with replacement) of the tag matrix rows
##'   \item Calculates the normalized signal (proportions) for each resample
##'   \item Computes percentile-based confidence intervals for each position
##'   \item Returns lower and upper bounds for each position
##' }
##'
##' @param tagMatrix matrix, rows represent peaks/regions and columns represent
##'   genomic positions. Values are tag counts or signal intensities
##' @param conf numeric, confidence level (e.g., 0.95 for 95% CI). Default is 0.95
##' @param resample integer, number of bootstrap resamples. More resamples provide
##'   more accurate CIs but take longer. Default is 500
##' @param ncpus integer, number of CPU cores to use for parallel processing.
##'   Default is \code{detectCores()-1}
##' @return A matrix with 2 rows ("Lower" and "Upper") and columns corresponding
##'   to genomic positions, containing the confidence interval bounds
##' @importFrom boot boot
##' @importFrom boot boot.ci
##' @importFrom parallel detectCores
##' @noRd
getTagCiMatrix <- function(tagMatrix, conf = 0.95, resample=500, ncpus=detectCores()-1){
    RESAMPLE_TIME <- resample
    trackLen <- ncol(tagMatrix)
    if (Sys.info()[1] == "Windows") {
        tagMxBoot <- boot(data = tagMatrix, statistic = getSgn, R = RESAMPLE_TIME)
    } else {
        tagMxBoot <- boot(data = tagMatrix, statistic = getSgn, R = RESAMPLE_TIME,
                          parallel = "multicore", ncpus = ncpus)
    }
    cat(">> Running bootstrapping for tag matrix...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
    tagMxBootCi <- sapply(seq_len(trackLen), function(i) {
                        bootCiToken <- boot.ci(tagMxBoot, type = "perc", index = i)
                        ## parse boot.ci results
                        return(parseBootCiPerc(bootCiToken))
                        }
                    )
    row.names(tagMxBootCi) <- c("Lower", "Upper")
    return(tagMxBootCi)
}

##' Get normalized tag counts from tag matrix
##'
##' This function calculates normalized tag counts (proportions) from a tag matrix,
##' optionally with confidence intervals.
##'
##' @description
##' The function sums tag counts across all peaks/regions for each genomic
##' position and normalizes by the total to get proportions. Optionally calculates
##' confidence intervals using bootstrapping.
##'
##' @param tagMatrix matrix, rows represent peaks/regions and columns represent
##'   genomic positions
##' @param xlim numeric vector of length 2, genomic range (start, end) to extract
##' @param conf numeric, confidence level for bootstrapping. If missing or NA,
##'   confidence intervals are not calculated
##' @param ... additional arguments passed to \code{getTagCiMatrix()}
##' @return A data.frame with columns:
##'   \itemize{
##'     \item \code{pos}: genomic positions
##'     \item \code{value}: normalized tag counts (proportions)
##'     \item \code{Lower}, \code{Upper}: confidence interval bounds (if \code{conf}
##'           is provided)
##'   }
##' @noRd
getTagCount <- function(tagMatrix, xlim, conf, ...) {
    ss <- colSums(tagMatrix)
    ss <- ss/sum(ss)
    ## plot(1:length(ss), ss, type="l", xlab=xlab, ylab=ylab)
    pos <- value <- NULL
    dd <- data.frame(pos=c(xlim[1]:xlim[2]), value=ss)
    if (!(missingArg(conf) || is.na(conf))){
        tagCiMx <- getTagCiMatrix(tagMatrix, conf = conf, ...)
        dd$Lower <- tagCiMx["Lower", ]
        dd$Upper <- tagCiMx["Upper", ]
    }
    return(dd)
}


##' Convert transcript IDs to gene IDs
##'
##' This function converts transcript IDs to gene IDs, with options to return
##' either transcript-gene pairs or gene IDs only.
##'
##' @description
##' The function maps transcript IDs to their corresponding gene IDs from the
##' cached TxDb object. It can return either formatted transcript/gene pairs
##' (e.g., "uc001aed.3/126789") or just the gene IDs.
##'
##' @param txid character vector of transcript IDs to convert
##' @param geneIdOnly logical, if TRUE returns only gene IDs. If FALSE, returns
##'   transcript/gene pairs in format "transcript_name/gene_id". Default is FALSE
##' @return Character vector of converted IDs. If \code{geneIdOnly=FALSE}, returns
##'   transcript/gene pairs. If \code{geneIdOnly=TRUE}, returns only gene IDs
##' @noRd
TXID2EG <- function(txid, geneIdOnly=FALSE) {
    txid <- as.character(txid)
    if (geneIdOnly == TRUE) {
        res <- TXID2EGID(txid)
    } else {
        res <- TXID2TXEG(txid)
    }
    return(res)
}

##' Convert transcript IDs to transcript/gene pairs
##'
##' This function converts transcript IDs to formatted strings containing both
##' transcript names and gene IDs (e.g., "uc001aed.3/126789").
##'
##' @param txid character vector of transcript IDs
##' @return Character vector of transcript/gene pairs in format "transcript_name/gene_id"
##' @importFrom GenomicFeatures transcripts
##' @importFrom yulab.utils get_cache_element
##' @importFrom yulab.utils update_cache_item
##' @noRd
TXID2TXEG <- function(txid) {
    # ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

    txid2geneid <- get_cache_element(item = ChIPseekerCache, elements = "txid2geneid")

    if(is.null(txid2geneid)){
        txdb <- get_cache_element(item = ChIPseekerCache, elements = "TXDB")
        txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
        idx <- which(sapply(txidinfo$gene_id, length) == 0)
        txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
        txid2geneid <- paste(mcols(txidinfo)[["tx_name"]],
                             mcols(txidinfo)[["gene_id"]],
                             sep="/")
        txid2geneid <- sub("/NA", "", txid2geneid)

        names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
        update_cache_item(item = ChIPseekerCache, list("txid2geneid" = txid2geneid))
    }

    # if (exists("txid2geneid", envir=ChIPseekerEnv, inherits=FALSE)) {
    #     txid2geneid <- get("txid2geneid", envir=ChIPseekerEnv)
    # } else {
    #     txdb <- get("TXDB", envir=ChIPseekerEnv)
    #     txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
    #     idx <- which(sapply(txidinfo$gene_id, length) == 0)
    #     txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
    #     txid2geneid <- paste(mcols(txidinfo)[["tx_name"]],
    #                          mcols(txidinfo)[["gene_id"]],
    #                          sep="/")
    #     txid2geneid <- sub("/NA", "", txid2geneid)

    #     names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
    #     assign("txid2geneid", txid2geneid, envir=ChIPseekerEnv)
    # }
    return(as.character(txid2geneid[txid]))
}

##' Convert transcript IDs to gene IDs only
##'
##' This function converts transcript IDs directly to their corresponding gene IDs.
##'
##' @param txid character vector of transcript IDs
##' @return Character vector of gene IDs corresponding to the input transcript IDs
##' @importFrom yulab.utils get_cache_element
##' @importFrom yulab.utils update_cache_item
##' @noRd
TXID2EGID <- function(txid) {
    # ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

    txid2geneid <- get_cache_element(item = ChIPseekerCache, elements = "txid2eg")

    if(is.null(txid2geneid)){
        txdb <- get_cache_element(item = ChIPseekerCache, elements = "TXDB")
        txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
        idx <- which(sapply(txidinfo$gene_id, length) == 0)
        txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
        txid2geneid <- as.character(mcols(txidinfo)[["gene_id"]])

        names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
        update_cache_item(item = ChIPseekerCache, list("txid2eg" = txid2geneid))
    }

    # if (exists("txid2eg", envir=ChIPseekerEnv, inherits=FALSE)) {
    #     txid2geneid <- get("txid2eg", envir=ChIPseekerEnv)
    # } else {
    #     txdb <- get("TXDB", envir=ChIPseekerEnv)
    #     txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
    #     idx <- which(sapply(txidinfo$gene_id, length) == 0)
    #     txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
    #     txid2geneid <- as.character(mcols(txidinfo)[["gene_id"]])

    #     names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
    #     assign("txid2eg", txid2geneid, envir=ChIPseekerEnv)
    # }
    return(as.character(txid2geneid[txid]))
}

##' Get index of first occurrence for each unique value
##'
##' This function finds the index of the first occurrence of each unique value
##' in a vector, which is useful for handling duplicate entries in annotation results.
##'
##' @description
##' When multiple annotations map to the same feature (e.g., multiple peaks
##' mapping to the same gene), this function identifies the first occurrence
##' to avoid duplicate processing.
##'
##' @param x vector (typically character or numeric) with potentially duplicate values
##' @return Integer vector of indices corresponding to the first occurrence of
##'   each unique value in \code{x}
##' @references Based on solution by Hervé Pagès:
##'   https://support.bioconductor.org/p/70432/#70545
##' @noRd
getFirstHitIndex <- function(x) {
    ## sapply(unique(x), function(i) which(x == i)[1])
    which(!duplicated(x))
}

##' Calculate overlap matrix for sets
##'
##' This function calculates the overlap matrix for multiple sets, which is
##' useful for creating Venn diagrams and understanding set relationships.
##'
##' @description
##' The function computes the number of elements in each possible intersection
##' of the input sets. It handles all combinations of set intersections using
##' an inclusion-exclusion principle approach.
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Generates all possible combinations of set inclusion/exclusion
##'   \item For each combination, calculates the intersection size
##'   \item Uses inclusion-exclusion principle to handle overlapping intersections
##'   \item Returns a data frame with binary indicators for each set and the
##'         intersection size
##' }
##'
##' The function is generic and works with any objects that support \code{intersect()}
##' and \code{length()} methods, including GRanges objects.
##'
##' @param Sets list of objects (e.g., GRanges objects, character vectors, etc.)
##'   to calculate overlaps for. Each element should support \code{intersect()}
##'   and \code{length()} methods
##' @return A data.frame with:
##'   \itemize{
##'     \item One column per set (named after the sets)
##'     \item Binary values (0 or 1) indicating which sets are included in each
##'           intersection
##'     \item A "Weight" column containing the number of elements in each intersection
##'   }
##'   Rows represent all possible combinations of set intersections
##' @importFrom gtools permutations
##' @export
##' @author G Yu
overlap <- function(Sets) {
    ## this function is very generic.
    ## it call the getIntersectLength function to calculate
    ## the number of the intersection.
    ## if it fail, take a look at the object type were supported by getIntersectLength function.

    nn <- names(Sets)
    w <- t(apply(permutations(2,length(Sets),0:1, repeats.allowed=TRUE), 1 , rev))
    rs <- rowSums(w)
    wd <- as.data.frame(w)
    wd$n <- NA
    for (i in length(nn):0) {
        idx <- which(rs == i)
        if (i == length(nn)) {
            len <- getIntersectLength(Sets, as.logical(w[idx,]))
            wd$n[idx] <- len
        } else if (i == 0) {
            wd$n[idx] <- 0
        } else {
            for (ii in idx) {
                ##print(ii)
                len <- getIntersectLength(Sets, as.logical(w[ii,]))
                ww = w[ii,]
                jj <- which(ww == 0)
                pp <- permutations(2, length(jj), 0:1, repeats.allowed=TRUE)

                for (aa in 2:nrow(pp)) {
                    ## 1st row is all 0, abondoned
                    xx <- jj[as.logical(pp[aa,])]
                    ww[xx] =ww[xx] +1
                    bb <-  t(apply(w, 1, function(i) i == ww))
                    wd$n[rowSums(bb) == length(ww) ]
                         ww <- w[ii,]
                    len <- len - wd$n[rowSums(bb) == length(ww) ]
                    ww <- w[ii,]
                }
                wd$n[ii] <- len
            }
        }
    }
    colnames(wd) = c(names(Sets), "Weight")
    return(wd)
}


##' Calculate intersection length for selected sets
##'
##' This is a helper function that calculates the number of elements in the
##' intersection of selected sets from a list.
##'
##' @description
##' The function takes a list of sets and a logical index vector, then calculates
##' the intersection of the sets indicated by TRUE values in the index.
##'
##' @param Sets list of objects to intersect
##' @param idx logical vector of length equal to \code{Sets}, indicating which
##'   sets to include in the intersection
##' @return Integer, the number of elements in the intersection of the selected sets
##' @noRd
getIntersectLength <- function(Sets, idx) {
    ## only use intersect and length methods in this function
    ## works fine with GRanges object
    ## and easy to extend to other objects.
    ss= Sets[idx]
    ol <- ss[[1]]

    if (sum(idx) == 1) {
        return(length(ol))
    }

    for (j in 2:length(ss)) {
        ol <-  intersect(ol, ss[[j]])
    }
    return(length(ol))
}

##' Load peak data from file or GRanges object
##'
##' This function loads peak data, accepting either a file path or a GRanges object.
##'
##' @description
##' The function provides a unified interface for loading peaks, handling both
##' file paths and GRanges objects. If a file path is provided, it calls
##' \code{readPeakFile()} to load the data.
##'
##' @param peak character, path to peak file (BED format), or GRanges object
##'   containing peaks
##' @param verbose logical, whether to print progress messages. Default is FALSE
##' @return GRanges object containing the peaks
##' @seealso \code{\link{readPeakFile}} for reading peak files
##' @noRd
loadPeak <- function(peak, verbose=FALSE) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        if (verbose)
            cat(">> loading peak file...\t\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
        peak.gr <- readPeakFile(peak, as="GRanges")
    } else {
        stop("peak should be GRanges object or a peak file...")
    }
    return(peak.gr)
}

##' Get TxDb object with default fallback
##'
##' This function returns the provided TxDb object, or uses a default human
##' genome annotation if NULL is provided. It does NOT load from files.
##'
##' @description
##' The function simply returns the provided TxDb/EnsDb object if it is not NULL.
##' If NULL is provided, it uses \code{TxDb.Hsapiens.UCSC.hg19.knownGene} as a
##' default (with a warning).
##'
##' \strong{Note:} This function does not load TxDb objects from files (sqlite,
##' GTF, GFF, etc.). The TxDb object must already be loaded in R. To create a
##' TxDb object from a file, use functions like:
##' \itemize{
##'   \item \code{GenomicFeatures::makeTxDbFromGFF()} for GTF/GFF files
##'   \item \code{GenomicFeatures::makeTxDbFromUCSC()} for UCSC genomes
##'   \item \code{AnnotationDbi::loadDb()} for sqlite database files
##'   \item Pre-built annotation packages (e.g., \code{TxDb.Hsapiens.UCSC.hg19.knownGene})
##' }
##'
##' @param TxDb TxDb or EnsDb annotation object (already loaded in R), or NULL
##'   to use default. Must be an already-loaded R object, not a file path
##' @return TxDb or EnsDb object. If input is NULL, returns the default human
##'   genome annotation (\code{TxDb.Hsapiens.UCSC.hg19.knownGene})
##' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
##' @noRd
loadTxDb <- function(TxDb) {
    if ( is.null(TxDb) ) {
        warning(">> TxDb is not specified, use 'TxDb.Hsapiens.UCSC.hg19.knownGene' by default...")
        TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    return(TxDb)
}

##' Get gene or transcript features from TxDb
##'
##' This function retrieves gene or transcript features from a TxDb object,
##' using caching to improve performance.
##'
##' @description
##' The function extracts gene or transcript features from the TxDb object and
##' caches them for future use. It initializes the ChIPseeker cache environment
##' and retrieves features from cache if available, or computes them if not.
##'
##' @param TxDb TxDb or EnsDb annotation object
##' @param by character, one of "gene" or "transcript". Determines whether to
##'   return gene-level or transcript-level features. Default is "gene"
##' @return GRanges object containing gene or transcript features:
##'   \itemize{
##'     \item If \code{by="gene"}: Returns genes with gene_id metadata
##'     \item If \code{by="transcript"}: Returns transcripts with transcript
##'           information
##'   }
##' @importFrom AnnotationDbi get
##' @importFrom GenomicFeatures genes
##' @importFrom GenomicFeatures transcriptsBy
##' @importFrom yulab.utils get_cache_element
##' @importFrom yulab.utils update_cache_item
##' @noRd
getGene <- function(TxDb, by="gene") {
    .ChIPseekerEnv(TxDb, item = ChIPseekerCache)
    # ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

    by <- match.arg(by, c("gene", "transcript"))

    if (by == "gene") {

        features <- get_cache_element(item = ChIPseekerCache, elements = "Genes")

        if(is.null(features)){
            features <- suppressMessages(genes(TxDb))
            update_cache_item(item = ChIPseekerCache, list("Genes" = features))
        }

        # if ( exists("Genes", envir=ChIPseekerEnv, inherits=FALSE) ) {
        #     features <- get("Genes", envir=ChIPseekerEnv)
        # } else {
        #     features <- suppressMessages(genes(TxDb))
        #     assign("Genes", features, envir=ChIPseekerEnv)
        # }
    } else {

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
    }

    return(features)
}


##' Get sample peak files included in the package
##'
##' This function returns a list of sample ChIP-seq peak files included in the
##' package for demonstration and testing purposes.
##'
##' @description
##' The function scans the package's extdata directory for sample peak files
##' and returns their paths. The files are named based on the protein/condition
##' they represent (extracted from the filename).
##'
##' @return A named list of file paths. Names are extracted from filenames
##'   (typically protein/condition names), and values are full file paths to the
##'   sample BED files
##' @export
##' @author G Yu
getSampleFiles <- function() {
    dir <- system.file("extdata", "GEO_sample_data", package="ChIPseeker")
    files <- list.files(dir)
    ## protein <- sub("GSM\\d+_", "", files)
    ## protein <- sub("_.+", "", protein)
    protein <- gsub(pattern='GSM\\d+_(\\w+_\\w+)_.*', replacement='\\1',files)
    protein <- sub("_Chip.+", "", protein)
    res <- paste(dir, files, sep="/")
    res <- as.list(res)
    names(res) <- protein
    return(res)
}
## @importFrom RCurl getURL
## getDirListing <- function (url) {
##     ## from GEOquery
##     print(url)
##     a <- getURL(url)
##     b <- textConnection(a)
##     d <- read.table(b, header = FALSE)
##     close(b)
##     return(d)
## }


##' Check if a path is a directory
##'
##' @param dir character, path to check
##' @return Logical, TRUE if the path exists and is a directory, FALSE otherwise
##' @noRd
is.dir <- function(dir) {
    if (file.exists(dir) == FALSE)
        return(FALSE)
    return(file.info(dir)$isdir)
}


##' Parse target peak parameter
##'
##' This function parses the targetPeak parameter, handling both file paths and
##' directories containing peak files.
##'
##' @description
##' The function accepts either:
##' \itemize{
##'   \item A single directory path: scans for BED files in the directory
##'   \item A single file path: returns the file path if it exists
##'   \item Multiple file paths: returns existing files
##' }
##'
##' @param targetPeak character vector, either a directory path, a single file
##'   path, or multiple file paths
##' @return Character vector of file paths to peak files
##' @noRd
parse_targetPeak_Param <- function(targetPeak) {
    if (length(targetPeak) == 1) {
        if (is.dir(targetPeak)) {
            files <- list.files(path=targetPeak)
            idx <- unlist(sapply(c("bed", "bedGraph", "Peak"), grep, x=files))
            idx <- sort(unique(idx))
            files <- files[idx]
            targetPeak <- sub("/$", "", targetPeak)
            res <- paste(targetPeak, files, sep="/")
        } else {
            if (!file.exists(targetPeak)) {
                stop("bed file is not exists...")
            } else {
                res <- targetPeak
            }
        }
    } else {
        if (is.dir(targetPeak[1])) {
            stop("targetPeak should be a vector of bed file names or a folder containing bed files...")
        } else {
            res <- targetPeak[file.exists(targetPeak)]
            if (length(res) == 0) {
                stop("targetPeak file not exists...")
            }
        }
    }
    return(res)
}


##' Get gene ID type from TxDb metadata
##'
##' This function extracts the gene ID type from TxDb metadata, which is used
##' to determine how to query annotation databases.
##'
##' @description
##' The function searches the TxDb metadata for the "Type of Gene ID" field and
##' returns its value. This is used to determine whether gene IDs are Entrez IDs,
##' Ensembl IDs, or other types.
##'
##' @param TxDb TxDb or EnsDb annotation object
##' @return Character string indicating the gene ID type (e.g., "Entrez Gene ID",
##'   "Ensembl Gene ID")
##' @noRd
IDType <- function(TxDb) {
    ##
    ## IDType <- metadata(TxDb)[8,2]
    ##
    ## update: 2015-10-27
    ## now IDType change from metadata(TxDb)[8,2] to metadata(TxDb)[9,2]
    ## it may change in future too
    ##
    ## it's safe to extract via grep

    md <- metadata(TxDb)
    md[grep("Type of Gene ID", md[,1]), 2]
}

##' Convert list of data frames to a single data frame
##'
##' This function combines a list of data frames into a single data frame,
##' adding an identifier column to track the source of each row.
##'
##' @description
##' The function:
##' \itemize{
##'   \item Adds a \code{.id} column containing the list element names
##'   \item Handles missing columns by adding NA values
##'   \item Combines all data frames using \code{rbind()}
##'   \item Sets \code{.id} as a factor with levels in reverse order of list names
##' }
##'
##' @param dataList list of data frames to combine. If named, names are used
##'   as identifiers
##' @return A data.frame combining all input data frames with an additional
##'   \code{.id} column identifying the source
##' @noRd
list_to_dataframe <- function(dataList) {
    if (is.null(names(dataList)))
        return(do.call('rbind', dataList))

    cn <- lapply(dataList, colnames) %>% unlist %>% unique
    cn <- c('.id', cn)
    dataList2 <- lapply(seq_along(dataList), function(i) {
        data = dataList[[i]]
        data$.id = names(dataList)[i]
        idx <- ! cn %in% colnames(data)
        if (sum(idx) > 0) {
            for (i in cn[idx]) {
                data[, i] <- NA
            }
        }
        return(data[,cn])
    })
    res <- do.call('rbind', dataList2)
    res$.id <- factor(res$.id, levels=rev(names(dataList)))
    return(res)
}

##' @importFrom GenomicRanges GRangesList
##' @export
GenomicRanges::GRangesList

## . function was from plyr package
##' capture name of variable
##'
##' @rdname dotFun
##' @export
##' @title .
##' @param ... expression
##' @param .env environment
##' @return expression
##' @examples
##' x <- 1
##' eval(.(x)[[1]])
. <- function (..., .env = parent.frame()) {
    structure(as.list(match.call()[-1]), env = .env, class = "quoted")
}


##' Check validity of upstream and downstream parameters
##'
##' This function validates that upstream and downstream parameters are properly
##' specified and compatible for use in genomic range calculations.
##'
##' @description
##' The function performs several validation checks:
##' \itemize{
##'   \item Ensures upstream and downstream have the same type (both numeric or
##'         both rel objects)
##'   \item Verifies they are numeric or NULL (not other types)
##'   \item For rel objects, checks values are in (0, 1)
##'   \item For numeric values, checks they are >= 1
##' }
##'
##' @param upstream numeric or rel object, upstream distance parameter
##' @param downstream numeric or rel object, downstream distance parameter
##' @return Invisibly returns NULL if validation passes, otherwise stops with
##'   an error message
##' @importFrom ggplot2 rel
##' @noRd
check_upstream_and_downstream <- function(upstream, downstream){

    ## upstream and downstream should be the same type
    if(class(upstream) != class(downstream)){
        stop("the type of upstream and downstream should be the same...")
    }

    ## downstream and upstream parameter should be numeric or NULL
    if(!is.numeric(upstream) && !is.null(upstream)){
        stop("upstream and downstream parameter should be numeric or NULL...")
    }

    ## the value of rel object should be in (0,1)
    if(inherits(upstream, 'rel')){
        if(as.numeric(upstream) < 0 || as.numeric(upstream) >1 ){
            stop('the value of rel object should be in (0,1)...')
        }
    }

    ## check actual number
    if(is.numeric(upstream) && !inherits(upstream, 'rel')){
        if(upstream < 1 | downstream < 1){
            stop('if upstream or downstream is integer, the value of it should be greater than 1...')
        }
    }
}


##' @importFrom ggplot2 rel
##'
##' @export
ggplot2::rel


##' Generate labels for genomic feature plots
##'
##' This function generates appropriate labels for plotting genomic features,
##' handling different feature types and plot types.
##'
##' @description
##' The function creates labels based on the feature type and plot type:
##' \itemize{
##'   \item For gene/transcript/exon/intron: uses "TSS" and "TTS" (transcription
##'         start/termination site)
##'   \item For UTRs: uses feature-specific labels (e.g., "5UTR_SS", "3UTR_TS")
##'   \item For body plots: returns both start and end site labels
##'   \item For start_site or end_site: returns single label
##' }
##'
##' @param type character, one of "start_site", "end_site", or "body"
##' @param by character, feature type. One of 'gene', 'transcript', 'exon',
##'   'intron', '3UTR', '5UTR', or 'UTR'
##' @return Character vector of labels. For "body" type, returns a vector of
##'   length 2 (start and end labels). Otherwise returns a single label
##' @noRd
make_label <- function(type, by){

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

    return(label)
}

##' Get genome version from cached TxDb
##'
##' This function extracts the genome version information from the cached TxDb
##' object metadata.
##'
##' @description
##' The function retrieves the genome version (e.g., "hg19", "hg38", "mm10")
##' from the TxDb metadata stored in the cache. This is used for display
##' purposes to show which genome is being used.
##'
##' @return Character string containing the genome version (e.g., "hg19")
##' @importFrom yulab.utils get_cache_item
##' @noRd
get_env_genome <- function(){

    current_env <- get_cache_item(item = ChIPseekerCache)

    env_txdb <- current_env$TXDB
    env_txdb_meta <- S4Vectors::metadata(env_txdb)
    env_txdb_version <- env_txdb_meta[grep("Genome",env_txdb_meta[,1]),2]

    return(env_txdb_version)
}