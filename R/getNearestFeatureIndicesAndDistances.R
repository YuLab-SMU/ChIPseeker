##' Get index of nearest features to peaks and calculate distances
##'
##' This function identifies the nearest gene/transcript feature to each peak and
##' calculates the distance from the peak to the feature's transcription start site (TSS).
##'
##' @description
##' The function finds the closest feature to each peak by considering both upstream
##' and downstream directions. It calculates distances from both the peak start and
##' peak end to the feature's TSS, then selects the feature with the minimum absolute
##' distance. The function handles overlaps between peaks and features, with different
##' behaviors depending on the \code{overlap} parameter.
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Resizes features to width=1 (TSS points) for efficient distance calculation
##'   \item Finds nearest upstream features using \code{follow()} and nearest
##'         downstream features using \code{precede()}
##'   \item Calculates distances from both peak start and peak end to each feature's TSS
##'   \item Selects the feature with the minimum absolute distance (or enforces
##'         direction based on \code{ignoreUpstream}/\code{ignoreDownstream})
##'   \item If overlaps are allowed (\code{ignoreOverlap=FALSE}), identifies
##'         overlapping features and assigns them with distance=0 (for \code{overlap="TSS"})
##'         or calculates distance from the closer peak end (for \code{overlap="all"})
##'   \item Returns indices, distances, and filtered peaks (excluding peaks with no
##'         nearest feature)
##' }
##'
##' Distance calculation is strand-aware: positive values indicate downstream,
##' negative values indicate upstream. Overlapping peaks get distance=0 if overlap is "TSS".
##'
##' @param peaks GRanges object containing genomic ranges of peaks
##' @param features GRanges object containing genomic features (genes or
##'   transcripts) to search for nearest features. Typically obtained from
##'   \code{getGene(TxDb)} or similar functions
##' @param sameStrand logical, whether to only consider features on the same
##'   strand as the peak when finding nearest features. If FALSE, searches both
##'   strands. Default is FALSE
##' @param ignoreOverlap logical, whether to ignore overlaps between peaks and
##'   features when finding the nearest feature. If FALSE, overlapping features
##'   will be prioritized and assigned distance=0 (for TSS overlaps) or calculated
##'   distance (for full feature overlaps). Default is FALSE
##' @param ignoreUpstream logical, if TRUE, only considers features downstream
##'   (3' end) of the peak. This restricts annotation to genes that come after the
##'   peak. Default is FALSE
##' @param ignoreDownstream logical, if TRUE, only considers features upstream
##'   (5' end) of the peak. This restricts annotation to genes that come before
##'   the peak. Default is FALSE
##' @param overlap character, one of "TSS" or "all". Determines how overlaps are
##'   detected and handled:
##'   \itemize{
##'     \item "TSS": Only considers overlaps with the TSS point (after features
##'           are resized to width=1). Overlapping peaks get distance=0
##'     \item "all": Considers overlaps with any part of the full feature range
##'           (before resizing). For overlapping peaks, calculates distance from
##'           the closer peak end (start or end) to the feature's TSS
##'   }
##'   Default is "TSS"
##' @return A list with three components:
##'   \itemize{
##'     \item \code{index}: Integer vector of feature indices (1-based) for the
##'       nearest feature to each peak. Length equals the number of peaks with
##'       valid nearest features
##'     \item \code{distance}: Numeric vector of distances from each peak to the
##'       TSS of its nearest feature. Positive values indicate downstream, negative
##'       values indicate upstream, and 0 indicates overlap. Length equals the
##'       number of peaks with valid nearest features
##'     \item \code{peak}: GRanges object containing only the peaks that have
##'       valid nearest features (peaks with no nearest feature in either direction
##'       are excluded)
##'   }
##' @import BiocGenerics IRanges GenomicRanges
##' @author G Yu
getNearestFeatureIndicesAndDistances <- function(peaks, features,
                                                 sameStrand = FALSE,
                                                 ignoreOverlap=FALSE,
                                                 ignoreUpstream=FALSE,
                                                 ignoreDownstream=FALSE,
                                                 overlap = "TSS") {

    overlap <- match.arg(overlap, c("TSS", "all"))

   ### find overlap between peaks and features
    if (!ignoreOverlap && overlap == "all") {
        overlap_hit <- findOverlaps(peaks, unstrand(features))
    }

    ## peaks only conatin all peak records, in GRanges object
    ## feature is the annotation in GRanges object

    ## only keep start position based on strand
    ## start(features) <- end(features) <- ifelse(strand(features) == "+", start(features), end(features))
    features <- resize(features, width=1) # faster

    ## add dummy NA feature for peaks that are at the last or first feature
    ## suggested by Michael Kluge
    features.bak <- features
    seqlevels(features) <- c(seqlevels(features), "chrNA")
    dummy <- GRanges("chrNA", IRanges(1,1))

    ## dummy$tx_id <- -1
    ## dummy$tx_name <- "NA"

    cns <- names(mcols(features))
    for (cn in cns) {
        if (grepl('id', cn)) {
            mcols(dummy)[[cn]] <- -1
        } else {
            mcols(dummy)[[cn]] <- NA
        }
    }

    features <- append(features, dummy)
    dummyID <- length(features)


    ### nearest upstream and downstream features, but ignore overlap with peak
    if (sameStrand) {
        ## nearest from peak start
        ps.idx <- follow(peaks, features)

        ## nearest from peak end
        pe.idx <- precede(peaks, features)
    } else {
        ps.idx <- follow(peaks, unstrand(features))
        pe.idx <- precede(peaks, unstrand(features))
    }

    na.idx <- is.na(ps.idx) & is.na(pe.idx)
    if (sum(na.idx) > 0) { ## suggested by Thomas Schwarzl
        ps.idx <- ps.idx[!na.idx]
        pe.idx <- pe.idx[!na.idx]
        ##peaks <- peaks[!na.idx]
    }

    # set NA values to dummy value if only one entry is affected
    ps.idx[is.na(ps.idx)] <- dummyID
    pe.idx[is.na(pe.idx)] <- dummyID

    ## features from nearest peak start
    psF <- features[ps.idx]

    ## feature distances from peak start
    psD <- ifelse(strand(psF) == "+", 1, -1) *
        (start(peaks[!na.idx]) - start(psF))
    psD[ps.idx == dummyID] <- Inf # ensure that there is even no match if a seq with name "chrNA" exists

    ## features from nearest peak end
    peF <- features[pe.idx]
    ## feature distances from peak end
    peD <- ifelse(strand(peF) == "+", 1, -1) *
        (end(peaks[!na.idx]) - start(peF))
    peD[pe.idx == dummyID] <- Inf # ensure that there is even no match if a seq with name "chrNA" exists

    ## restore the old feature object
    features <- features.bak

    pse <- data.frame(ps=psD, pe=peD)
    if (ignoreUpstream) {
        j <- rep(2, nrow(pse))
    } else if (ignoreDownstream) {
        j <- rep(1, nrow(pse))
    } else {
        j <- apply(pse, 1, function(i) which.min(abs(i)))
    }

    ## index
    idx <- ps.idx
    idx[j==2] <- pe.idx[j==2]

    ## distance
    dd <- psD
    dd[j==2] <- peD[j==2]

    index <- distanceToTSS <- rep(NA, length(peaks))
    distanceToTSS[!na.idx] <- dd
    index[!na.idx] <- idx

    if (!ignoreOverlap) {
        if (overlap == "all") {
            hit <- overlap_hit
            if ( length(hit) != 0 ) {
                qh <- queryHits(hit)
                hit.idx <- getFirstHitIndex(qh)
                hit <- hit[hit.idx]
                peakIdx <- queryHits(hit)
                featureIdx <- subjectHits(hit)

                index[peakIdx] <- featureIdx
                distance_both_end <- data.frame(start=start(peaks[peakIdx]) - start(features[featureIdx]),
                                          end = end(peaks[peakIdx]) - start(features[featureIdx]))
                distance_idx <- apply(distance_both_end, 1, function(i) which.min(abs(i)))
                distance_minimal <- distance_both_end[,1]
                distance_minimal[distance_idx == 2] <- distance_both_end[distance_idx==2, 2]

                distanceToTSS[peakIdx] <- distance_minimal * ifelse(strand(features[featureIdx]) == "+", 1, -1)

            }
        } else {
            ## overlap == "TSS": find overlaps with TSS points (resized features of width 1, TSS sites only)
            hit <- findOverlaps(peaks, unstrand(features))

            if ( length(hit) != 0 ) {
                qh <- queryHits(hit)
                hit.idx <- getFirstHitIndex(qh)
                hit <- hit[hit.idx]
                peakIdx <- queryHits(hit)
                featureIdx <- subjectHits(hit)

                index[peakIdx] <- featureIdx
                distanceToTSS[peakIdx] <- 0
            }
        }

    }

    j <- is.na(distanceToTSS) | is.na(index)

    res <- list(index=index[!j],
                distance=distanceToTSS[!j],
                peak=peaks[!j])

    return(res)
}

isPeakFeatureOverlap <- function(peak, feature) {
    peakRange <- ranges(peak)
    featureRange <- ranges(feature)
    x <- intersect(peakRange, featureRange)
    return(length(x) != 0)
}
