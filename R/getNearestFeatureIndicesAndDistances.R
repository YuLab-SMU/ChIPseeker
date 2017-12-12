##' get index of features that closest to peak and calculate distance
##'
##'
##' @title getNearestFeatureIndicesAndDistances
##' @param peaks peak in GRanges
##' @param features features in GRanges
##' @param sameStrand logical, whether find nearest gene in the same strand
##' @param ignoreOverlap logical, whether ignore overlap of TSS with peak
##' @param ignoreUpstream logical, if True only annotate gene at the 3' of the peak.
##' @param ignoreDownstream logical, if True only annotate gene at the 5' of the peak.
##' @param overlap one of "TSS" or "all"
##' @return list
##' @import BiocGenerics IRanges GenomicRanges
##' @author G Yu
getNearestFeatureIndicesAndDistances <- function(peaks, features,
                                                 sameStrand = FALSE,
                                                 ignoreOverlap=FALSE,
                                                 ignoreUpstream=FALSE,
                                                 ignoreDownstream=FALSE,
                                                 overlap = "TSS") {

    overlap <- match.arg(overlap, c("TSS", "all"))

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
        ## hit <- findOverlaps(peaks, unstrand(features))

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
        }

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
