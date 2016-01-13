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

    if (!ignoreOverlap) {
        overlap_hit_TSS <- findOverlaps(peaks, unstrand(features))
    }
    
    if (sameStrand) {
        ## nearest from peak start
        ps.idx <- follow(peaks, features)
        
        ## nearest from peak end
        pe.idx <- precede(peaks, features)
    } else {
        ps.idx <- follow(peaks, unstrand(features))
        pe.idx <- precede(peaks, unstrand(features))
    }
    
    na.idx <- is.na(ps.idx) | is.na(pe.idx)
    ## if (sum(na.idx) > 1) {
    if (sum(na.idx) > 0) { ## suggested by Thomas Schwarzl
        ps.idx <- ps.idx[!na.idx]
        pe.idx <- pe.idx[!na.idx]
        peaks <- peaks[!na.idx]
    }
    
    ## features from nearest peak start
    psF <- features[ps.idx]
    ## feature distances from peak start
    psD <- ifelse(strand(psF) == "+", 1, -1) *
        (start(peaks) - start(psF))
    
    ## features from nearest peak end
    peF <- features[pe.idx]
    ## feature distances from peak end
    peD <- ifelse(strand(peF) == "+", 1, -1) *
        (end(peaks) - start(peF))
    
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
                
                idx[peakIdx] <- featureIdx
                distance_both_end <- data.frame(start=start(peaks) - start(features[featureIdx]),
                                          end = end(peaks) - start(features[featureIdx]))
                distance_idx <- apply(distance_both_end, 1, function(i) which.min(abs(i)))
                distance_minimal <- distance_both_end[,1]
                distance_minimal[distance_idx == 2] <- distance_both_end[distance_idx==2, 2]

                dd[peakIdx] <- distance_minimal
            }
        }
        
        hit <- overlap_hit_TSS
        if ( length(hit) != 0 ) {
            qh <- queryHits(hit)
            hit.idx <- getFirstHitIndex(qh)
            hit <- hit[hit.idx]
            peakIdx <- queryHits(hit)
            featureIdx <- subjectHits(hit)
            
            idx[peakIdx] <- featureIdx
            dd[peakIdx] <- 0
        }
        
    }
    
    ## pn.idx <- nearest(peaks, features)
    ## isOverlap <- sapply(1:length(pn.idx), function(i) {
    ##     isPeakFeatureOverlap(peaks[i], features2[pn.idx[i]])
    ## })
    ## isOverlap <- unlist(isOverlap)

    ## if(sum(isOverlap) > 0) {
    ##     idx[isOverlap] <- pn.idx[isOverlap]
    ##     dd[isOverlap] <- 0
    ## }
    
    res <- list(index=idx, distance=dd, peak=peaks)
    
    return(res)
}

isPeakFeatureOverlap <- function(peak, feature) {
    peakRange <- ranges(peak)
    featureRange <- ranges(feature)
    x <- intersect(peakRange, featureRange)
    return(length(x) != 0)
}
