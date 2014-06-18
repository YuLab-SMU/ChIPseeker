##' get index of features that closest to peak and calculate distance
##'
##' 
##' @title getNearestFeatureIndicesAndDistances 
##' @param peaks peak in GRanges 
##' @param features features in GRanges
##' @return data.frame
##' @importFrom IRanges precede
##' @importFrom IRanges follow
##' @importFrom IRanges start
##' @importFrom IRanges end
##' @importFrom BiocGenerics strand
## @importMethodsFrom GenomicRanges strand
##' @author G Yu
getNearestFeatureIndicesAndDistances <- function(peaks, features) {
    ## peaks only conatin all peak records, in GRanges object
    ## feature is the annotation in GRanges object

    ## only keep start position based on strand
    start(features) <- end(features) <- ifelse(strand(features) == "+", start(features), end(features))
    
    ## nearest from peak start
    ps.idx <- precede(peaks, features)
    ## nearest from peak end
    pe.idx <- follow(peaks, features)
    
    ## features from nearest peak start
    psF <- features[ps.idx]
    ## feature distances from peak start
    psD <- ifelse(strand(psF) == "+",
                 start(peaks) - start(psF),
                 end(psF)-end(peaks))

    ## features from nearest peak end
    peF <- features[pe.idx]
    ## feature distances from peak end
    peD <- ifelse(strand(peF) == "+",
                  end(peaks) - start(peF),
                  end(peF)-start(peaks))

    pse <- data.frame(ps=psD, pe=peD)
    j <- apply(pse, 1, function(i) which.min(abs(i)))

    ## index
    idx <- ps.idx
    idx[j==2] <- pe.idx[j==2]
    
    ## distance
    dd <- psD
    dd[j==2] <- peD[j==2]

    
    hit <- findOverlaps(peaks, features)
    if ( length(hit) != 0 ) {
        qh <- queryHits(hit)
        hit.idx <- getFirstHitIndex(qh)
        hit <- hit[hit.idx]
        peakIdx <- queryHits(hit)
        featureIdx <- subjectHits(hit)

        idx[peakIdx] <- featureIdx
        dd[peakIdx] <- 0
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
    
    res <- data.frame(index=idx, distance=dd)
    return(res)
}

isPeakFeatureOverlap <- function(peak, feature) {
    peakRange <- ranges(peak)
    featureRange <- ranges(feature)
    x <- intersect(peakRange, featureRange)
    return(length(x) != 0)
}
