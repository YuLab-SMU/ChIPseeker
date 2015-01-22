
##' @importFrom plyr ddply
##' @importFrom IRanges ranges
##' @importFrom IRanges ranges<-
##' @importFrom IRanges start
##' @importFrom IRanges start<-
##' @importFrom IRanges end
##' @importFrom IRanges end<-
getAllFlankingGene <- function(peak.gr, features, distance=5000) {
    peak.gr2 <- peak.gr
    start(ranges(peak.gr)) = start(ranges(peak.gr)) - distance
    end(ranges(peak.gr)) = end(ranges(peak.gr)) + distance
    hit <- findOverlaps(peak.gr, features)
    qh <- queryHits(hit)
    sh <- subjectHits(hit)
    featureHit <- features[sh]
    eg <- TXID2EG(featureHit$tx_id, geneIdOnly=TRUE)
    names(featureHit)=NULL

    hitInfo <- as.data.frame(featureHit)
    hitInfo$geneId <- eg
    hitInfo$peakIdx <- qh

    overlapHit <- findOverlaps(peak.gr2, featureHit)
    hitInfo$distance <- NA
    hitInfo$distance[subjectHits(overlapHit)] <- 0

    psD <- ifelse(strand(featureHit) == "+",
                  start(peak.gr2[qh]) - start(featureHit),
                  end(featureHit)-end(peak.gr2[qh]))
    
    peD <- ifelse(strand(featureHit) == "+",
                  end(peak.gr2[qh]) - start(featureHit),
                  end(featureHit)-start(peak.gr2[qh]))

    idx <- abs(psD) > abs(peD)
    dd <- psD
    dd[idx] <- peD[idx]

    ii <- is.na(hitInfo$distance)
    hitInfo$distance[ii] <- dd[ii]

    peakIdx <- tx_name <- geneId <- distance <- NULL
    hitInfo2 <- ddply(hitInfo, .(peakIdx), transform,
                      flank_txIds=paste(tx_name, collapse=";"),
                      flank_geneIds=paste(geneId, collapse=";"),
                      flank_gene_distances=paste(distance, collapse=";"))

    res <- hitInfo2[,c("peakIdx", "flank_txIds", "flank_geneIds", "flank_gene_distances")]
    res <- unique(res)

    res$flank_txIds <- as.character(res$flank_txIds)
    res$flank_geneIds <- as.character(res$flank_geneIds)
    res$flank_gene_distances <- as.character(res$flank_gene_distances)
    
    return(res)
    
}
