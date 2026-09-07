##' Get all flanking genes for peaks
##'
##' This function identifies all genes or transcripts that are within a specified
##' distance from each peak, extending the peak range on both sides to search for
##' nearby features.
##'
##' @description
##' For each peak, the function extends the peak range by the specified distance
##' on both upstream and downstream sides, then identifies all genes/transcripts
##' that overlap with this extended region. It calculates the distance from the
##' original peak to each flanking feature and aggregates the results by peak.
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Extends each peak by \code{distance} base pairs on both sides
##'   \item Finds all features (genes/transcripts) that overlap with the extended
##'         peak regions using \code{findOverlaps}
##'   \item For features that overlap the original peak, sets distance to 0
##'   \item For non-overlapping features, calculates the distance from the peak to
##'         the feature's TSS (transcription start site), using the closer of peak
##'         start or peak end
##'   \item Groups results by peak and concatenates all flanking gene IDs and
##'         distances with semicolons
##' }
##'
##' Distance calculation is strand-aware and uses the closer peak end as reference:
##' \itemize{
##'   \item For positive strand features: calculates distance from both peak start
##'         and peak end to the feature's TSS (start position), then selects the
##'         one with smaller absolute value
##'   \item For negative strand features: calculates distance from both peak start
##'         and peak end to the feature's TSS (end position), then selects the
##'         one with smaller absolute value
##' }
##' This ensures the reported distance represents the shortest distance from any
##' point of the peak to the feature's transcription start site.
##'
##' @param peak.gr GRanges object containing genomic ranges of peaks
##' @param features GRanges object containing genomic features (genes or
##'   transcripts) to search for flanking genes. Typically obtained from
##'   \code{getGene(TxDb)} or similar functions
##' @param level character, one of "transcript" or "gene". Determines whether to
##'   work at transcript level or gene level. When "transcript", transcript IDs
##'   are included in the output. Default is "transcript"
##' @param distance numeric, distance in base pairs to extend on both sides of
##'   each peak when searching for flanking genes. Default is 5000 (5kb)
##' @return A data.frame with the following columns:
##'   \itemize{
##'     \item \code{peakIdx}: Index of the peak (1-based)
##'     \item \code{flank_geneIds}: Semicolon-separated list of all flanking
##'       gene IDs within the extended region
##'     \item \code{flank_gene_distances}: Semicolon-separated list of distances
##'       from the peak to each flanking gene. Distance of 0 indicates overlap
##'       with the original peak. Positive values indicate downstream, negative
##'       values indicate upstream
##'     \item \code{flank_txIds}: (Only when \code{level="transcript"})
##'       Semicolon-separated list of all flanking transcript IDs
##'   }
##'   The data.frame contains one row per peak that has at least one flanking gene.
##' @import IRanges
##' @importFrom dplyr mutate
##' @importFrom dplyr group_by
##' @author G Yu
##' @seealso \code{\link{annotatePeak}} which uses this function when
##'   \code{addFlankGeneInfo=TRUE}
##' @noRd
getAllFlankingGene <- function(peak.gr, features,
                              level="transcript", distance=5000) {
    peak.gr2 <- peak.gr

    ### extend the peak range by the distance
    start(ranges(peak.gr)) = start(ranges(peak.gr)) - distance
    end(ranges(peak.gr)) = end(ranges(peak.gr)) + distance

    ### find overlap between the extended peak range and features
    hit <- findOverlaps(peak.gr, unstrand(features))
    qh <- queryHits(hit)
    sh <- subjectHits(hit)

    featureHit <- features[sh]
    names(featureHit)=NULL
    hitInfo <- as.data.frame(featureHit)

    ### get the gene ID of the features
    if (level == "transcript") {
        eg <- TXID2EG(featureHit$tx_id, geneIdOnly=TRUE)
        hitInfo$geneId <- eg
    } else {
        cn <- colnames(hitInfo)
        colnames(hitInfo)[cn == "gene_id"] <- "geneId"
    }


    hitInfo$peakIdx <- qh

    ### find overlap between the peak range and the features
    ### Features that overlap the original peak get distance = 0

    overlapHit <- findOverlaps(peak.gr2, unstrand(featureHit))
    hitInfo$distance <- NA
    hitInfo$distance[subjectHits(overlapHit)] <- 0

    ### calculate the distance between the non-overlapping peak and the features
    psD <- ifelse(strand(featureHit) == "+",
                  start(peak.gr2[qh]) - start(featureHit),
                  end(featureHit)-end(peak.gr2[qh]))

    ### calculate the distance between the non-overlapping peak and the features
    peD <- ifelse(strand(featureHit) == "+",
                  end(peak.gr2[qh]) - start(featureHit),
                  end(featureHit)-start(peak.gr2[qh]))

    idx <- abs(psD) > abs(peD)
    dd <- psD
    dd[idx] <- peD[idx]

    ii <- is.na(hitInfo$distance)
    hitInfo$distance[ii] <- dd[ii]

    peakIdx <- tx_name <- geneId <- distance <- NULL

    if (level == "transcript") {
        hitInfo2 <- group_by(hitInfo, peakIdx) %>%
            mutate(flank_txIds=paste(tx_name, collapse=";"),
                   flank_geneIds=paste(geneId, collapse=";"),
                   flank_gene_distances=paste(distance, collapse=";"))
        res <- hitInfo2[,c("peakIdx", "flank_txIds", "flank_geneIds", "flank_gene_distances")]
        res$flank_txIds <- as.character(res$flank_txIds)
    } else {
        hitInfo2 <- group_by(hitInfo, peakIdx) %>%
            mutate(flank_geneIds=paste(geneId, collapse=";"),
                   flank_gene_distances=paste(distance, collapse=";"))
        res <- hitInfo2[,c("peakIdx", "flank_geneIds", "flank_gene_distances")]
    }

    res <- unique(res)
    res$flank_geneIds <- as.character(res$flank_geneIds)
    res$flank_gene_distances <- as.character(res$flank_gene_distances)

    return(res)
}
