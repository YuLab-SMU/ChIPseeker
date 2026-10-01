##' Align the reported feature with the feature that supplied the annotation
##'
##' `annotationFeatureId` holds the id (transcript or gene, depending on
##' `idColumn`) of the region that produced the exon/intron/UTR annotation, see
##' `getGenomicAnnotation.internal()`.  The index and the distance of those
##' entries are moved to that feature and the distance to its TSS is
##' recalculated, so that annotation and the reported feature columns describe
##' the same feature (issue #252).
##'
##' @param peak.gr GRanges of the annotated peaks
##' @param features GRanges of the features reported by annotatePeak()
##' @param index integer vector of nearest feature indices
##' @param distance numeric vector of distances to the nearest TSS
##' @param annotationFeatureId character vector of feature ids, `NA` where the
##'   annotation does not come from an exon/intron/UTR hit
##' @param idColumn character, name of the mcols column of `features` holding
##'   the feature ids (`"tx_id"` or `"gene_id"`)
##' @return list with the aligned `index` and `distance`
##' @noRd
##' @author G Yu
.alignAnnotationFeature <- function(peak.gr, features, index, distance,
                                    annotationFeatureId, idColumn) {
    if (is.null(annotationFeatureId) ||
        !idColumn %in% colnames(mcols(features))) {
        return(list(index = index, distance = distance))
    }

    featureIdx <- match(
        as.character(annotationFeatureId),
        as.character(mcols(features)[[idColumn]])
    )
    align_idx <- which(!is.na(featureIdx))
    if (length(align_idx) == 0) {
        return(list(index = index, distance = distance))
    }

    index[align_idx] <- featureIdx[align_idx]

    ## recompute the distance to the TSS of the aligned feature, using the
    ## same strand-aware convention as
    ## getNearestFeatureIndicesAndDistances(): positive values are downstream.
    selected <- features[featureIdx[align_idx]]
    selected_tss <- resize(selected, width = 1)
    selected_peaks <- peak.gr[align_idx]
    selected_strand <- as.character(strand(selected))
    d_start <- ifelse(
        selected_strand == "+",
        start(selected_peaks) - start(selected_tss),
        start(selected_tss) - start(selected_peaks)
    )
    d_end <- ifelse(
        selected_strand == "+",
        end(selected_peaks) - start(selected_tss),
        start(selected_tss) - end(selected_peaks)
    )
    distance[align_idx] <- ifelse(abs(d_start) <= abs(d_end), d_start, d_end)

    list(index = index, distance = distance)
}

##' Align the reported transcript with the annotated transcript
##'
##' Transcript level wrapper of [`.alignAnnotationFeature()`] (`idColumn =
##' "tx_id"`).
##'
##' @inheritParams .alignAnnotationFeature
##' @return list with the aligned `index` and `distance`
##' @noRd
##' @author G Yu
.alignTranscriptAnnotation <- function(peak.gr, features, index, distance,
                                    annotationFeatureId) {
    .alignAnnotationFeature(peak.gr, features, index, distance,
                            annotationFeatureId, idColumn = "tx_id")
}
