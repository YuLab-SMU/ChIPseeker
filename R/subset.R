##' Subset csAnno object
##'
##' This function subsets a \code{csAnno} object by filtering peaks based on
##' logical conditions, and automatically updates the associated annotation
##' statistics and detail annotations.
##'
##' @description
##' The function allows subsetting of \code{csAnno} objects using the same syntax
##' as subsetting GRanges objects. It filters the peaks in the \code{anno} slot
##' and automatically updates:
##' \itemize{
##'   \item The \code{detailGenomicAnnotation} data frame to match the filtered peaks
##'   \item The \code{annoStat} annotation statistics
##'   \item The \code{peakNum} total peak count
##' }
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Creates an index based on chromosome, start, and end positions of
##'         all peaks
##'   \item Subsets the GRanges object in the \code{anno} slot using the provided
##'         conditions
##'   \item Filters the \code{detailGenomicAnnotation} data frame to keep only
##'         rows corresponding to the remaining peaks
##'   \item Recalculates annotation statistics using \code{getGenomicAnnoStat()}
##'   \item Updates the peak count
##' }
##'
##' The \code{tssRegion}, \code{level}, and \code{hasGenomicAnnotation} slots
##' remain unchanged.
##'
##' @param x \code{csAnno} object to subset
##' @param ... logical expressions indicating which peaks to keep. Can use any
##'   column names from the GRanges metadata (e.g., \code{annotation},
##'   \code{distanceToTSS}, \code{geneId}, etc.). See \code{\link[base]{subset}}
##'   for details
##' @return A \code{csAnno} object containing only the filtered peaks, with
##'   updated annotation statistics and detail annotations
##' @importFrom S4Vectors subset
##' @importFrom BiocGenerics start
##' @importFrom BiocGenerics end
##' @method subset csAnno
##' @export
##' @examples
##' \dontrun{
##' ## Subset peaks by annotation type
##' peakAnno_promoter <- subset(peakAnno, annotation == "Promoter")
##'
##' ## Subset peaks by distance to TSS
##' peakAnno_near <- subset(peakAnno, abs(distanceToTSS) < 5000)
##'
##' ## Subset by multiple conditions
##' peakAnno_filtered <- subset(peakAnno,
##'                              annotation == "Promoter" &
##'                              abs(distanceToTSS) < 3000)
##' }
##' @author G Yu
subset.csAnno <- function(x, ... ){

  index <- paste(seqnames(x@anno),start(x@anno),end(x@anno), sep = "_")
  # subset the GRanges
  x@anno <- subset(x@anno, ...)
  index2 <- paste(seqnames(x@anno),start(x@anno),end(x@anno), sep = "_")

  # the tssRgion, level, hsaGenomicAnnotation keep unchanged

  # change the detailGenomicAnnotation
  x@detailGenomicAnnotation <- x@detailGenomicAnnotation[index %in% index2,]

  # change the annotation stat
  x@annoStat <- getGenomicAnnoStat(x@anno)

  # change peak number
  x@peakNum <-  length(x@anno)

  return(x)

}
