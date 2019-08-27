##' @importFrom S4Vectors subset
##' @importFrom S4Vectors length
##' @importFrom BiocGenerics start
##' @importFrom BiocGenerics end
##' @export
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
