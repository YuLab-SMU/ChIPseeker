##' read peak file and store in data.frame or GRanges object
##'
##' 
##' @title readPeakFile
##' @param peakfile peak file
##' @param as output format, one of GRanges or data.frame
##' @return peak information, in GRanges or data.frame object
##' @importFrom IRanges IRanges
##' @importFrom GenomicRanges GRanges
##' @export
##' @examples
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peak.gr <- readPeakFile(peakfile, as="GRanges")
##' peak.gr
##' @author G Yu
readPeakFile <- function(peakfile, as="GRanges") {
    as <- match.arg(as, c("GRanges", "data.frame"))
    peak.df <- peak2DF(peakfile)
    if (as == "data.frame")
        return(peak.df)
    peak.gr <- peakDF2GRanges(peak.df)
    return(peak.gr)
}

peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}

peak2DF <- function(peakfile) {
    ## determine file format
    if (isBedFile(peakfile)) {
        peak.df <- read.delim(peakfile, header=FALSE, comment.char="#")        
    } else {
        peak.df <- read.delim(peakfile, header=TRUE, comment.char="#")
    }
    return(peak.df)
}

isBedFile <- function(peakfile) {
    ## peakfile is a peak file name
    grepl("\\.bed$", peakfile) || grepl("\\.bed.gz$", peakfile) ||
        grepl("\\Peak.gz$", peakfile) || grepl("\\.bedGraph.gz$", peakfile)
}
