##' read peak file and store in data.frame or GRanges object
##'
##' 
##' @title readPeakFile
##' @param peakfile peak file
##' @param as output format, one of GRanges or data.frame
##' @param ... additional parameter
##' @return peak information, in GRanges or data.frame object
##' @import IRanges GenomicRanges
##' @export
##' @examples
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peak.gr <- readPeakFile(peakfile, as="GRanges")
##' peak.gr
##' @author G Yu
readPeakFile <- function(peakfile, as="GRanges", ...) {
    as <- match.arg(as, c("GRanges", "data.frame"))
    peak.df <- peak2DF(peakfile, ...)
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

##' @importFrom utils read.delim
peak2DF <- function(peakfile, header, ...) {
    if (missing(header)) {
        ## determine file format
        if (isBedFile(peakfile)) {
            header <- FALSE
        } else {
            header <- TRUE
        }
    }
    peak.df <- read.delim(peakfile, header=header, comment.char="#", ...)
    ## coordinate system in BED file is start at 0
    ## refer to http://asia.ensembl.org/info/website/upload/bed.html?redirect=no
    ## The chromEnd base is not included in the display of the feature.
    ## For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100,
    ## and span the bases numbered 0-99.
    ## so chromEnd, peak.df[,3], is not needed to +1
    peak.df[,2] <- peak.df[,2] + 1
    return(peak.df)
}

isBedFile <- function(peakfile) {
    ## peakfile is a peak file name
    grepl("\\.bed$", peakfile) || grepl("\\.bed.gz$", peakfile) || 
    grepl("\\Peak.gz$", peakfile) || grepl("\\.bedGraph.gz$", peakfile) || 
    grepl("\\.narrowPeak$", peakfile) || grepl("\\.broadPeak$",peakfile) ||
    grepl("\\..gappedPeak$", peakfile)
}
