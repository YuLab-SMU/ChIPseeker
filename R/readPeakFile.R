##' Read peak file and convert to GRanges or data.frame
##'
##' This function reads peak files in various formats (BED, narrowPeak, broadPeak,
##' etc.) and converts them to either a GRanges object or data.frame. The function
##' automatically handles BED file coordinate conversion (0-based to 1-based).
##'
##' @description
##' The function reads tab-delimited peak files and converts them to the specified
##' output format. It automatically detects BED format files based on file extension
##' and handles the coordinate system conversion (BED files use 0-based start
##' coordinates, which are converted to 1-based for R compatibility). Additional
##' columns beyond chromosome, start, and end positions are preserved in the output.
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Detects file format based on extension (BED files: .bed, .bed.gz,
##'         .narrowPeak, .broadPeak, .gappedPeak, .bedGraph.gz, or files ending
##'         with Peak.gz)
##'   \item Reads the file using \code{utils::read.delim()} with appropriate
##'         header settings (BED files have no header, other formats may have headers)
##'   \item Converts BED coordinates from 0-based to 1-based by adding 1 to the
##'         start position (BED end position is already exclusive, so no adjustment
##'         needed)
##'   \item Converts to the requested output format (GRanges or data.frame)
##'   \item Preserves all additional columns beyond the first three (chromosome,
##'         start, end) as metadata columns in GRanges or as additional columns in
##'         data.frame
##' }
##'
##' Supported file formats include standard BED files, narrowPeak, broadPeak,
##' gappedPeak, bedGraph, and other tab-delimited formats with chromosome, start,
##' and end positions in the first three columns.
##'
##' @param peakfile character, path to the peak file. Can be a local file path or
##'   a URL. Supported formats include BED (.bed, .bed.gz), narrowPeak (.narrowPeak),
##'   broadPeak (.broadPeak), gappedPeak (.gappedPeak), bedGraph (.bedGraph.gz),
##'   or any tab-delimited file with chromosome, start, and end positions
##' @param as character, output format. One of "GRanges" (default) or "data.frame".
##'   If "GRanges", returns a GRanges object with additional columns as metadata.
##'   If "data.frame", returns a data.frame
##' @param ... additional parameters passed to \code{utils::read.delim()} for
##'   reading the file. Useful parameters include \code{sep}, \code{quote},
##'   \code{stringsAsFactors}, etc. Note: \code{header} is automatically
##'   determined based on file format, but can be explicitly specified to override
##' @return Depending on the \code{as} parameter:
##'   \itemize{
##'     \item If \code{as="GRanges"}: A GRanges object with:
##'       \itemize{
##'         \item \code{seqnames}: Chromosome names (from column 1)
##'         \item \code{ranges}: IRanges with start and end positions (from columns 2 and 3)
##'         \item Additional metadata columns: All columns beyond the first three
##'               are preserved as metadata columns with their original names
##'       }
##'     \item If \code{as="data.frame"}: A data.frame with all columns from the
##'       original file, with BED coordinates converted to 1-based if applicable
##'   }
##' @import IRanges GenomicRanges
##' @export
##' @examples
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peak.gr <- readPeakFile(peakfile, as="GRanges")
##' peak.gr
##'
##' ## Read a BED file
##' ## peak.gr <- readPeakFile("peaks.bed", as="GRanges")
##'
##' ## Read with custom parameters
##' ## peak.df <- readPeakFile("peaks.txt", as="data.frame", sep="\t")
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
##' @noRd
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
    grepl("\\.gappedPeak$", peakfile)
}
