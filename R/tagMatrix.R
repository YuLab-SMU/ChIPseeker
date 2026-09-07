##' Prepare promoter regions from TxDb
##'
##' This function extracts promoter regions (TSS regions) from a TxDb or EnsDb
##' annotation object. It is a convenience wrapper around \code{getBioRegion()}
##' for extracting promoter regions.
##'
##' @description
##' The function extracts promoter regions defined as regions around transcription
##' start sites (TSS). For each gene or transcript, it creates a region extending
##' upstream and downstream from the TSS. The function calls \code{getBioRegion()}
##' with \code{type="start_site"}.
##'
##' @param TxDb TxDb or EnsDb annotation object. If NULL, will attempt to
##'   determine from the data
##' @param upstream numeric, distance in base pairs to extend upstream from TSS.
##'   Default is 1000
##' @param downstream numeric, distance in base pairs to extend downstream from
##'   TSS. Default is 1000
##' @param by character, one of "gene" or "transcript". Determines whether to
##'   extract promoters at gene level or transcript level. Default is "gene"
##' @return GRanges object containing promoter regions with attributes:
##'   \itemize{
##'     \item \code{type}: "start_site"
##'     \item \code{by}: the value of the \code{by} parameter
##'     \item \code{label}: generated label for the regions
##'     \item \code{upstream}: upstream distance
##'     \item \code{downstream}: downstream distance
##'   }
##' @seealso \code{\link{getBioRegion}} for more flexible region extraction
##' @export
getPromoters <- function(TxDb=NULL,
                         upstream=1000,
                         downstream=1000,
                         by = "gene") {

  getBioRegion(TxDb = TxDb,
               upstream = upstream,
               downstream = downstream,
               by = by,
               type = "start_site")
}


##' Prepare a biological region of selected feature from TxDb
##'
##' This function extracts biological regions (promoters, gene bodies, UTRs, etc.)
##' from a TxDb or EnsDb annotation object. It provides a unified interface for
##' extracting different types of genomic regions.
##'
##' @description
##' The function extracts genomic regions based on feature type and region type.
##' It supports three region types:
##' \itemize{
##'   \item \code{start_site}: Region around the transcription start site (TSS).
##'         For a transcript at chr1:1000-1400, with upstream=100 and downstream=100,
##'         this extracts chr1:900-1100 (TSS +/- 100bp)
##'   \item \code{end_site}: Region around the transcription termination site (TTS).
##'         For a transcript at chr1:1000-1400, with upstream=100 and downstream=100,
##'         this extracts chr1:1300-1500 (TTS +/- 100bp)
##'   \item \code{body}: The full feature body. For a transcript at chr1:1000-1400,
##'         this extracts chr1:1000-1400 (the entire transcript)
##' }
##'
##' @details
##' The function supports multiple feature types:
##' \itemize{
##'   \item \code{gene}: Gene-level regions
##'   \item \code{transcript}: Transcript-level regions
##'   \item \code{exon}: Exon regions
##'   \item \code{intron}: Intron regions
##'   \item \code{3UTR}, \code{5UTR}, \code{UTR}: UTR regions
##' }
##'
##' For \code{start_site} and \code{end_site} types, the function creates regions
##' extending upstream and downstream from the site. For \code{body} type, it
##' returns the full feature body. The function handles strand information
##' correctly (upstream/downstream are relative to the feature's strand).
##'
##' @param TxDb TxDb or EnsDb annotation object. If NULL, will attempt to
##'   determine from the data
##' @param upstream numeric, distance in base pairs to extend upstream from start
##'   site or end site. For \code{body} type, this can extend the 5' flank.
##'   Default is 1000
##' @param downstream numeric, distance in base pairs to extend downstream from
##'   start site or end site. For \code{body} type, this can extend the 3' flank.
##'   Default is 1000
##' @param by character, one of 'gene', 'transcript', 'exon', 'intron', '3UTR',
##'   '5UTR', 'UTR'. Determines which feature type to extract. Default is "gene"
##' @param type character, one of "start_site", "end_site", "body". Determines
##'   which part of the feature to extract. Default is "start_site"
##' @return GRanges object containing the extracted regions with attributes:
##'   \itemize{
##'     \item \code{type}: The region type ("start_site", "end_site", or "body")
##'     \item \code{by}: The feature type used
##'     \item \code{label}: Generated label for the regions
##'     \item \code{upstream}: Upstream distance (for start_site/end_site)
##'     \item \code{downstream}: Downstream distance (for start_site/end_site)
##'   }
##' @seealso \code{\link{getPromoters}} for a convenience function to get promoters,
##'   \code{\link{makeBioRegionFromGranges}} for creating regions from custom GRanges
##' @import BiocGenerics IRanges GenomicRanges
##' @importFrom yulab.utils get_cache_item
##' @author Guangchuang Yu, Ming L
##' @export
##' @import BiocGenerics IRanges GenomicRanges
##' @importFrom yulab.utils get_cache_item
##' @author Guangchuang Yu, Ming L
##' @export
getBioRegion <- function(TxDb=NULL,
                         upstream=1000,
                         downstream=1000,
                         by="gene",
                         type="start_site"){

  by <- match.arg(by, c('gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR','UTR'))
  type <- match.arg(type, c("start_site", "end_site", "body"))

  TxDb <- loadTxDb(TxDb)
  .ChIPseekerEnv(TxDb, item = ChIPseekerCache)
  # ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

  label <- make_label(type = type, by = by)


  if(by == 'gene' || by == 'transcript'){
    regions <- getGene(TxDb, by)
  }

  if (by == "exon") {
    # exonList <- get_exonList(ChIPseekerEnv)
    exonList <- get_exonList(item = ChIPseekerCache)
    regions <-  unlist(exonList)
  }

  if (by == "intron") {
    # intronList <- get_intronList(ChIPseekerEnv)
    intronList <- get_intronList(item = ChIPseekerCache)
    regions <- unlist(intronList)
  }

  if (by == "3UTR") {
    threeUTRList <- threeUTRsByTranscript(TxDb)
    regions <- unlist(threeUTRList)
  }

  if (by == "5UTR") {
    fiveUTRList <- fiveUTRsByTranscript(TxDb)
    regions <- unlist(fiveUTRList)
  }

  if (by == 'UTR'){
    three_URT <- threeUTRsByTranscript(TxDb)
    three_UTR_regions <- unlist(three_URT)
    five_UTR <- fiveUTRsByTranscript(TxDb)
    five_UTR_regions <- unlist(five_UTR)
    regions <- c(three_UTR_regions,five_UTR_regions)
  }

  if(type == "start_site"){
    coordinate<- ifelse(strand(regions) == "+", start(regions), end(regions))
  }else if(type == "end_site"){
    coordinate<- ifelse(strand(regions) == "+", end(regions), start(regions))
  }else{
    ## assign attribute
    attr(regions, 'type') = type
    attr(regions, 'by') = by
    attr(regions, 'label') = label

    return(regions)
  }

  ## issue and code obtained from Chen Ting(NIH/NCI)
  start_site <- ifelse(strand(regions) == "+",coordinate-upstream, coordinate-downstream)
  end_site <- ifelse(strand(regions) == "+", coordinate+downstream, coordinate+upstream)

  bioRegion <- GRanges(seqnames=seqnames(regions),
                       ranges=IRanges(start_site, end_site),
                       strand=strand(regions))
  bioRegion <- unique(bioRegion)

  ## assign attribute
  attr(bioRegion, 'type') = type
  attr(bioRegion, 'by') = by

  ## different region have different label to be added to the figures
  ## so we attach label to the Granges object
  attr(bioRegion, 'label') = label

  attr(bioRegion, 'upstream') = upstream
  attr(bioRegion, 'downstream') = downstream

  return(bioRegion)
}

##' Make biological regions from user-provided GRanges object
##'
##' This function creates biological regions from a user-provided GRanges object,
##' allowing extraction of regions around start sites, end sites, or full bodies
##' for custom genomic features not available in TxDb objects (e.g., enhancers,
##' insulators).
##'
##' @description
##' The function extracts regions from a custom GRanges object, similar to
##' \code{getBioRegion()} but working with user-defined features instead of TxDb
##' annotations. This is useful for features like enhancers, insulators, or
##' other regulatory elements that are not part of standard gene annotations.
##'
##' @details
##' The function supports three region types:
##' \itemize{
##'   \item \code{start_site}: Region around the start of each range. For an
##'         enhancer at chr1:1000-1400 with upstream=100 and downstream=100,
##'         extracts chr1:900-1100 (start +/- 100bp, strand-aware)
##'   \item \code{end_site}: Region around the end of each range. For an enhancer
##'         at chr1:1000-1400 with upstream=100 and downstream=100, extracts
##'         chr1:1300-1500 (end +/- 100bp, strand-aware)
##'   \item \code{body}: The full range body. For an enhancer at chr1:1000-1400,
##'         extracts chr1:1000-1400. For this type, \code{upstream} and
##'         \code{downstream} can be NULL
##' }
##'
##' The function handles strand information correctly: for negative strand features,
##' "upstream" and "downstream" are reversed relative to the genomic coordinates.
##'
##' @param gr GRanges object containing regions of interest (e.g., enhancers,
##'   insulators, custom regulatory elements)
##' @param by character, user-specified label for the feature type (e.g., "gene",
##'   "insulator", "enhancer"). Used for generating region labels. Cannot be omitted
##' @param type character, one of "start_site", "end_site", "body". Determines
##'   which part of each range to extract. Cannot be omitted
##' @param upstream numeric or NULL, distance in base pairs to extend upstream
##'   from start site or end site. Can be NULL if \code{type == 'body'}.
##'   Default is 1000
##' @param downstream numeric or NULL, distance in base pairs to extend downstream
##'   from start site or end site. Can be NULL if \code{type == 'body'}.
##'   Default is 1000
##' @return GRanges object containing the extracted regions with attributes:
##'   \itemize{
##'     \item \code{type}: The region type
##'     \item \code{by}: The feature type label
##'     \item \code{label}: Generated label for the regions
##'     \item \code{upstream}: Upstream distance (for start_site/end_site)
##'     \item \code{downstream}: Downstream distance (for start_site/end_site)
##'   }
##' @seealso \code{\link{getBioRegion}} for extracting regions from TxDb objects
##' @import BiocGenerics IRanges GenomicRanges
##' @export
##' @import BiocGenerics IRanges GenomicRanges
##' @export
makeBioRegionFromGranges <- function(gr,
                                     by,
                                     type,
                                     upstream=1000,
                                     downstream=1000){

  if (!is(gr, "GRanges")) {
    stop("windows should be a GRanges object...")
  }

  type <- match.arg(type, c("start_site", "end_site", "body"))

  label <- make_label(type = type, by = by)
  regions <- gr

  if(type == "start_site"){
    coordinate<- ifelse(strand(regions) == "+", start(regions), end(regions))
  }else if(type == "end_site"){
    coordinate<- ifelse(strand(regions) == "+", end(regions), start(regions))
  }else{
    ## assign attribute
    attr(regions, 'type') = type
    attr(regions, 'by') = by
    attr(regions, 'label') = label

    return(regions)
  }

  ## issue and code obtained from Chen Ting(NIH/NCI)
  start_site <- ifelse(strand(regions) == "+",coordinate-upstream, coordinate-downstream)
  end_site <- ifelse(strand(regions) == "+", coordinate+downstream, coordinate+upstream)

  bioRegion <- GRanges(seqnames=seqnames(regions),
                       ranges=IRanges(start_site, end_site),
                       strand=strand(regions))
  bioRegion <- unique(bioRegion)

  ## assign attribute
  attr(bioRegion, 'type') = type
  attr(bioRegion, 'by') = by
  attr(bioRegion, 'label') = label
  attr(bioRegion, 'upstream') = upstream
  attr(bioRegion, 'downstream') = downstream

  return(bioRegion)

}


##' Calculate tag matrix for peak coverage visualization
##'
##' This function calculates a tag matrix representing peak coverage across genomic
##' regions (e.g., promoters, gene bodies). The matrix can be used for visualization
##' with functions like \code{plotPeakProf()} or \code{plotAvgProf()}.
##'
##' @description
##' The function computes peak coverage across a set of genomic regions (windows).
##' It supports two methods:
##' \itemize{
##'   \item \strong{Direct method} (when \code{nbin=NULL}): For regions of equal
##'         size (e.g., promoters), calculates coverage at each base pair position
##'   \item \strong{Binning method} (when \code{nbin} is specified): For regions
##'         of variable size (e.g., gene bodies), divides each region into a fixed
##'         number of bins and calculates average coverage per bin
##' }
##'
##' @details
##' The function can work in two modes:
##' \enumerate{
##'   \item \strong{With pre-made windows}: Provide a \code{windows} GRanges object
##'         created by \code{getPromoters()}, \code{getBioRegion()}, or
##'         \code{makeBioRegionFromGranges()}
##'   \item \strong{Without windows}: Provide \code{TxDb} and region parameters
##'         (\code{type}, \code{by}, \code{upstream}, \code{downstream}), and the
##'         function will create windows automatically
##' }
##'
##' For \code{upstream} and \code{downstream} parameters:
##' \itemize{
##'   \item If \code{windows} is provided: For \code{type='body'}, these extend
##'         the flanks; for \code{type='start_site'/'end_site'}, they are ignored
##'         (taken from window attributes)
##'   \item If \code{windows} is missing: These define the region size for
##'         \code{type='start_site'/'end_site'}, or extend flanks for \code{type='body'}
##' }
##'
##' The \code{weightCol} parameter allows weighting peaks by a metadata column
##' (e.g., peak score or signal value), useful for signal-weighted coverage.
##'
##' @param peak peak file (BED format) or GRanges object containing ChIP-seq peaks
##' @param upstream numeric or \code{rel()} object, distance to extend upstream.
##'   For \code{type='body'}, can extend 5' flank. For \code{type='start_site'/'end_site'},
##'   only used when \code{windows} is missing. Can be NULL for body regions
##' @param downstream numeric or \code{rel()} object, distance to extend downstream.
##'   For \code{type='body'}, can extend 3' flank. For \code{type='start_site'/'end_site'},
##'   only used when \code{windows} is missing. Can be NULL for body regions
##' @param windows GRanges object containing regions of interest. Should be created
##'   by \code{getPromoters()}, \code{getBioRegion()}, or \code{makeBioRegionFromGranges()}.
##'   If missing, will be created from \code{TxDb} and other parameters
##' @param type character, one of "start_site", "end_site", "body". Required if
##'   \code{windows} is missing
##' @param by character, one of 'gene', 'transcript', 'exon', 'intron', '3UTR',
##'   '5UTR', or user-specified label. Required if \code{windows} is missing
##' @param TxDb TxDb/EnsDb object or GRanges object. If TxDb, used to create windows.
##'   If GRanges, treated as custom regions and passed to \code{makeBioRegionFromGranges()}.
##'   Required if \code{windows} is missing
##' @param weightCol character, name of a metadata column in the peak GRanges to
##'   use as weights for coverage calculation. If NULL, all peaks have equal weight.
##'   Default is NULL
##' @param nbin integer, number of bins for binning method. Required for
##'   \code{type='body'}, optional for other types. If NULL, uses direct method
##'   (requires equal-width windows). Default is NULL
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param ignore_strand logical, whether to ignore strand information when
##'   calculating coverage. If FALSE, negative strand regions are reverse-complemented.
##'   Default is FALSE
##' @return A matrix (tagMatrix) with:
##'   \itemize{
##'     \item Rows: One per region (window)
##'     \item Columns: One per position/bin
##'     \item Values: Peak coverage/signal at each position
##'   }
##'   The matrix has attributes:
##'   \itemize{
##'     \item \code{upstream}: Upstream distance
##'     \item \code{downstream}: Downstream distance
##'     \item \code{type}: Region type
##'     \item \code{label}: Region label
##'     \item \code{is.binning}: Whether binning method was used
##'   }
##' @seealso \code{\link{getPromoters}}, \code{\link{getBioRegion}},
##'   \code{\link{makeBioRegionFromGranges}} for creating windows,
##'   \code{\link{plotPeakProf}}, \code{\link{plotAvgProf}} for visualization
##' @importFrom ggplot2 rel
##' @export
getTagMatrix <- function(peak,
                         upstream,
                         downstream,
                         windows,
                         type,
                         by,
                         TxDb=NULL,
                         weightCol = NULL,
                         nbin = NULL,
                         verbose = TRUE,
                         ignore_strand= FALSE){

  is_GRanges_of_TxDb <- FALSE
  if (is(TxDb, "GRanges")) {
    is_GRanges_of_TxDb <- TRUE
    message("#\n#.. 'TxDb' is a self-defined 'GRanges' object...\n#")
  }

  if(missingArg(windows)){

    if(is_GRanges_of_TxDb){

      ## make windows from self-made granges object
      windows <- makeBioRegionFromGranges(gr=TxDb,
                                          by=by,
                                          type=type,
                                          upstream=upstream,
                                          downstream=downstream)

    }else{

      ## make windows from txdb object
      windows <- getBioRegion(TxDb=TxDb,
                              upstream=upstream,
                              downstream=downstream,
                              by=by,
                              type=type)


    }

  }else{

    if (!is(windows, "GRanges")) {
      stop("windows should be a GRanges object...")
    }

    if(is.null(attr(windows,'type'))){
      stop("windows should be made from getPromoters()/getBioRegion()/makeBioRegionFromGranges()")
    }

    type <- attr(windows, 'type')
    by <- attr(windows, 'by')

  }

  # check the upstream and downstream parameter
  if(type == "body"){
    if(missingArg(upstream)){
      upstream <- NULL
    }

    if(missingArg(downstream)){
      downstream <- NULL
    }

  }else{
    upstream <- attr(windows, 'upstream')
    downstream <- attr(windows, 'downstream')
  }

  ## check upstream and downstream parameter
  check_upstream_and_downstream(upstream = upstream, downstream = downstream)

  if(type != 'body'){
    if(inherits(upstream, 'rel') || is.null(upstream)){
      stop("upstream and downstream for site region should be actual number...")
    }
  }

  ## check nbin parameters
  if(!is.null(nbin) && !is.numeric(nbin)){
    stop('nbin should be NULL or numeric...')
  }

  if(type == 'body' && is.null(nbin)){
    stop('plotting body region should set the nbin parameter...')
  }

  ## check nbin parameter
  if(!is.null(nbin)){
    cat(">> binning method is used...",
        format(Sys.time(), "%Y-%m-%d %X"), "\n",sep = "")

    is.binning <- TRUE
  }else{

    is.binning <- FALSE
  }

  if (verbose) {
    cat(">> preparing ",type," regions"," by ",by,"... ",
        format(Sys.time(), "%Y-%m-%d %X"), "\n",sep = "")
  }


  if(is.binning){

    if (verbose) {
      cat(">> preparing tag matrix by binning... ",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }

    tagMatrix <- getTagMatrix.binning.internal(peak = peak,
                                               weightCol = weightCol,
                                               windows = windows,
                                               nbin = nbin,
                                               upstream = upstream,
                                               downstream = downstream,
                                               ignore_strand = ignore_strand)
  }else{

    if (verbose) {
      cat(">> preparing tag matrix... ",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }

    tagMatrix <- getTagMatrix.internal(peak=peak,
                                       weightCol=weightCol,
                                       windows=windows,
                                       ignore_strand=ignore_strand)
  }

  ## assign attribute
  attr(tagMatrix, 'upstream') = upstream
  attr(tagMatrix, 'downstream') = downstream
  attr(tagMatrix, 'type') = attr(windows, 'type')
  attr(tagMatrix, 'label') = attr(windows, 'label')
  attr(tagMatrix, "is.binning") <- is.binning

  return(tagMatrix)
}


##' Calculate tag matrix for equal-width regions (internal function)
##'
##' This is an internal function that calculates peak coverage across regions of
##' equal width. It is called by \code{getTagMatrix()} when \code{nbin=NULL}
##' (direct method, not binning).
##'
##' @description
##' The function calculates coverage at each base pair position across all windows.
##' It requires that all windows have the same width. For negative strand regions,
##' it reverse-complements the coverage to align all regions in the same orientation
##' (5' to 3').
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Loads peaks and calculates coverage (optionally weighted)
##'   \item Filters windows to those within chromosome boundaries
##'   \item Extracts coverage views for each window
##'   \item Creates a matrix with one row per window and one column per position
##'   \item Reverse-complements negative strand regions (if \code{ignore_strand=FALSE})
##'   \item Removes windows with zero coverage
##' }
##'
##' @param peak peak file (BED format) or GRanges object containing ChIP-seq peaks
##' @param weightCol character, name of a metadata column to use as weights.
##'   If NULL, all peaks have equal weight. Default is NULL
##' @param windows GRanges object containing regions. All regions must have the
##'   same width (e.g., promoter regions of fixed size)
##' @param ignore_strand logical, whether to ignore strand information. If FALSE,
##'   negative strand regions are reverse-complemented. Default is FALSE
##' @return A matrix with:
##'   \itemize{
##'     \item Rows: One per window (regions with non-zero coverage)
##'     \item Columns: One per base pair position
##'     \item Values: Coverage/signal at each position
##'   }
##'   Row names correspond to the original window indices
##' @import BiocGenerics S4Vectors IRanges GenomeInfoDb GenomicRanges
##' @author G Yu
##' @noRd
getTagMatrix.internal <- function(peak,
                                  weightCol=NULL,
                                  windows,
                                  ignore_strand= FALSE) {
  peak.gr <- loadPeak(peak)

  if (! is(windows, "GRanges")) {
    stop("windows should be a GRanges object...")
  }
  if (length(unique(width(windows))) != 1) {
    stop("width of windows should be equal...")
  }

  ## if (!exists("ChIPseekerEnv", envir = .GlobalEnv)) {
  ##     assign("ChIPseekerEnv", new.env(), .GlobalEnv)
  ## }
  ## ChIPseekerEnv <- get("ChIPseekerEnv", envir = .GlobalEnv)

  ## if (exists("peak", envir=ChIPseekerEnv, inherits=FALSE) &&
  ##     exists("promoters", envir=ChIPseekerEnv, inherits=FALSE) &&
  ##     exists("weightCol", envir=ChIPseekerEnv, inherits=FALSE) &&
  ##     exists("tagMatrix", envir=ChIPseekerEnv, inherits=FALSE) ) {

  ##     pp <- get("peak", envir=ChIPseekerEnv)
  ##     promoters <- get("promoters", envir=ChIPseekerEnv)
  ##     w <- get("weightCol", envir=ChIPseekerEnv)

  ##     if (all(pp == peak)) {
  ##         if (all(windows == promoters)) {
  ##             if ( (is.null(w) && is.null(weightCol)) ||
  ##                 (!is.null(w) && !is.null(weightCol) && w == weightCol)) {
  ##                 tagMatrix <- get("tagMatrix", envir=ChIPseekerEnv)
  ##                 return(tagMatrix)
  ##             } else {
  ##                 assign("weightCol", weightCol, envir=ChIPseekerEnv)
  ##             }
  ##         } else {
  ##             assign("promoters", windows)
  ##             ## make sure it is not conflict with getPromoters
  ##             if ( exists("upstream", envir=ChIPseekerEnv, inherits=FALSE))
  ##                 rm("upstream", envir=ChIPseekerEnv)
  ##         }
  ##     } else {
  ##         assign("peak", peak, envir=ChIPseekerEnv)
  ##     }

  ## }

  ## if ( !exists("peak", envir=ChIPseekerEnv, inherits=FALSE)) {
  ##     assign("peak", peak, envir=ChIPseekerEnv)
  ## }

  ## if ( !exists("promoters", envir=ChIPseekerEnv, inherits=FALSE)) {
  ##     assign("promoters", windows, envir=ChIPseekerEnv)
  ## }

  ## if (!exists("weightCol", envir=ChIPseekerEnv, inherits=FALSE)) {
  ##     assign("weightCol", weightCol, envir=ChIPseekerEnv)
  ## }
  if (is.null(weightCol)) {
    peak.cov <- coverage(peak.gr)
  } else {
    weight <- mcols(peak.gr)[[weightCol]]
    peak.cov <- coverage(peak.gr, weight=weight)
  }
  cov.len <- elementNROWS(peak.cov)
  cov.width <- GRanges(seqnames=names(cov.len),
                       IRanges(start=rep(1, length(cov.len)),
                               end=cov.len))
  windows <- subsetByOverlaps(windows, cov.width,
                              type="within", ignore.strand=FALSE)

  chr.idx <- intersect(names(peak.cov),
                       unique(as.character(seqnames(windows))))

  peakView <- Views(peak.cov[chr.idx], as(windows, "IntegerRangesList")[chr.idx])
  tagMatrixList <- lapply(peakView, function(x) t(viewApply(x, as.vector)))
  tagMatrix <- do.call("rbind", tagMatrixList)

  ## get the index of windows, that are reorganized by as(windows, "IntegerRangesList")
  idx.list <- split(1:length(windows),  as.factor(seqnames(windows)))
  idx <- unlist(idx.list[chr.idx], use.names=FALSE)
  
  rownames(tagMatrix) <- idx
  tagMatrix <- tagMatrix[order(idx),]

  ## minus strand
  if (!ignore_strand) {
    minus.idx <- which(as.character(strand(windows)) == "-")
    tagMatrix[minus.idx,] <- tagMatrix[minus.idx, ncol(tagMatrix):1]
  }

  tagMatrix <- tagMatrix[rowSums(tagMatrix)!=0,]
  ## assign("tagMatrix", tagMatrix, envir=ChIPseekerEnv)
  return(tagMatrix)
}


##' Calculate tag matrix using binning method (internal function)
##'
##' This is an internal function that calculates peak coverage using a binning
##' approach, inspired by deeptools computeMatrix. It is called by \code{getTagMatrix()}
##' when \code{nbin} is specified.
##'
##' @description
##' The function divides each region into a fixed number of bins and calculates
##' average coverage per bin. This allows handling regions of variable size (e.g.,
##' gene bodies of different lengths) by normalizing them to the same number
##' of bins. The approach is similar to deeptools computeMatrix
##' (\url{https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html}).
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Loads peaks and calculates coverage (optionally weighted)
##'   \item Filters windows to those within chromosome boundaries
##'   \item Optionally extends windows by \code{upstream} and \code{downstream}
##'         (can be relative values using \code{rel()})
##'   \item For each window, divides it into \code{nbin} equal-sized bins
##'   \item Calculates average coverage per bin
##'   \item Creates a matrix with one row per window and one column per bin
##'   \item Reverse-complements negative strand regions (if \code{ignore_strand=FALSE})
##'   \item Filters out windows shorter than the minimum required length
##' }
##'
##' The binning method is essential for visualizing coverage across variable-length
##' regions (e.g., gene bodies), as it normalizes all regions to the same number
##' of data points.
##'
##' @param peak peak file (BED format) or GRanges object containing ChIP-seq peaks
##' @param weightCol character, name of a metadata column to use as weights.
##'   If NULL, all peaks have equal weight. Default is NULL
##' @param windows GRanges object containing regions. Regions can have variable
##'   sizes (e.g., gene bodies of different lengths)
##' @param nbin integer, number of bins to divide each region into. Should not
##'   exceed the minimum region length. Default is 800
##' @param upstream numeric, \code{rel()} object, or NULL. Distance to extend
##'   upstream. For \code{type='body'}, can extend 5' flank. If \code{rel()},
##'   extends by a fraction of region width. Default is NULL
##' @param downstream numeric, \code{rel()} object, or NULL. Distance to extend
##'   downstream. For \code{type='body'}, can extend 3' flank. If \code{rel()},
##'   extends by a fraction of region width. Default is NULL
##' @param ignore_strand logical, whether to ignore strand information. If FALSE,
##'   negative strand regions are reverse-complemented. Default is FALSE
##' @return A matrix with:
##'   \itemize{
##'     \item Rows: One per window (regions meeting minimum length requirement)
##'     \item Columns: One per bin (\code{nbin} columns)
##'     \item Values: Average coverage/signal per bin
##'   }
##' @import BiocGenerics S4Vectors IRanges GenomeInfoDb GenomicRanges
##' @importFrom ggplot2 rel
getTagMatrix.binning.internal <- function(peak,
                                          weightCol = NULL,
                                          windows,
                                          nbin = 800,
                                          upstream = NULL,
                                          downstream = NULL,
                                          ignore_strand = FALSE){

  min_body_length <- filter_length <- nbin
  peak.gr <- loadPeak(peak)
  type <- attr(windows, 'type')


  if (!is(windows, "GRanges")) {
    stop("windows should be a GRanges object...")
  }

  if (is.null(weightCol)) {
    peak.cov <- coverage(peak.gr)
  } else {
    weight <- mcols(peak.gr)[[weightCol]]
    peak.cov <- coverage(peak.gr, weight=weight)
  }


  cov.len <- elementNROWS(peak.cov)
  cov.width <- GRanges(seqnames=names(cov.len),
                       IRanges(start=rep(1, length(cov.len)),
                               end=cov.len))

  windows <- subsetByOverlaps(windows,
                              cov.width,
                              type="within",
                              ignore.strand=FALSE)

  ## extend the windows by rel object
  if(inherits(upstream, 'rel')){

    windows1 <- windows

    if(!ignore_strand){

      positive_index <- which(as.character(strand(windows1)) == "+")
      negative_index <- which(as.character(strand(windows1)) == "-")
      start(windows1)[positive_index] <- suppressWarnings(start(windows1)[positive_index] - floor(width(windows)[positive_index]*as.numeric(upstream)))
      end(windows1)[positive_index] <- suppressWarnings(end(windows1)[positive_index] + floor(width(windows)[positive_index]*as.numeric(downstream)))

      start(windows1)[negative_index] <- suppressWarnings(start(windows1)[negative_index] - floor(width(windows)[negative_index]*as.numeric(downstream)))
      end(windows1)[negative_index] <- suppressWarnings(end(windows1)[negative_index] + floor(width(windows)[negative_index]*as.numeric(upstream)))

    }else{

      start(windows1) <- suppressWarnings(start(windows1) - floor(width(windows)*as.numeric(upstream)))
      end(windows1) <- suppressWarnings(end(windows1) + floor(width(windows)*as.numeric(downstream)))

    }

    windows <- windows1
    nbin <- floor(nbin*(1+as.numeric(downstream)+as.numeric(upstream)))
    min_body_length <- min_body_length*(1+as.numeric(upstream)+as.numeric(downstream))

    cat(">> preparing matrix with extension from (",attr(windows,'label')[1],"-",
        100*as.numeric(upstream),"%)~(",attr(windows,'label')[2],"+",
        100*as.numeric(downstream),"%)... ",
        format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
  }

  ## do not extend
  if(is.null(upstream)){
    if(attr(windows, 'type') == 'body'){
      cat(">> preparing matrix for ",attr(windows, 'type')," region with no flank extension... ",
          format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
    }else{
      cat(">> preparing matrix for ",attr(windows,'type')," region... ",
          format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
    }
  }

  ## extend the windows by actual number
  if(!is.null(upstream) && !inherits(upstream, 'rel') && attr(windows, 'type')== 'body'){

    windows1 <- windows

    if(!ignore_strand){

      positive_index <- which(as.character(strand(windows1)) == "+")
      negative_index <- which(as.character(strand(windows1)) == "-")

      start(windows1)[positive_index] <- suppressWarnings(start(windows1)[positive_index] - upstream)
      end(windows1)[positive_index] <- suppressWarnings(end(windows1)[positive_index] + downstream)

      start(windows1)[negative_index] <- suppressWarnings(start(windows1)[negative_index] - downstream)
      end(windows1)[negative_index] <- suppressWarnings(end(windows1)[negative_index] + upstream)

    }else{

      start(windows1) <- suppressWarnings(start(windows1) - upstream)
      end(windows1) <- suppressWarnings(end(windows1) + downstream)

    }

    windows <- windows1
    upstreamPer <- floor(upstream/1000)*0.1
    downstreamPer <- floor(downstream/1000)*0.1
    nbin <- floor(nbin*(1+upstreamPer+downstreamPer))
    min_body_length <- min_body_length+upstream+downstream

    cat(">> preparing matrix with flank extension from (",attr(windows,'label')[1],"-",
        upstream,"bp)~(",attr(windows,'label')[2],"+",downstream,"bp)... ",
        format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")
  }

  chr.idx <- intersect(names(peak.cov),
                       unique(as.character(seqnames(windows))))

  idx.list <- split(seq_len(length(windows)), as.factor(seqnames(windows)))
  idx.list <- idx.list[chr.idx]
  
  windows <- as(windows, "IntegerRangesList")[chr.idx]
  attr(windows,'type') <- type

  peakView <- Views(peak.cov[chr.idx],
                    windows)

  ## remove the gene that has no binding proteins
  for (i in 1:length(peakView)) {

    index <- viewSums(peakView[[i]])!= 0
    peakView[[i]] <- peakView[[i]][index]
    windows[[i]] <- windows[[i]][index]
    idx.list[[i]] <- idx.list[[i]][index]
  } 
  
  tagMatrixList <- lapply(peakView, function(x) viewApply(x, as.vector))

  if(!attr(windows, 'type') == 'body'){

    tagMatrixList <- lapply(tagMatrixList, function(x) t(x))

    # to remove the chromosome that do not bind protein
    index <- vapply(tagMatrixList, function(x) length(x)>0, FUN.VALUE = logical(1))
    tagMatrixList <- tagMatrixList[index]
    windows <- windows[index]
    idx.list <- idx.list[index]
    
    ## create a matrix to receive binning results
    tagMatrix <- list()

    ## this circulation is to deal with different chromosomes
    for (i in 1:length(tagMatrixList)) {

      tagMatrix[[i]] <- matrix(nrow = nrow(tagMatrixList[[i]]),ncol = nbin)

      ## this circulation is to deal with different genes
      for (j in 1:nrow(tagMatrixList[[i]])) {

        ## seq is the distance between different bins
        seq <- floor(length(tagMatrixList[[i]][j,])/nbin)

        ## cursor record the position of calculation
        cursor <- 1

        ## the third circulation is to calculate the binding strength
        ## it has two parts
        ## the first part is to for the nbin(1:nbin-1)
        ## because the seq is not derived from exact division
        ## the second part is to compensate the loss of non-exact-division

        ## this the first part for 1:(nbin-1)
        for (k in 1:(nbin-1)) {

          read <- 0

          for (z in cursor:(cursor+seq-1)) {
            read <- read + tagMatrixList[[i]][j,z]
          }

          tagMatrix[[i]][j,k] <- read/seq

          cursor <- cursor+seq
        }

        ## this the second part to to compensate the loss of non-exact-division
        read <- 0
        for (z in cursor:length(tagMatrixList[[i]][j,])) {
          read <- read+tagMatrixList[[i]][j,z]
        }

        tagMatrix[[i]][j,nbin] <- read/(length(tagMatrixList[[i]][j,])-cursor+1)
      }

      if(!ignore_strand){
        minus.idx <- which(as.character(mcols(windows[[i]])[["strand"]]) == "-")
        tagMatrix[[i]][minus.idx,] <- tagMatrix[[i]][minus.idx, ncol(tagMatrix[[i]]):1]
      }
    }

  }else{

    ## extend genebody by atual number
    if(!is.null(upstream) & !inherits(upstream, 'rel')){

      for (i in 1:length(tagMatrixList)) {
        if (length(class(tagMatrixList[[i]])) != 1) {
          sample <- tagMatrixList[[i]]
          tagMatrixList[[i]] <- lapply(seq_len(ncol(sample)), function(i) sample[,i])
        }
      }

      index <- vapply(tagMatrixList, function(x) length(x)>0, FUN.VALUE = logical(1))
      tagMatrixList <- tagMatrixList[index]
      windows <- windows[index]
      idx.list <- idx.list[index]
      
      ## count the amount before filtering
      pre_amount <- 0
      for(i in 1:length(tagMatrixList)){
        pre_amount <- pre_amount+length(tagMatrixList[[i]])
      }

      for (i in 1:length(tagMatrixList)) {

        index <- vapply(tagMatrixList[[i]], function(y) length(y)>min_body_length,FUN.VALUE = logical(1))
        tagMatrixList[[i]] <- tagMatrixList[[i]][index]
        windows[[i]] <- windows[[i]][index]
        idx.list[[i]] <- idx.list[[i]][index]
      }

      ## count the amount after filtering
      amount <- 0
      for(i in 1:length(tagMatrixList)){
        amount <- amount+length(tagMatrixList[[i]])
      }

      cat(">> ",pre_amount-amount," peaks(",100*((pre_amount-amount)/pre_amount),
          "%), having lengths smaller than ",filter_length,"bp, are filtered... ",
          format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")

      upstreamnbin <- floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer)))
      bodynbin <- floor(nbin*(1/(1+upstreamPer+downstreamPer)))
      downstreamnbin <- floor(nbin*(downstreamPer/(1+upstreamPer+downstreamPer)))

      tagMatrix <- list()

      for (i in 1:length(tagMatrixList)) {

        tagMatrix[[i]] <- matrix(nrow = length(tagMatrixList[[i]]),ncol = nbin)

        ## count the upstream
        for (j in 1:length(tagMatrixList[[i]])) {

          seq <- floor(upstream/upstreamnbin)
          cursor <- 1

          for (k in 1:(upstreamnbin-1)) {

            read <- 0

            for (z in cursor:(cursor+seq-1)) {
              read <- read + tagMatrixList[[i]][[j]][z]
            }

            tagMatrix[[i]][j,k] <- read/seq

            cursor <- cursor+seq
          }


          read <- 0
          for (z in cursor:upstream) {
            read <- read+tagMatrixList[[i]][[j]][z]
          }

          tagMatrix[[i]][j,upstreamnbin] <- read/(upstream-cursor)

        }

        ## count genebody
        for (j in 1:length(tagMatrixList[[i]])) {

          seq <- floor((length(tagMatrixList[[i]][[j]])-upstream-downstream)/bodynbin)
          cursor <- upstream+1

          for (k in (upstreamnbin+1):(upstreamnbin+bodynbin-1)) {

            read <- 0

            for (z in cursor:(cursor+seq-1)) {
              read <- read + tagMatrixList[[i]][[j]][z]
            }

            tagMatrix[[i]][j,k] <- read/seq

            cursor <- cursor+seq
          }

          read <- 0
          for (z in cursor:(length(tagMatrixList[[i]][[j]])-downstream)) {
            read <- read+tagMatrixList[[i]][[j]][z]
          }

          tagMatrix[[i]][j,bodynbin+upstreamnbin] <- read/(length(tagMatrixList[[i]][[j]])-downstream-cursor)
        }

        ## count downstream
        for (j in 1:length(tagMatrixList[[i]])) {

          seq <- floor(downstream/downstreamnbin)
          cursor <- length(tagMatrixList[[i]][[j]])-downstream+1

          for (k in (upstreamnbin+bodynbin+1):(nbin-1)) {

            read <- 0

            for (z in cursor:(cursor+seq-1)) {
              read <- read + tagMatrixList[[i]][[j]][z]
            }

            tagMatrix[[i]][j,k] <- read/seq

            cursor <- cursor+seq
          }

          read <- 0
          for (z in cursor:length(tagMatrixList[[i]][[j]])) {
            read <- read+tagMatrixList[[i]][[j]][z]
          }

          tagMatrix[[i]][j,nbin] <- read/(length(tagMatrixList[[i]][[j]])-cursor+1)
        }

        if(!ignore_strand){
          minus.idx <- which(as.character(mcols(windows[[i]])[["strand"]]) == "-")
          tagMatrix[[i]][minus.idx,] <- tagMatrix[[i]][minus.idx, ncol(tagMatrix[[i]]):1]
        }

      }

    }else{

      for (i in 1:length(tagMatrixList)) {
        if (length(class(tagMatrixList[[i]])) != 1) {
          sample <- tagMatrixList[[i]]
          tagMatrixList[[i]] <- lapply(seq_len(ncol(sample)), function(i) sample[,i])
        }
      }

      index <- vapply(tagMatrixList, function(x) length(x)>0, FUN.VALUE = logical(1))
      tagMatrixList <- tagMatrixList[index]
      windows <- windows[index]
      idx.list <- idx.list[index]
      
      ## count the amount before filtering
      pre_amount <- 0
      for(i in 1:length(tagMatrixList)){
        pre_amount <- pre_amount+length(tagMatrixList[[i]])
      }

      for (i in 1:length(tagMatrixList)) {

        index <- vapply(tagMatrixList[[i]], function(y) length(y)>min_body_length,FUN.VALUE = logical(1))
        tagMatrixList[[i]] <- tagMatrixList[[i]][index]
        windows[[i]] <- windows[[i]][index]
        idx.list[[i]] <- idx.list[[i]][index]
      }

      ## count the amount after filtering
      amount <- 0
      for(i in 1:length(tagMatrixList)){
        amount <- amount+length(tagMatrixList[[i]])
      }

      cat(">> ",pre_amount-amount," peaks(",100*((pre_amount-amount)/pre_amount),
          "%), having lengths smaller than ",filter_length,"bp, are filtered... ",
          format(Sys.time(), "%Y-%m-%d %X"),"\n",sep = "")

      tagMatrix <- list()

      for (i in 1:length(tagMatrixList)) {

        tagMatrix[[i]] <- matrix(nrow = length(tagMatrixList[[i]]),ncol = nbin)

        for (j in 1:length(tagMatrixList[[i]])) {

          seq <- floor(length(tagMatrixList[[i]][[j]])/nbin)
          cursor <- 1

          for (k in 1:(nbin-1)) {

            read <- 0

            for (z in cursor:(cursor+seq-1)) {
              read <- read + tagMatrixList[[i]][[j]][z]
            }

            tagMatrix[[i]][j,k] <- read/seq

            cursor <- cursor+seq
          }

          read <- 0
          for (z in cursor:length(tagMatrixList[[i]][[j]])) {
            read <- read+tagMatrixList[[i]][[j]][z]
          }

          tagMatrix[[i]][j,nbin] <- read/(length(tagMatrixList[[i]][[j]])-cursor+1)

        }

        if(!ignore_strand){
          minus.idx <- which(as.character(mcols(windows[[i]])[["strand"]]) == "-")
          tagMatrix[[i]][minus.idx,] <- tagMatrix[[i]][minus.idx, ncol(tagMatrix[[i]]):1]
        }
      }

    }

  }

  ## combine the results
  tagMatrix <- do.call("rbind",tagMatrix)
  idx <- unlist(idx.list, use.names=FALSE)
  tagMatrix <- tagMatrix[order(idx),,drop=FALSE]
  
  return(tagMatrix)
}


##' Nested function for getTagMatrix() to deal with multiple windows
##'
##' This is an internal function.
##' Calculate tag matrices for multiple window sets
##'
##' This function calculates tag matrices for multiple sets of windows simultaneously.
##' It is similar to \code{getTagMatrix()} but allows processing multiple feature
##' types or custom regions in a single call.
##'
##' @title getTagMatrix2
##'
##' @param peak peak file (BED format) or GRanges object containing ChIP-seq peaks
##' @param upstream numeric or \code{rel()} object, distance to extend upstream.
##'   See \code{\link{getTagMatrix}} for details
##' @param downstream numeric or \code{rel()} object, distance to extend downstream.
##'   See \code{\link{getTagMatrix}} for details
##' @param windows_name character vector, names for the window sets. Should have
##'   the same length as \code{by}. These names will be used to label the output
##'   matrices
##' @param type character, one of "start_site", "end_site", "body". Applied to
##'   all window sets
##' @param by character vector, feature types or labels. Each element can be:
##'   \itemize{
##'     \item A standard type: 'gene', 'transcript', 'exon', 'intron', '3UTR',
##'           '5UTR', 'UTR'
##'     \item A user-specified label for custom regions (must have corresponding
##'           GRanges in \code{TxDb})
##'   }
##' @param TxDb named list, where:
##'   \itemize{
##'     \item Names correspond to elements in \code{by}
##'     \item Values are TxDb/EnsDb objects (for standard types) or GRanges
##'           objects (for custom types)
##'   }
##' @param weightCol character, name of a metadata column to use as weights.
##'   If NULL, all peaks have equal weight. Default is NULL
##' @param nbin integer, number of bins for binning method. Required for
##'   \code{type='body'}. Default is NULL
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param ignore_strand logical, whether to ignore strand information.
##'   Default is FALSE
##' @return A named list of tag matrices, one per window set. Each matrix has
##'   the same structure as returned by \code{getTagMatrix()}, with attributes
##'   for upstream, downstream, type, label, and is.binning
##' @importFrom ggplot2 rel
##' @noRd
getTagMatrix2 <- function(peak,
                          upstream,
                          downstream,
                          windows_name,
                          type,
                          by,
                          TxDb=NULL,
                          weightCol = NULL,
                          nbin = NULL,
                          verbose = TRUE,
                          ignore_strand= FALSE){

  names(TxDb) <- by

  windows <- lapply(as.list(by), function(x){

    if(x %in% c('gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR', 'UTR')){

      result <- getBioRegion(TxDb=TxDb[[x]],
                             upstream=upstream,
                             downstream=downstream,
                             by=x,
                             type=type)
    }else{

      result <- makeBioRegionFromGranges(gr=TxDb[[x]],
                                         by=x,
                                         type=type,
                                         upstream=upstream,
                                         downstream=downstream)

    }

    return(result)

  })

  names(windows) <- windows_name

  # check the upstream and downstream parameter for body
  if(type == "body"){
    if(missingArg(upstream)){
      upstream <- NULL
    }

    if(missingArg(downstream)){
      downstream <- NULL
    }

  }else{
    upstream <- attr(windows[[1]], 'upstream')
    downstream <- attr(windows[[1]], 'downstream')
  }

  ## check upstream and downstream parameter
  check_upstream_and_downstream(upstream = upstream, downstream = downstream)

  if(type != 'body'){
    if(inherits(upstream, 'rel') || is.null(upstream)){
      stop("upstream and downstream for site region should be actual number...")
    }
  }

  ## check nbin parameters
  if(!is.null(nbin) && !is.numeric(nbin)){
    stop('nbin should be NULL or numeric...')
  }

  if(type == 'body' && is.null(nbin)){
    stop('plotting body region should set the nbin parameter...')
  }

  ## check nbin parameter
  if(!is.null(nbin)){
    cat(">> binning method is used...",
        format(Sys.time(), "%Y-%m-%d %X"), "\n",sep = "")

    is.binning <- TRUE
  }else{

    is.binning <- FALSE
  }

  if (verbose) {
    cat(">> preparing ",type," regions"," by ",paste(by,collapse = " "),"... ",
        format(Sys.time(), "%Y-%m-%d %X"), "\n",sep = "")
  }


  if(is.binning){

    if (verbose) {
      cat(">> preparing tag matrix by binning... ",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }

    tagMatrix <- getTagMatrix2.binning.internal(peak = peak,
                                                weightCol = weightCol,
                                                windows = windows,
                                                windows_name=windows_name,
                                                nbin = nbin,
                                                upstream = upstream,
                                                downstream = downstream,
                                                ignore_strand = ignore_strand)
  }else{

    if (verbose) {
      cat(">> preparing tag matrix... ",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }

    tagMatrix <- getTagMatrix2.internal(peak=peak,
                                        weightCol=weightCol,
                                        windows=windows,
                                        windows_name=windows_name,
                                        ignore_strand=ignore_strand)
  }

  names(tagMatrix) <- windows_name

  ## assign attribute
  tagMatrix <- lapply(tagMatrix, function(x){
    attr(x, 'upstream') = upstream
    attr(x, 'downstream') = downstream
    attr(x, 'type') = attr(windows[[1]], 'type')
    attr(x, 'label') = attr(windows[[1]], 'label')
    attr(x, "is.binning") <- is.binning
    return(x)
  })

  return(tagMatrix)

}

##' Calculate tag matrices for multiple window sets (internal function, direct method)
##'
##' This is an internal function that calculates tag matrices for multiple window
##' sets using the direct method (no binning). It is called by \code{getTagMatrix2()}
##' when \code{nbin=NULL}.
##'
##' @description
##' The function processes each window set separately using \code{getTagMatrix.internal()},
##' then returns a list of matrices. All window sets must contain equal-width regions.
##'
##' @param peak peak file (BED format) or GRanges object containing ChIP-seq peaks
##' @param windows named list of GRanges objects, one per window set. Names should
##'   match \code{windows_name}
##' @param windows_name character vector, names of the window sets
##' @param weightCol character, name of a metadata column to use as weights.
##'   If NULL, all peaks have equal weight. Default is NULL
##' @param ignore_strand logical, whether to ignore strand information.
##'   Default is FALSE
##' @return A named list of tag matrices, one per window set. Each matrix has
##'   the same structure as returned by \code{getTagMatrix.internal()}
##' @noRd
getTagMatrix2.internal <- function(peak,
                                   weightCol=NULL,
                                   windows,
                                   windows_name,
                                   ignore_strand= FALSE) {

  mt_list <- lapply(windows_name, function(x){

    windows_tmp <- windows[[x]]

    mt <- getTagMatrix.internal(peak=peak,
                                weightCol=weightCol,
                                windows=windows_tmp,
                                ignore_strand=ignore_strand)

    return(mt)
  })

  return(mt_list)
}

##' Calculate tag matrices for multiple window sets using binning (internal function)
##'
##' This is an internal function that calculates tag matrices for multiple window
##' sets using the binning method. It is called by \code{getTagMatrix2()} when
##' \code{nbin} is specified.
##'
##' @description
##' The function processes each window set separately using \code{getTagMatrix.binning.internal()},
##' then returns a list of matrices. This allows handling variable-length regions
##' across multiple feature types.
##'
##' @param peak peak file (BED format) or GRanges object containing ChIP-seq peaks
##' @param upstream numeric, \code{rel()} object, or NULL. Distance to extend
##'   upstream. See \code{\link{getTagMatrix.binning.internal}} for details
##' @param downstream numeric, \code{rel()} object, or NULL. Distance to extend
##'   downstream. See \code{\link{getTagMatrix.binning.internal}} for details
##' @param windows named list of GRanges objects, one per window set. Names should
##'   match \code{windows_name}
##' @param windows_name character vector, names of the window sets
##' @param weightCol character, name of a metadata column to use as weights.
##'   If NULL, all peaks have equal weight. Default is NULL
##' @param nbin integer, number of bins to divide each region into. Default is 800
##' @param ignore_strand logical, whether to ignore strand information.
##'   Default is FALSE
##' @return A named list of tag matrices, one per window set. Each matrix has
##'   the same structure as returned by \code{getTagMatrix.binning.internal()}
##' @noRd
getTagMatrix2.binning.internal <- function(peak,
                                           weightCol = NULL,
                                           windows,
                                           windows_name,
                                           nbin = 800,
                                           upstream = NULL,
                                           downstream = NULL,
                                           ignore_strand = FALSE){

  mt_list <- lapply(windows_name, function(x){

    windows_tmp <- windows[[x]]

    mt <- getTagMatrix.binning.internal(peak = peak,
                                        weightCol = weightCol,
                                        windows = windows_tmp,
                                        nbin = nbin,
                                        upstream = upstream,
                                        downstream = downstream,
                                        ignore_strand = ignore_strand)

    return(mt)
  })

  return(mt_list)
  
}
