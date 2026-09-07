##' Plot peak profiles around genomic features
##'
##' This function visualizes the distribution of ChIP-seq peaks around genomic
##' features (e.g., TSS, gene body, exons) by plotting average peak profiles
##' or heatmaps.
##'
##' @description
##' The function is a flexible entry point for plotting peak profiles. It can
##' work with either pre-computed tag matrices or raw peak files/GRanges objects.
##' When provided with peaks, it automatically generates tag matrices and plots
##' them. The function routes to different plotting functions based on the input
##' type and parameters.
##'
##' @details
##' The function behavior depends on the input:
##' \itemize{
##'   \item If \code{tagMatrix} is provided: plots directly from the tag matrix
##'         (supports both normal and binning modes)
##'   \item If \code{tagMatrix} is NULL: generates tag matrix from peaks and routes to:
##'     \itemize{
##'       \item \code{plotPeakProf2()} if \code{by} has length 1
##'       \item \code{plotPeakProf_MultiWindows()} if \code{by} has length > 1
##'     }
##' }
##'
##' The function supports:
##' \itemize{
##'   \item Single or multiple peak sets (list of peaks)
##'   \item Confidence intervals for average profiles
##'   \item Binning mode for large datasets
##'   \item Faceting options for comparing multiple datasets
##' }
##'
##' \code{TxDb} parameter can accept txdb object.
##' But many regions can not be obtained by txdb object. In this case,
##' Users can provide self-made granges served the same role
##' as txdb object and pass to \code{TxDb} object.
##'
##' \code{by} the features of interest.
##'
##' (1) if users use \code{txdb}, \code{by} can be one of 'gene', 'transcript', 'exon',
##' 'intron' , '3UTR' , '5UTR', 'UTR'. These features can be obtained by functions from txdb object.
##'
##' (2) if users use self-made granges object, \code{by} can be everything. Because this \code{by}
##' will not pass to functions to get features, which is different from the case of using
##' txdb object. This \code{by} is only used to made labels showed in picture.
##'
##' \code{type} means the property of the region. one of the "start site",
##' "end site" and "body".
##'
##' \code{upstream} and \code{downstream} parameter have different usages:
##'
##' (1) if \code{type == 'body'}, \code{upstream} and \code{downstream} can use to extend
##' the flank of body region.
##'
##' (2) if \code{type == 'start_site'/'end_site'}, \code{upstream} and \code{downstream} refer to
##' the upstream and downstream of the start_site or the end_site.
##'
##' \code{weightCol} refers to column in peak file. This column acts as a weight value. Details
##' see \url{https://github.com/YuLab-SMU/ChIPseeker/issues/15}
##'
##' \code{nbin} refers to the number of bins. \code{getTagMatrix()} provide a binning method
##' to get the tag matrix.
##'
##' There are two ways input a list of window.
##'
##' (1) Users can input a list of self-made granges objects
##'
##' (2) Users can input a list of \code{by} and only one \code{type}. In this way,
##' \code{plotPeakProf_MultiWindows()} can made a list of window from txdb object based on \code{by} and \code{type}.
##'
##' Warning:
##'
##' (1) All of these window should be the same type. It means users can only
##' compare a list of "start site"/"end site"/"body region" with the same upstream
##' and downstream.
##'
##' (2) So it will be only one \code{type} and several \code{by}.
##'
##' (3) Users can make window by txdb object or self-made granges object. Users can only
##' choose one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR' or 'UTR' in the
##' way of using txdb object. User can input any \code{by} in the way of using
##' self-made granges object.
##'
##' (4) Users can mingle the \code{by} designed for the two ways. \code{plotPeakProf_MultiWindows} can
##' accpet the hybrid \code{by}. But the above rules should be followed.
##'
##' \url{https://github.com/YuLab-SMU/ChIPseeker/issues/189}
##'
##' @title plotPeakProf_MultiWindows
##'
##' @param tagMatrix matrix or list of matrices, optional. Pre-computed tag matrix
##'   (from \code{getTagMatrix()}). If provided, peaks are not required. If NULL,
##'   tag matrix will be generated from peaks. Default is NULL
##' @param peak peak file (BED format) or GRanges object, or a list of peak
##'   files/GRanges objects for comparing multiple datasets. Required if
##'   \code{tagMatrix} is NULL
##' @param weightCol character, name of a metadata column in the GRanges object
##'   to use as weights for coverage calculation. If NULL, all peaks have equal
##'   weight. Default is NULL
##' @param TxDb TxDb or EnsDb annotation object, or a list of GRanges objects
##'   containing user-defined genomic features. If a list of GRanges is provided,
##'   each element should correspond to a window defined by \code{by}. Default is NULL
##' @param upstream numeric, upstream distance in base pairs. Interpretation
##'   depends on \code{type}: for "start_site"/"end_site", it's the distance
##'   upstream of the site; for "body", it's the extension beyond the feature body.
##'   Can be NULL if using binning mode. Default varies by function
##' @param downstream numeric, downstream distance in base pairs. Interpretation
##'   depends on \code{type}: for "start_site"/"end_site", it's the distance
##'   downstream of the site; for "body", it's the extension beyond the feature
##'   body. Can be NULL if using binning mode. Default varies by function
##' @param by character vector, feature(s) of interest. If using TxDb: one of
##'   'gene', 'transcript', 'exon', 'intron', '3UTR', '5UTR', 'UTR'. If using
##'   self-made GRanges: can be any label (used only for plot labels). For
##'   \code{plotPeakProf_MultiWindows}, can be a vector of length > 1 to compare
##'   multiple windows
##' @param type character, one of "start_site", "end_site", or "body". Determines
##'   which part of the feature to plot: transcription start site, transcription
##'   end site, or the entire feature body. For \code{plotPeakProf_MultiWindows},
##'   must be a single value (all windows must have the same type)
##' @param windows_name character vector, optional names for each window (used
##'   in \code{plotPeakProf_MultiWindows}). These names appear as labels in the
##'   plot. If NULL and \code{TxDb} is a named list, uses names from \code{TxDb}.
##'   If NULL and \code{by} is provided, uses \code{by} values. Default is NULL
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param conf numeric, confidence level (e.g., 0.95 for 95% CI). If provided,
##'   confidence intervals are calculated and displayed as ribbons around the
##'   average profile. If NA or missing, no confidence intervals are shown.
##'   Default varies by function
##' @param facet character, one of "none", "row", or "column". Controls how
##'   multiple datasets are displayed: "none" overlays them on the same plot,
##'   "row" creates row facets, "column" creates column facets. Default is "row"
##'   for \code{plotPeakProf}, "none" for others
##' @param free_y logical, if TRUE and faceting is used, y-axis scales
##'   independently for each facet. Default is TRUE
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param nbin integer, number of bins for binning mode. If provided, regions
##'   are divided into this many equal-sized bins instead of using absolute
##'   positions. Useful for comparing features of different lengths. If NULL,
##'   uses absolute positions. Default is NULL
##' @param ignore_strand logical, whether to ignore strand information when
##'   calculating coverage. If FALSE, peaks and features are matched by strand.
##'   Default is FALSE
##' @param ... additional arguments passed to ggplot2 functions for plot
##'   customization
##' @return ggplot object
##' @importFrom methods is
##' @importFrom methods as
##' @importFrom methods missingArg
##' @importFrom methods new
##' @export
plotPeakProf <- function(tagMatrix = NULL,
                         peak,
                         upstream,
                         downstream,
                         conf,
                         by,
                         type,
                         windows_name = NULL,
                         weightCol = NULL,
                         TxDb = NULL,
                         xlab = "Genomic Region (5'->3')",
                         ylab = "Peak Count Frequency",
                         facet = "row",
                         free_y = TRUE,
                         verbose = TRUE,
                         nbin = NULL,
                         ignore_strand = FALSE,
                         ...){

  if(is.null(tagMatrix)){

    conf <- if(missingArg(conf)) NA else conf
    upstream <- if(missingArg(upstream)) NULL else upstream
    downstream <- if(missingArg(downstream)) NULL else downstream

    if(length(by) == 1){

      plotPeakProf2(peak = peak,
                    upstream = upstream,
                    downstream = downstream,
                    conf = conf,
                    by = by,
                    type = type,
                    weightCol = weightCol,
                    TxDb = TxDb,
                    xlab = xlab,
                    ylab = ylab,
                    facet = facet,
                    free_y = free_y,
                    verbose = verbose,
                    nbin = nbin,
                    ignore_strand = ignore_strand,
                    ...)

    }else{

      if(is.null(windows_name) && !is.null(names(TxDb)))
        windows_name <- names(TxDb)

      plotPeakProf_MultiWindows(peak = peak,
                                upstream = upstream,
                                downstream = downstream,
                                conf = conf,
                                by = by,
                                type = type,
                                windows_name = windows_name,
                                weightCol = weightCol,
                                TxDb = TxDb,
                                xlab = xlab,
                                ylab = ylab,
                                facet = facet,
                                free_y = free_y,
                                verbose = verbose,
                                nbin = nbin,
                                ignore_strand = ignore_strand,
                                ...)

    }

  }else{

    if(is(tagMatrix, "list")){
      upstream <- attr(tagMatrix[[1]], 'upstream')
      downstream <- attr(tagMatrix[[1]], 'downstream')
      label <- attr(tagMatrix[[1]], 'label')
      attr(tagMatrix, 'type') <- attr(tagMatrix[[1]], 'type')
      attr(tagMatrix, 'is.binning') <- attr(tagMatrix[[1]], 'is.binning')

    }else{
      upstream <- attr(tagMatrix, 'upstream')
      downstream <- attr(tagMatrix, 'downstream')
      label <- attr(tagMatrix, 'label')
    }


    if(attr(tagMatrix, 'is.binning')){

      if (!(missingArg(conf) || is.na(conf))){

        plotAvgProf.binning(tagMatrix = tagMatrix,
                            xlab = xlab,
                            ylab = ylab,
                            conf = conf,
                            facet = facet,
                            free_y = free_y,
                            upstream = upstream,
                            downstream = downstream,
                            label = label,
                            ...)

      }else{

        plotAvgProf.binning(tagMatrix = tagMatrix,
                            xlab = xlab,
                            ylab = ylab,
                            facet = facet,
                            free_y = free_y,
                            upstream = upstream,
                            downstream = downstream,
                            label = label,
                            ...)

      }


    }else{

      xlim <- c(-upstream, downstream)

      if (!(missingArg(conf) || is.na(conf))){

        plotAvgProf (tagMatrix = tagMatrix,
                     xlim = xlim,
                     xlab = xlab,
                     ylab = ylab,
                     conf = conf,
                     facet = facet,
                     free_y = free_y,
                     origin_label = label,
                     ...)

      }else{

        plotAvgProf (tagMatrix = tagMatrix,
                     xlim = xlim,
                     xlab = xlab,
                     ylab = ylab,
                     facet = facet,
                     free_y = free_y,
                     origin_label = label,
                     ...)

      }


    }
  }

}


##' plot the profile of peaks
##'
##'
##' Plot average peak profile
##'
##' This function plots the average profile of peaks around genomic features
##' (e.g., TSS, gene body) as a line plot, optionally with confidence intervals.
##'
##' @description
##' The function creates a line plot showing the average signal intensity of
##' peaks across genomic regions. It can plot single or multiple peak sets,
##' with optional confidence intervals (error bars or ribbons) and faceting
##' options for comparing datasets.
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Calculates average signal across all features for each position
##'   \item Optionally calculates confidence intervals (if \code{conf} is provided)
##'   \item Creates a line plot with optional ribbon for confidence intervals
##'   \item Supports faceting for multiple datasets (row, column, or none)
##'   \item Allows independent y-axis scaling per facet
##' }
##'
##' For multiple datasets (list input), the function automatically creates
##' separate facets or overlays them based on the \code{facet} parameter.
##'
##' @param tagMatrix matrix or list of matrices. Each matrix should have features
##'   as rows and genomic positions as columns. Tag matrices are typically obtained
##'   from \code{getTagMatrix()}. If a list is provided, each element should be
##'   a named tag matrix
##' @param xlim numeric vector of length 2, x-axis limits (start, end positions).
##'   Should match the number of columns in tagMatrix. Default example: c(-3000, 3000)
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param conf numeric, confidence level (e.g., 0.95 for 95% CI). If provided,
##'   confidence intervals are calculated and displayed as ribbons. If NA or missing,
##'   no confidence intervals are shown
##' @param facet character, one of "none", "row", or "column". Controls how multiple
##'   datasets are displayed: "none" overlays them, "row" creates row facets,
##'   "column" creates column facets. Default is "none"
##' @param free_y logical, if TRUE and faceting is used, y-axis scales independently
##'   for each facet. Default is TRUE
##' @param origin_label character, label for the center point (e.g., "TSS", "TES").
##'   A vertical line is drawn at position 0 with this label. Default is "TSS"
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param ... additional arguments passed to ggplot2 functions
##' @return A ggplot2 object showing the average peak profile as a line plot,
##'   with optional confidence intervals and faceting
##' @author G Yu; Y Yan
##' @export
plotAvgProf <- function(tagMatrix, xlim,
                        xlab="Genomic Region (5'->3')",
                        ylab = "Peak Count Frequency",
                        conf,
                        facet="none",
                        free_y = TRUE,
                        origin_label = "TSS",
                        verbose = TRUE,
                        ...) {

  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)

  if (verbose) {
      cat(">> plotting figure...\t\t\t",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }

  conf <- if(missingArg(conf)) NA else conf

  if (!(missingArg(conf) || is.na(conf))){
    p <- plotAvgProf.internal(tagMatrix = tagMatrix,
                              conf = conf,
                              xlim = xlim,
                              xlab = xlab,
                              ylab = ylab,
                              facet = facet,
                              free_y = free_y,
                              origin_label = origin_label,
                              ...)
  } else {
    p <- plotAvgProf.internal(tagMatrix,
                              xlim = xlim,
                              xlab = xlab,
                              ylab = ylab,
                              facet = facet,
                              free_y = free_y,
                              origin_label = origin_label,
                              ...)
  }
  return(p)
}


##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_ribbon
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 facet_grid
##' @noRd
plotAvgProf.internal <- function(tagMatrix, conf,
                                 xlim = c(-3000,3000),
                                 xlab = "Genomic Region (5'->3')",
                                 ylab = "Peak Count Frequency",
                                 facet="none",
                                 free_y = TRUE,
                                 origin_label,
                                 ...) {

  listFlag <- FALSE
  if (is(tagMatrix, "list")) {
    if ( is.null(names(tagMatrix)) ) {
      nn <- paste0("peak", seq_along(tagMatrix))
      warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
      names(tagMatrix) <- nn
      ## stop("tagMatrix should be a named list...")
    }
    listFlag <- TRUE
  }

  if ( listFlag ) {
    facet <- match.arg(facet, c("none", "row", "column"))
    if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix[[1]]) ) {
      stop("please specify appropreate xcoordinations...")
    }
  } else {
    if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix) ) {
      stop("please specify appropreate xcoordinations...")
    }
  }

  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)
  ##
  conf <- if(missingArg(conf)) NA else conf

  pos <- value <- .id <- Lower <- Upper <- NULL

  if ( listFlag ) {
    tagCount <- lapply(tagMatrix, function(x) getTagCount(x, xlim = xlim, conf = conf, ...))
    tagCount <- list_to_dataframe(tagCount)
    tagCount$.id <- factor(tagCount$.id, levels=names(tagMatrix))
    p <- ggplot(tagCount, aes(pos, group=.id, color=.id))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = .id),
                           linetype = 0, alpha = 0.2)
    }
  } else {
    tagCount <- getTagCount(tagMatrix, xlim = xlim, conf = conf, ...)
    p <- ggplot(tagCount, aes(pos))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper),
                           linetype = 0, alpha = 0.2)
    }
  }

  p <- p + geom_line(aes(y = value))

  if ( 0 > xlim[1] && 0 < xlim[2] ) {
    p <- p + geom_vline(xintercept=0,
                        linetype="longdash")
    p <- p + scale_x_continuous(breaks=c(xlim[1], floor(xlim[1]/2),
                                         0,
                                         floor(xlim[2]/2), xlim[2]),
                                labels=c(paste0(xlim[1],"bp"), paste0(floor(xlim[1]/2),"bp"),
                                         origin_label,
                                         paste0(floor(xlim[2]/2),"bp"), paste0(xlim[2], "bp")))
  }

  if (listFlag) {
    cols <- getCols(length(tagMatrix))
    p <- p + scale_color_manual(values=cols)
    if (facet == "row") {
      if (free_y) {
        p <- p + facet_grid(.id ~ ., scales = "free_y")
      } else {
        p <- p + facet_grid(.id ~ .)
      }
    } else if (facet == "column") {
      if (free_y) {
        p <-  p + facet_grid(. ~ .id, scales = "free_y")
      } else {
        p <-  p + facet_grid(. ~ .id)
      }
    }
  }
  p <- p+xlab(xlab)+ylab(ylab)
  p <- p + theme_bw() + theme(legend.title=element_blank())
  if(facet != "none") {
    p <- p + theme(legend.position="none")
  }
  return(p)
}

##' plot the profile of peaks that align to flank sequences of TSS
##'
##' This function is the old function of \code{plotPeakProf2}. It can
##' only plot the start site region of gene.
##'
##' @title plotAvgProf2
##' @param peak peak file (BED format) or GRanges object, or a list of peak
##'   files/GRanges objects for comparing multiple datasets
##' @param weightCol character, name of a metadata column in the GRanges object
##'   to use as weights for coverage calculation. If NULL, all peaks have equal
##'   weight. Default is NULL
##' @param TxDb TxDb or EnsDb annotation object. Default is NULL
##' @param upstream numeric, upstream distance in base pairs from the start site.
##'   Default is 1000
##' @param downstream numeric, downstream distance in base pairs from the start
##'   site. Default is 1000
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param conf numeric, confidence level (e.g., 0.95 for 95% CI). If provided,
##'   confidence intervals are displayed. If NA or missing, no confidence
##'   intervals are shown
##' @param facet character, one of "none", "row", or "column". Controls how
##'   multiple datasets are displayed. Default is "none"
##' @param free_y logical, if TRUE and faceting is used, y-axis scales
##'   independently for each facet. Default is TRUE
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param ignore_strand logical, whether to ignore strand information when
##'   calculating coverage. Default is FALSE
##' @param ... additional arguments passed to ggplot2 functions
##' @return ggplot object
##' @export
##' @author G Yu, Ming L
plotAvgProf2 <- function(peak, weightCol = NULL, TxDb = NULL,
                         upstream = 1000, downstream = 1000,
                         xlab = "Genomic Region (5'->3')",
                         ylab = "Peak Count Frequency",
                         conf,
                         facet = "none",
                         free_y = TRUE,
                         verbose = TRUE,
                         ignore_strand = FALSE,
                         ...) {

  plotPeakProf2(peak = peak,
                upstream = upstream,
                downstream = downstream,
                conf,
                by = "gene",
                type = "start_site",
                weightCol = weightCol,
                TxDb = TxDb,
                xlab = xlab,
                ylab = ylab,
                facet = facet,
                free_y = free_y,
                verbose = verbose,
                ignore_strand = ignore_strand,
                ...)

}

##' Plot average peak profile using binning method
##'
##' This function plots average peak profiles when tag matrices were generated
##' using binning (dividing regions into equal-sized bins). It handles variable
##' region lengths by normalizing to a fixed number of bins.
##'
##' @description
##' The binning method is useful for comparing peaks across features of different
##' lengths (e.g., genes of varying sizes). Instead of absolute positions, the
##' x-axis represents relative positions within each feature, divided into bins.
##'
##' @details
##' The function supports flexible upstream/downstream specifications:
##' \itemize{
##'   \item \code{rel()} objects: specify as percentage of feature length
##'         (e.g., \code{rel(0.2)} for 20% of feature length)
##'   \item Integers: specify as absolute base pairs
##'   \item NULL: no extension (gene body only)
##' }
##'
##' @param tagMatrix matrix or list of matrices generated with binning (from
##'   \code{getTagMatrix()} with \code{nbin} parameter). Each matrix has features
##'   as rows and bins as columns
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param conf numeric, confidence level for confidence intervals. If NA or
##'   missing, no confidence intervals are shown
##' @param facet character, one of "none", "row", or "column" for faceting.
##'   Default is "none"
##' @param free_y logical, if TRUE, y-axis scales independently per facet.
##'   Default is TRUE
##' @param upstream numeric, \code{rel()} object, or NULL. Upstream extension:
##'   \code{rel(0.2)} = 20% of feature length, integer = absolute bp, NULL = no
##'   extension. Default is NULL
##' @param downstream numeric, \code{rel()} object, or NULL. Downstream extension:
##'   \code{rel(0.2)} = 20% of feature length, integer = absolute bp, NULL = no
##'   extension. Default is NULL
##' @param label character, label for the center point (e.g., "TSS", "TES").
##'   A vertical line marks this position
##' @param ... additional arguments passed to ggplot2 functions
##' @return A ggplot2 object showing the average peak profile using binning,
##'   with x-axis representing relative positions within features
##' @importFrom ggplot2 rel
##' @noRd
plotAvgProf.binning <- function(tagMatrix,
                                xlab = "Genomic Region (5'->3')",
                                ylab = "Peak Count Frequency",
                                conf,
                                facet ="none",
                                free_y = TRUE,
                                upstream = NULL,
                                downstream = NULL,
                                label,
                                ...) {

  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)
  conf <- if(missingArg(conf)) NA else conf

  if (!(missingArg(conf) || is.na(conf))){
    p <- plotAvgProf.binning.internal(tagMatrix ,
                                      conf = conf,
                                      xlab = xlab,
                                      ylab = ylab,
                                      facet = facet,
                                      free_y = free_y,
                                      upstream = upstream,
                                      downstream = downstream,
                                      label = label,
                                      ...)
  } else {
    p <- plotAvgProf.binning.internal(tagMatrix ,
                                      xlab = xlab,
                                      ylab = ylab,
                                      facet = facet,
                                      free_y = free_y,
                                      upstream = upstream,
                                      downstream = downstream,
                                      label = label,
                                      ...)
  }
  return(p)
}


##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_ribbon
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 facet_grid
##' @importFrom ggplot2 rel
##' @noRd
plotAvgProf.binning.internal <- function(tagMatrix,
                                         conf,
                                         xlab = "Genomic Region (5'->3')",
                                         ylab = "Peak Count Frequency",
                                         facet="none",
                                         free_y = TRUE,
                                         upstream = NULL,
                                         downstream = NULL,
                                         label,
                                         ...) {

  listFlag <- FALSE
  if (is(tagMatrix, "list")) {
    if ( is.null(names(tagMatrix )) ) {
      nn <- paste0("peak", seq_along(tagMatrix ))
      warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
      names(tagMatrix) <- nn
      ## stop("tagMatrix should be a named list...")
    }
    listFlag <- TRUE
  }

  if(listFlag){
    nbin <- dim(tagMatrix[[1]])[2]
  }else{
    nbin <- dim(tagMatrix)[2]
  }
  xlim <- c(1,nbin)

  if ( listFlag ) {
    facet <- match.arg(facet, c("none", "row", "column"))
  }

  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)
  ##
  conf <- if(missingArg(conf)) NA else conf

  pos <- value <- .id <- Lower <- Upper <- NULL

  if ( listFlag ) {
    tagCount <- lapply(tagMatrix , function(x) getTagCount(x, xlim = xlim, conf = conf, ...))
    tagCount <- list_to_dataframe(tagCount)
    tagCount$.id <- factor(tagCount$.id, levels=names(tagMatrix ))
    p <- ggplot(tagCount, aes(pos, group=.id, color=.id))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = .id),
                           linetype = 0, alpha = 0.2)
    }
  } else {
    tagCount <- getTagCount(tagMatrix , xlim = xlim, conf = conf, ...)
    p <- ggplot(tagCount, aes(pos))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper),
                           linetype = 0, alpha = 0.2)
    }
  }

  p <- p + geom_line(aes(y = value))

  ## x_scale for genebody
  if(attr(tagMatrix, 'type') == 'body'){
    ## x_scale for gene body with no flank extension
    if(is.null(upstream)){
      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*0.25),
                                           floor(nbin*0.5),
                                           floor(nbin*0.75),
                                           nbin),
                                  labels=c(label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2]))
    }


    ## x_scale for flank extension by relative value
    if(inherits(upstream, 'rel')){

      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*(as.numeric(upstream)*100/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+25)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+50)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+75)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+100)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           nbin),
                                  labels=c(paste0("-",as.numeric(upstream)*100,"%"),
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0("+",as.numeric(downstream)*100,"%")))
      p <- p + geom_vline(xintercept=floor(nbin*(as.numeric(upstream)*100/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                          linetype="longdash")

      p <- p + geom_vline(xintercept=floor(nbin*((as.numeric(upstream)*100+100)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                          linetype="longdash")
    }

    ## x_scale for flank extension by absolute value
    if(!is.null(upstream) & !inherits(upstream, 'rel')){

      upstreamPer <- floor(upstream/1000)*0.1
      downstreamPer <- floor(downstream/1000)*0.1

      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.25)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.5)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.75)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))),
                                           nbin),
                                  labels=c(paste0("-",upstream,"bp"),
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0(downstream,"bp")))
      p <- p + geom_vline(xintercept=floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
                          linetype="longdash")

      p <- p + geom_vline(xintercept=floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))),
                          linetype="longdash")
    }
  }


  ## x_scale for start region
  if(attr(tagMatrix, 'type') != 'body'){

    p <- p + scale_x_continuous(breaks=c(1,
                                         floor(nbin*0.25),
                                         floor(nbin*0.5),
                                         floor(nbin*0.75),
                                         nbin),
                                labels=c(paste0("-",upstream,"bp"),
                                         paste0("-",floor(upstream*0.5),"bp"),
                                         label,
                                         paste0(floor(downstream*0.5),"bp"),
                                         paste0(downstream,"bp")))

    p <- p + geom_vline(xintercept=floor(nbin*0.5),
                        linetype="longdash")
  }


  if (listFlag) {
    cols <- getCols(length(tagMatrix))
    p <- p + scale_color_manual(values=cols)
    if (facet == "row") {
      if (free_y) {
        p <- p + facet_grid(.id ~ ., scales = "free_y")
      } else {
        p <- p + facet_grid(.id ~ .)
      }
    } else if (facet == "column") {
      if (free_y) {
        p <-  p + facet_grid(. ~ .id, scales = "free_y")
      } else {
        p <-  p + facet_grid(. ~ .id)
      }
    }
  }
  p <- p+xlab(xlab)+ylab(ylab)
  p <- p + theme_bw() + theme(legend.title=element_blank())
  if(facet != "none") {
    p <- p + theme(legend.position="none")
  }
  return(p)
}


##' plot the profile of peaks automatically
##'
##' \code{peak} stands for the peak file.
##'
##' \code{by} the features of interest.
##'
##' (1) if users use \code{txdb}, \code{by} can be one of 'gene', 'transcript', 'exon',
##' 'intron' , '3UTR' , '5UTR', 'UTR'. These features can be obtained by functions from txdb object.
##'
##' (2) if users use self-made granges object, \code{by} can be everything. Because this \code{by}
##' will not pass to functions to get features, which is different from the case of using
##' txdb object. This \code{by} is only used to made labels showed in picture.
##'
##' \code{type} means the property of the region. one of the "start site",
##' "end site" and "body".
##'
##' \code{upstream} and \code{downstream} parameter have different usages:
##'
##' (1) if \code{type == 'body'}, \code{upstream} and \code{downstream} can use to extend
##' the flank of body region.
##'
##' (2) if \code{type == 'start_site'/'end_site'}, \code{upstream} and \code{downstream} refer to
##' the upstream and downstream of the start_site or the end_site.
##'
##' \code{weightCol} refers to column in peak file. This column acts as a weight vaule. Details
##' see \url{https://github.com/YuLab-SMU/ChIPseeker/issues/15}
##'
##' \code{nbin} refers to the number of bins, providing a binning method
##' to get the tag matrix.
##'
##' \code{TxDb} parameter can accept txdb object.
##' But many regions can not be obtained by txdb object. In this case,
##' Users can provide self-made granges served the same role
##' as txdb object and pass to \code{TxDb} object.
##'
##' \code{plotPeakProf2()} is different from the \code{plotPeakProf()}. \code{plotPeakProf2()} do not
##' need to provide \code{window} parameter, which means \code{plotPeakProf2()} will call relevent
##' functions to make \code{window} automatically.
##'
##' @title plotPeakProf2
##' @param peak peak file (BED format) or GRanges object, or a list of peak
##'   files/GRanges objects for comparing multiple datasets
##' @param upstream numeric, upstream distance in base pairs. Interpretation
##'   depends on \code{type}: for "start_site"/"end_site", it's the distance
##'   upstream of the site; for "body", it's the extension beyond the feature body
##' @param downstream numeric, downstream distance in base pairs. Interpretation
##'   depends on \code{type}: for "start_site"/"end_site", it's the distance
##'   downstream of the site; for "body", it's the extension beyond the feature
##'   body
##' @param conf numeric, confidence level (e.g., 0.95 for 95% CI). If provided,
##'   confidence intervals are displayed. If NA or missing, no confidence
##'   intervals are shown
##' @param by character, feature of interest. If using TxDb: one of 'gene',
##'   'transcript', 'exon', 'intron', '3UTR', '5UTR', 'UTR'. If using self-made
##'   GRanges: can be any label (e.g., "enhancer", "promoter") used only for plot
##'   labels
##' @param type character, one of "start_site", "end_site", or "body". Determines
##'   which part of the feature to plot
##' @param weightCol character, name of a metadata column in the GRanges object
##'   to use as weights for coverage calculation. If NULL, all peaks have equal
##'   weight. Default is NULL
##' @param TxDb TxDb or EnsDb annotation object, or a GRanges object containing
##'   user-defined genomic features. Default is NULL
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param facet character, one of "none", "row", or "column". Controls how
##'   multiple datasets are displayed. Default is "none"
##' @param free_y logical, if TRUE and faceting is used, y-axis scales
##'   independently for each facet. Default is TRUE
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param nbin integer, number of bins for binning mode. If provided, regions
##'   are divided into this many equal-sized bins. If NULL, uses absolute
##'   positions. Default is NULL
##' @param ignore_strand logical, whether to ignore strand information when
##'   calculating coverage. Default is FALSE
##' @param ... additional arguments passed to ggplot2 functions
##' @return ggplot object
##' @export
##' @author G Yu, Ming Li
plotPeakProf2 <- function(peak,
                          upstream,
                          downstream,
                          conf,
                          by,
                          type,
                          weightCol = NULL,
                          TxDb = NULL,
                          xlab = "Genomic Region (5'->3')",
                          ylab = "Peak Count Frequency",
                          facet = "none",
                          free_y = TRUE,
                          verbose = TRUE,
                          nbin = NULL,
                          ignore_strand = FALSE,
                          ...){

  conf <- if(missingArg(conf)) NA else conf
  upstream <- if(missingArg(upstream)) NULL else upstream
  downstream <- if(missingArg(downstream)) NULL else downstream

  if ( is(peak, "list") ) {
    tagMatrix <- lapply(peak, getTagMatrix,
                        upstream = upstream,
                        downstream = downstream,
                        type = type,
                        TxDb = TxDb,
                        by = by,
                        weightCol = weightCol,
                        nbin = nbin,
                        verbose = verbose,
                        ignore_strand = ignore_strand)
  } else {
    tagMatrix <- getTagMatrix(peak = peak,
                              upstream = upstream,
                              downstream = downstream,
                              type = type,
                              by = by,
                              TxDb = TxDb,
                              weightCol = weightCol,
                              nbin = nbin,
                              verbose = verbose,
                              ignore_strand = ignore_strand)
  }


  if (!(missingArg(conf) || is.na(conf))){
    p <- plotPeakProf(tagMatrix = tagMatrix,
                      conf = conf,
                      xlab = xlab,
                      ylab = ylab,
                      facet = facet,
                      free_y = free_y,
                      ...)

  } else {
    p <- plotPeakProf(tagMatrix = tagMatrix,
                      xlab = xlab,
                      ylab = ylab,
                      facet= facet,
                      free_y = free_y,
                      ...)
  }
  return(p)

}


##' plot the profile of peaks in two or more windows
##'
##'
##' This function comes from \url{https://github.com/YuLab-SMU/ChIPseeker/issues/189}
##'`
##' \code{plotPeakProf_MultiWindows()} is almost the same as \code{plotPeakProf2()}, having
##' the main difference of accepting two or more granges objects. Accepting more
##' granges objects can help compare the same peaks in different windows.
##'
##' \code{TxDb} parameter can accept txdb object.
##' But many regions can not be obtained by txdb object. In this case,
##' Users can provide self-made granges served the same role
##' as txdb object and pass to \code{TxDb} object.
##'
##' \code{by} the features of interest.
##'
##' (1) if users use \code{txdb}, \code{by} can be one of 'gene', 'transcript', 'exon',
##' 'intron' , '3UTR' , '5UTR', 'UTR'. These features can be obtained by functions from txdb object.
##'
##' (2) if users use self-made granges object, \code{by} can be everything. Because this \code{by}
##' will not pass to functions to get features, which is different from the case of using
##' txdb object. This \code{by} is only used to made labels showed in picture.
##'
##' \code{type} means the property of the region. one of the "start site",
##' "end site" and "body".
##'
##' \code{upstream} and \code{downstream} parameter have different usages:
##'
##' (1) if \code{type == 'body'}, \code{upstream} and \code{downstream} can use to extend
##' the flank of body region.
##'
##' (2) if \code{type == 'start_site'/'end_site'}, \code{upstream} and \code{downstream} refer to
##' the upstream and downstream of the start_site or the end_site.
##'
##' \code{weightCol} refers to column in peak file. This column acts as a weight value. Details
##' see \url{https://github.com/YuLab-SMU/ChIPseeker/issues/15}
##'
##' \code{nbin} refers to the number of bins. \code{getTagMatrix()} provide a binning method
##' to get the tag matrix.
##'
##' There are two ways input a list of window.
##'
##' (1) Users can input a list of self-made granges objects
##'
##' (2) Users can input a list of \code{by} and only one \code{type}. In this way,
##' \code{plotPeakProf_MultiWindows()} can made a list of window from txdb object based on \code{by} and \code{type}.
##'
##' Warning:
##'
##' (1) All of these window should be the same type. It means users can only
##' compare a list of "start site"/"end site"/"body region" with the same upstream
##' and downstream.
##'
##' (2) So it will be only one \code{type} and several \code{by}.
##'
##' (3) Users can make window by txdb object or self-made granges object. Users can only
##' choose one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR' or 'UTR' in the
##' way of using txdb object. User can input any \code{by} in the way of using
##' self-made granges object.
##'
##' (4) Users can mingle the \code{by} designed for the two ways. \code{plotPeakProf_MultiWindows} can
##' accpet the hybrid \code{by}. But the above rules should be followed.
##'
##'
##' @title plotPeakProf_MultiWindows
##' @param peak peak file (BED format) or GRanges object, or a list of peak
##'   files/GRanges objects for comparing multiple datasets
##' @param upstream numeric, upstream distance in base pairs. Interpretation
##'   depends on \code{type}: for "start_site"/"end_site", it's the distance
##'   upstream of the site; for "body", it's the extension beyond the feature body
##' @param downstream numeric, downstream distance in base pairs. Interpretation
##'   depends on \code{type}: for "start_site"/"end_site", it's the distance
##'   downstream of the site; for "body", it's the extension beyond the feature
##'   body
##' @param conf numeric, confidence level (e.g., 0.95 for 95% CI). If provided,
##'   confidence intervals are displayed. If NA or missing, no confidence
##'   intervals are shown
##' @param by character vector, features of interest (length > 1). If using TxDb:
##'   one or more of 'gene', 'transcript', 'exon', 'intron', '3UTR', '5UTR', 'UTR'.
##'   If using self-made GRanges: can be any labels. All windows must have the
##'   same \code{type}
##' @param type character, one of "start_site", "end_site", or "body". Must be
##'   a single value (all windows must have the same type). Determines which part
##'   of the feature to plot
##' @param windows_name character vector, optional names for each window. These
##'   names appear as labels in the plot. If NULL and \code{TxDb} is a named list,
##'   uses names from \code{TxDb}. If NULL, uses \code{by} values. Length must match
##'   \code{by}. Default is NULL
##' @param weightCol character, name of a metadata column in the GRanges object
##'   to use as weights for coverage calculation. If NULL, all peaks have equal
##'   weight. Default is NULL
##' @param TxDb TxDb or EnsDb annotation object, or a list of GRanges objects
##'   containing user-defined genomic features. If a list, each element should
##'   correspond to a window defined by \code{by}. Default is NULL
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param facet character, one of "none", "row", or "column". Controls how
##'   multiple windows are displayed. Default is "row"
##' @param free_y logical, if TRUE and faceting is used, y-axis scales
##'   independently for each facet. Default is TRUE
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param nbin integer, number of bins for binning mode. If provided, regions
##'   are divided into this many equal-sized bins. If NULL, uses absolute
##'   positions. Default is NULL
##' @param ignore_strand logical, whether to ignore strand information when
##'   calculating coverage. Default is FALSE
##' @param ... additional arguments passed to ggplot2 functions
##' @return ggplot object
##' @noRd
plotPeakProf_MultiWindows <- function(peak,
                                      upstream,
                                      downstream,
                                      conf,
                                      by,
                                      type,
                                      windows_name = NULL,
                                      weightCol = NULL,
                                      TxDb = NULL,
                                      xlab = "Genomic Region (5'->3')",
                                      ylab = "Peak Count Frequency",
                                      facet = "row",
                                      free_y = TRUE,
                                      verbose = TRUE,
                                      nbin = NULL,
                                      ignore_strand = FALSE,
                                      ...){

  conf <- if(missingArg(conf)) NA else conf
  upstream <- if(missingArg(upstream)) NULL else upstream
  downstream <- if(missingArg(downstream)) NULL else downstream

  ## check type
  if(length(type) != 1){
    stop("It should be only one type...")
  }

  ## make the window name
  if (is.null(windows_name)) {
    nn <- by
    warning("set the name automatically to ", paste(nn, collapse=' '))
    windows_name <- nn
  }else{
    if (length(windows_name) != length(by)) {
      stop("the length of the window name and the by should be equal...")
    }
  }


  if ( is(peak, "list") ) {
    tagMatrix <- lapply(peak, getTagMatrix2,
                        upstream=upstream,
                        downstream=downstream,
                        windows_name=windows_name,
                        type=type,
                        by=by,
                        TxDb=TxDb,
                        weightCol = weightCol,
                        nbin = nbin,
                        verbose = verbose,
                        ignore_strand= ignore_strand)
  } else {
    tagMatrix <- getTagMatrix2(peak=peak,
                               upstream=upstream,
                               downstream=downstream,
                               windows_name=windows_name,
                               type=type,
                               by=by,
                               TxDb=TxDb,
                               weightCol = weightCol,
                               nbin = nbin,
                               verbose = verbose,
                               ignore_strand= ignore_strand)
  }

  if (!(missingArg(conf) || is.na(conf))){
    p <- plotMultiProf(tagMatrix = tagMatrix,
                       conf = conf,
                       xlab = xlab,
                       ylab = ylab,
                       facet = facet,
                       free_y = free_y,
                       ...)

  } else {
    p <- plotMultiProf(tagMatrix = tagMatrix,
                       xlab = xlab,
                       ylab = ylab,
                       facet= facet,
                       free_y = free_y,
                       ...)
  }
  return(p)

}


##' Plot multiple peak profiles (internal function for plotPeakProf_MultiWindows)
##'
##' @param tagMatrix list of tag matrices, where each element corresponds to a
##'   different window. Tag matrices should have attributes 'upstream',
##'   'downstream', 'type', and 'is.binning'
##' @param conf numeric, confidence level (e.g., 0.95 for 95% CI). If provided,
##'   confidence intervals are displayed. If NA or missing, no confidence
##'   intervals are shown
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param facet character, one of "none", "row", or "column". Controls how
##'   multiple windows are displayed. Default is "none"
##' @param free_y logical, if TRUE and faceting is used, y-axis scales
##'   independently for each facet. Default is TRUE
##' @param ... additional arguments passed to ggplot2 functions
##' @noRd
plotMultiProf <- function(tagMatrix,
                          conf,
                          xlab="Genomic Region (5'->3')",
                          ylab = "Peak Count Frequency",
                          facet="none",
                          free_y = TRUE,
                          ...){


  if(is(tagMatrix[[1]][[1]],"matrix")){
    upstream <- attr(tagMatrix[[1]][[1]], 'upstream')
    downstream <- attr(tagMatrix[[1]][[1]], 'downstream')
    # attr(tagMatrix, 'type') <- attr(tagMatrix[[1]][[1]], 'type')
    # attr(tagMatrix, 'is.binning') <- attr(tagMatrix[[1]][[1]], 'is.binning')
    binFlag <- attr(tagMatrix[[1]][[1]], 'is.binning')
    type <- attr(tagMatrix[[1]][[1]], 'type')

  }else{
    upstream <- attr(tagMatrix[[1]], 'upstream')
    downstream <- attr(tagMatrix[[1]], 'downstream')
    binFlag <- attr(tagMatrix[[1]], 'is.binning')
    type <- attr(tagMatrix[[1]], 'type')
  }

  if(type == "body"){

    label <- c("SS","TS")

  }else if(type == "start_site"){

    label <- "SS"

  }else{

    label <- "TS"

  }


  if(binFlag){

    if (!(missingArg(conf) || is.na(conf))){

      plotMultiProf.binning(tagMatrix = tagMatrix,
                            xlab = xlab,
                            ylab = ylab,
                            conf = conf,
                            facet = facet,
                            free_y = free_y,
                            upstream = upstream,
                            downstream = downstream,
                            label = label,
                            ...)

    }else{

      plotMultiProf.binning(tagMatrix = tagMatrix,
                            xlab = xlab,
                            ylab = ylab,
                            facet = facet,
                            free_y = free_y,
                            upstream = upstream,
                            downstream = downstream,
                            label = label,
                            ...)
    }


  }else{

    xlim <- c(-upstream, downstream)

    if (!(missingArg(conf) || is.na(conf))){

      plotMultiProf.normal(tagMatrix = tagMatrix,
                           xlim = xlim,
                           xlab = xlab,
                           ylab = ylab,
                           conf = conf,
                           facet = facet,
                           free_y = free_y,
                           origin_label = label,
                           ...)

    }else{

      plotMultiProf.normal(tagMatrix = tagMatrix,
                           xlim = xlim,
                           xlab = xlab,
                           ylab = ylab,
                           facet = facet,
                           free_y = free_y,
                           origin_label = label,
                           ...)
    }
  }

}

##' Plot multiple peak profiles using normal (non-binning) method
##'
##' @param tagMatrix list of tag matrices, where each element corresponds to a
##'   different window. Tag matrices should have features as rows and genomic
##'   positions as columns
##' @param xlim numeric vector of length 2, x-axis limits (start, end positions).
##'   Should match the number of columns in tagMatrix
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param conf numeric, confidence level (e.g., 0.95 for 95% CI). If provided,
##'   confidence intervals are displayed. If NA or missing, no confidence
##'   intervals are shown
##' @param facet character, one of "none", "row", or "column". Controls how
##'   multiple windows are displayed. Default is "none"
##' @param free_y logical, if TRUE and faceting is used, y-axis scales
##'   independently for each facet. Default is TRUE
##' @param origin_label character vector, label(s) for the center point(s).
##'   For "body" type: c("SS", "TS"); for "start_site": "SS"; for "end_site": "TS".
##'   Default is "TSS"
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param ... additional arguments passed to ggplot2 functions
##' @noRd
plotMultiProf.normal <- function(tagMatrix, xlim,
                                 xlab="Genomic Region (5'->3')",
                                 ylab = "Peak Count Frequency",
                                 conf,
                                 facet="none",
                                 free_y = TRUE,
                                 origin_label = "TSS",
                                 verbose = TRUE,
                                 ...) {

  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)

  if (verbose) {
    cat(">> plotting figure...\t\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }

  conf <- if(missingArg(conf)) NA else conf

  if (!(missingArg(conf) || is.na(conf))){

    p <- plotMultiProf.normal.internal(tagMatrix = tagMatrix,
                                       conf = conf,
                                       xlim = xlim,
                                       xlab = xlab,
                                       ylab = ylab,
                                       facet = facet,
                                       free_y = free_y,
                                       origin_label = origin_label,
                                       ...)


  } else {

    p <- plotMultiProf.normal.internal(tagMatrix,
                                       xlim = xlim,
                                       xlab = xlab,
                                       ylab = ylab,
                                       facet = facet,
                                       free_y = free_y,
                                       origin_label = origin_label,
                                       ...)
  }
  return(p)
}

##' Plot multiple peak profiles using normal method (internal function)
##'
##' @param tagMatrix list of tag matrices, where each element corresponds to a
##'   different window. Tag matrices should have features as rows and genomic
##'   positions as columns
##' @param conf numeric, confidence level (e.g., 0.95 for 95% CI). If provided,
##'   confidence intervals are displayed. If NA or missing, no confidence
##'   intervals are shown
##' @param xlim numeric vector of length 2, x-axis limits (start, end positions).
##'   Should match the number of columns in tagMatrix. Default is c(-3000, 3000)
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param facet character, one of "none", "row", or "column". Controls how
##'   multiple windows are displayed. Default is "none"
##' @param free_y logical, if TRUE and faceting is used, y-axis scales
##'   independently for each facet. Default is TRUE
##' @param origin_label character vector, label(s) for the center point(s).
##'   For "body" type: c("SS", "TS"); for "start_site": "SS"; for "end_site": "TS"
##' @param ... additional arguments passed to ggplot2 functions
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_ribbon
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 facet_grid
##' @noRd
plotMultiProf.normal.internal <- function(tagMatrix, conf,
                                          xlim = c(-3000,3000),
                                          xlab = "Genomic Region (5'->3')",
                                          ylab = "Peak Count Frequency",
                                          facet="row",
                                          free_y = TRUE,
                                          origin_label,
                                          ...) {

  listFlag <- FALSE
  if (is.null(attr(tagMatrix[[1]],'upstream'))) {
    if ( is.null(names(tagMatrix)) ) {
      nn <- paste0("peak", seq_along(tagMatrix))
      warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
      names(tagMatrix) <- nn
      ## stop("tagMatrix should be a named list...")
    }
    listFlag <- TRUE
  }

  if ( listFlag ) {
    facet <- match.arg(facet, c("none", "row", "column"))
    if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix[[1]][[1]]) ) {
      stop("please specify appropreate xcoordinations...")
    }
  } else {
    if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix[[1]]) ) {
      stop("please specify appropreate xcoordinations...")
    }
  }

  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)
  ##
  conf <- if(missingArg(conf)) NA else conf

  pos <- value <- .id <- Lower <- Upper <- NULL

  if ( listFlag ) {

    tagCount <- lapply(as.list(names(tagMatrix)), function(x){

      tmp <- tagMatrix[[x]]
      tagCount_tmp <- lapply(as.list(names(tmp)),function(x){
        result <- getTagCount(tmp[[x]], xlim = xlim, conf = conf, ...)
        result$type <- x

        return(result)
      })
      tagCount_tmp <- list_to_dataframe(tagCount_tmp)
      return(tagCount_tmp)

    })

    names(tagCount) <- names(tagMatrix)
    tagCount <- list_to_dataframe(tagCount)
    tagCount$.id <- factor(tagCount$.id, levels=names(tagMatrix))
    p <- ggplot(tagCount, aes(pos, group=type, color=type))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = type),
                           linetype = 0, alpha = 0.2)
    }

  } else {

    tagCount <- lapply(as.list(names(tagMatrix)), function(x){

      result <- getTagCount(tagMatrix[[x]], xlim = xlim, conf = conf, ...)
      result$type <- x

      return(result)
    })

    tagCount <- do.call("rbind",tagCount)

    p <- ggplot(tagCount, aes(x = pos))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper,fill = type),
                           linetype = 0, alpha = 0.2)
    }
  }

  p <- p + geom_line(aes(y = value,color = type))

  if ( 0 > xlim[1] && 0 < xlim[2] ) {
    p <- p + geom_vline(xintercept=0,
                        linetype="longdash")
    p <- p + scale_x_continuous(breaks=c(xlim[1], floor(xlim[1]/2),
                                         0,
                                         floor(xlim[2]/2), xlim[2]),
                                labels=c(paste0(xlim[1],"bp"), paste0(floor(xlim[1]/2),"bp"),
                                         origin_label,
                                         paste0(floor(xlim[2]/2),"bp"), paste0(xlim[2], "bp")))
  }

  if (listFlag) {
    # cols <- getCols(length(tagMatrix[[1]]))
    # p <- p + scale_color_manual(values=cols)
    if (facet == "row") {
      if (free_y) {
        p <- p + facet_grid(.id ~ ., scales = "free_y")
      } else {
        p <- p + facet_grid(.id ~ .)
      }
    } else if (facet == "column") {
      if (free_y) {
        p <-  p + facet_grid(. ~ .id, scales = "free_y")
      } else {
        p <-  p + facet_grid(. ~ .id)
      }
    }
  }

  p <- p+xlab(xlab)+ylab(ylab)
  p <- p + theme_bw() + theme(legend.title=element_blank())

  # if(facet != "none") {
  #   p <- p + theme(legend.position="none")
  # }

  return(p)
}

##' internal function
##'
##' @param tagMatrix tagMatrix
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param upstream the upstream extension
##' @param downstream the downstream extension
##' @param label the label of the center
##' @param ... additional parameter
##' @noRd
plotMultiProf.binning <- function(tagMatrix,
                                  xlab = "Genomic Region (5'->3')",
                                  ylab = "Peak Count Frequency",
                                  conf,
                                  facet ="none",
                                  free_y = TRUE,
                                  upstream = NULL,
                                  downstream = NULL,
                                  label,
                                  ...) {

  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)
  conf <- if(missingArg(conf)) NA else conf

  if (!(missingArg(conf) || is.na(conf))){
    p <- plotMultiProf.binning.internal(tagMatrix ,
                                        conf = conf,
                                        xlab = xlab,
                                        ylab = ylab,
                                        facet = facet,
                                        free_y = free_y,
                                        upstream = upstream,
                                        downstream = downstream,
                                        label = label,
                                        ...)
  } else {
    p <- plotMultiProf.binning.internal(tagMatrix ,
                                        xlab = xlab,
                                        ylab = ylab,
                                        facet = facet,
                                        free_y = free_y,
                                        upstream = upstream,
                                        downstream = downstream,
                                        label = label,
                                        ...)
  }
  return(p)
}

##' Plot multiple peak profiles using binning method (internal function)
##'
##' @param tagMatrix list of tag matrices generated with binning, where each
##'   element corresponds to a different window. Tag matrices should have features
##'   as rows and bins as columns
##' @param conf numeric, confidence level for confidence intervals. If NA or
##'   missing, no confidence intervals are shown
##' @param xlab character, x-axis label. Default is "Genomic Region (5'->3')"
##' @param ylab character, y-axis label. Default is "Peak Count Frequency"
##' @param facet character, one of "none", "row", or "column" for faceting.
##'   Default is "none"
##' @param free_y logical, if TRUE, y-axis scales independently per facet.
##'   Default is TRUE
##' @param upstream numeric, \code{rel()} object, or NULL. Upstream extension:
##'   \code{rel(0.2)} = 20% of feature length, integer = absolute bp, NULL = no
##'   extension. Default is NULL
##' @param downstream numeric, \code{rel()} object, or NULL. Downstream extension:
##'   \code{rel(0.2)} = 20% of feature length, integer = absolute bp, NULL = no
##'   extension. Default is NULL
##' @param label character vector, label(s) for the center point(s). For "body"
##'   type: c("SS", "TS"); for "start_site": "SS"; for "end_site": "TS"
##' @param ... additional arguments passed to ggplot2 functions
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_ribbon
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 facet_grid
##' @importFrom ggplot2 rel
##' @noRd
plotMultiProf.binning.internal <- function(tagMatrix,
                                           conf,
                                           xlab = "Genomic Region (5'->3')",
                                           ylab = "Peak Count Frequency",
                                           facet="none",
                                           free_y = TRUE,
                                           upstream = NULL,
                                           downstream = NULL,
                                           label,
                                           ...) {

  listFlag <- FALSE
  if (is(tagMatrix[[1]][[1]],"matrix")) {
    if ( is.null(names(tagMatrix)) ) {
      nn <- paste0("peak", seq_along(tagMatrix))
      warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
      names(tagMatrix) <- nn
      ## stop("tagMatrix should be a named list...")
    }
    listFlag <- TRUE
  }

  if(listFlag){
    nbin <- dim(tagMatrix[[1]][[1]])[2]
    type <- attr(tagMatrix[[1]][[1]], 'type')
  }else{
    nbin <- dim(tagMatrix[[1]])[2]
    type <- attr(tagMatrix[[1]], 'type')
  }
  xlim <- c(1,nbin)

  if ( listFlag ) {
    facet <- match.arg(facet, c("none", "row", "column"))
  }

  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)
  ##
  conf <- if(missingArg(conf)) NA else conf

  pos <- value <- .id <- Lower <- Upper <- NULL

  if ( listFlag ) {

    tagCount <- lapply(as.list(names(tagMatrix)), function(x){

      tmp <- tagMatrix[[x]]
      tagCount_tmp <- lapply(as.list(names(tmp)),function(x){
        result <- getTagCount(tmp[[x]], xlim = xlim, conf = conf, ...)
        result$type <- x

        return(result)
      })
      tagCount_tmp <- list_to_dataframe(tagCount_tmp)
      return(tagCount_tmp)

    })

    names(tagCount) <- names(tagMatrix)
    tagCount <- list_to_dataframe(tagCount)
    tagCount$.id <- factor(tagCount$.id, levels=names(tagMatrix))
    p <- ggplot(tagCount, aes(pos, group=type, color=type))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = type),
                           linetype = 0, alpha = 0.2)
    }

  } else {

    tagCount <- lapply(as.list(names(tagMatrix)), function(x){

      result <- getTagCount(tagMatrix[[x]], xlim = xlim, conf = conf, ...)
      result$type <- x

      return(result)
    })

    tagCount <- do.call("rbind",tagCount)

    p <- ggplot(tagCount, aes(pos,group=type,color=type))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper,fill = type),
                           linetype = 0, alpha = 0.2)
    }
  }

  p <- p + geom_line(aes(y = value,color = type))

  ## x_scale for genebody
  if(type == 'body'){
    ## x_scale for gene body with no flank extension
    if(is.null(upstream)){
      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*0.25),
                                           floor(nbin*0.5),
                                           floor(nbin*0.75),
                                           nbin),
                                  labels=c(label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2]))
    }


    ## x_scale for flank extension by relative value
    if(inherits(upstream, 'rel')){

      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*(as.numeric(upstream)*100/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+25)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+50)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+75)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+100)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           nbin),
                                  labels=c(paste0("-",as.numeric(upstream)*100,"%"),
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0("+",as.numeric(downstream)*100,"%")))
      p <- p + geom_vline(xintercept=floor(nbin*(as.numeric(upstream)*100/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                          linetype="longdash")

      p <- p + geom_vline(xintercept=floor(nbin*((as.numeric(upstream)*100+100)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                          linetype="longdash")
    }

    ## x_scale for flank extension by absolute value
    if(!is.null(upstream) & !inherits(upstream, 'rel')){

      upstreamPer <- floor(upstream/1000)*0.1
      downstreamPer <- floor(downstream/1000)*0.1

      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.25)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.5)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.75)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))),
                                           nbin),
                                  labels=c(paste0("-",upstream,"bp"),
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0(downstream,"bp")))
      p <- p + geom_vline(xintercept=floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
                          linetype="longdash")

      p <- p + geom_vline(xintercept=floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))),
                          linetype="longdash")
    }
  }


  ## x_scale for start region
  if(type != 'body'){

    p <- p + scale_x_continuous(breaks=c(1,
                                         floor(nbin*0.25),
                                         floor(nbin*0.5),
                                         floor(nbin*0.75),
                                         nbin),
                                labels=c(paste0("-",upstream,"bp"),
                                         paste0("-",floor(upstream*0.5),"bp"),
                                         label,
                                         paste0(floor(downstream*0.5),"bp"),
                                         paste0(downstream,"bp")))

    p <- p + geom_vline(xintercept=floor(nbin*0.5),
                        linetype="longdash")
  }


  if (listFlag) {

    if (facet == "row") {
      if (free_y) {
        p <- p + facet_grid(.id ~ ., scales = "free_y")
      } else {
        p <- p + facet_grid(.id ~ .)
      }
    } else if (facet == "column") {
      if (free_y) {
        p <-  p + facet_grid(. ~ .id, scales = "free_y")
      } else {
        p <-  p + facet_grid(. ~ .id)
      }
    }
  }
  p <- p+xlab(xlab)+ylab(ylab)
  p <- p + theme_bw() + theme(legend.title=element_blank())
  # if(facet != "none") {
  #   p <- p + theme(legend.position="none")
  # }
  return(p)
}


##' plot the heatmap of tagMatrix
##'
##'
##' @title tagHeatmap
##' @param tagMatrix matrix or list of matrices. Tag matrix with features as rows
##'   and genomic positions as columns, typically obtained from \code{getTagMatrix()}.
##'   If a list is provided, creates separate heatmaps for each matrix
##' @param xlab character, x-axis label. Default is "" (empty)
##' @param ylab character, y-axis label. Default is "" (empty)
##' @param title character, plot title. If NULL, no title is shown. Default is NULL
##' @param palette character, color palette name from RColorBrewer for the heatmap.
##'   See \code{\link[ggplot2]{scale_colour_brewer}} for options. Default is "RdBu"
##' @param nrow integer, number of rows for arranging multiple heatmaps when
##'   \code{tagMatrix} is a list. If NULL, automatically determined. Default is NULL
##' @param ncol integer, number of columns for arranging multiple heatmaps when
##'   \code{tagMatrix} is a list. If NULL, automatically determined. Default is NULL
##' @return figure
##' @export
##' @author G Yu
tagHeatmap <- function(tagMatrix,
                       xlab="",
                       ylab="",
                       title=NULL,
                       palette="RdBu",
                       nrow = NULL,
                       ncol = NULL) {
  listFlag <- FALSE
  if (is(tagMatrix, "list")) {
    listFlag <- TRUE
  }
  peakHeatmap.internal2(tagMatrix = tagMatrix,
                        listFlag = listFlag,
                        palette = palette,
                        xlab = xlab,
                        ylab = ylab,
                        title = title,
                        ncol = ncol,
                        nrow = nrow)
}

##' Plot heatmap of peaks around genomic features
##'
##' This function creates a heatmap showing the distribution of ChIP-seq peaks
##' around genomic features (e.g., TSS, gene body). Each row represents a feature
##' and each column represents a genomic position.
##'
##' @description
##' The function generates tag matrices from peaks and visualizes them as heatmaps.
##' It supports both absolute positioning and binning modes. Features are typically
##' sorted by signal strength, allowing identification of features with strong
##' peak enrichment.
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Generates tag matrices from peaks using \code{getTagMatrix()}
##'   \item Optionally uses custom windows (if provided) or generates them from
##'         TxDb based on \code{by} and \code{type}
##'   \item Creates heatmap visualization with color-coded signal intensities
##'   \item Supports single or multiple peak sets (list input)
##' }
##'
##' @param peak peak file (BED format) or GRanges object, or a named list of
##'   peak files/GRanges objects for comparing multiple datasets
##' @param weightCol character, name of metadata column to use as weights. If NULL,
##'   all peaks have equal weight. Default is NULL
##' @param TxDb TxDb or EnsDb annotation object, or a GRanges object containing
##'   custom genomic features
##' @param upstream numeric, upstream distance in base pairs. Default is 1000
##' @param downstream numeric, downstream distance in base pairs. Default is 1000
##' @param xlab character, x-axis label. Default is "" (empty)
##' @param ylab character, y-axis label. Default is "" (empty)
##' @param title character, plot title. If NULL, no title is shown. Default is NULL
##' @param palette character, color palette name from RColorBrewer. If NULL, uses
##'   default. See \code{\link[ggplot2]{scale_colour_brewer}} for options.
##'   Default is NULL
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param by character, feature type. With TxDb: one of "gene", "transcript",
##'   "exon", "intron", "3UTR", "5UTR", "UTR". Default is "gene"
##' @param type character, one of "start_site", "end_site", or "body". Default
##'   is "start_site"
##' @param nbin integer, number of bins for binning mode. If NULL, uses absolute
##'   positioning. Default is NULL
##' @param ignore_strand logical, whether to ignore strand information. Default is FALSE
##' @param windows GRanges object, custom genomic regions to use instead of
##'   generating from TxDb. If missing, regions are generated from TxDb.
##'   Default is missing
##' @param nrow integer, number of rows for arranging multiple heatmaps. Default is NULL
##' @param ncol integer, number of columns for arranging multiple heatmaps. Default is NULL
##' @return A ggplot2 heatmap object showing peak signals across features and
##'   positions. For multiple peak sets, returns a grid of heatmaps
##' @export
##' @author G Yu
peakHeatmap <- function(peak, weightCol=NULL, TxDb=NULL,
                        upstream=1000, downstream=1000,
                        xlab="", ylab="", title=NULL,
                        palette=NULL, verbose=TRUE,
                        by="gene", type="start_site",
                        nbin = NULL,ignore_strand = FALSE,
                        windows,ncol = NULL, nrow = NULL) {
  listFlag <- FALSE
  if ( is(peak, "list") ) {
    listFlag <- TRUE
    if (is.null(names(peak)))
      stop("peak should be a peak file or a name list of peak files...")
  }

  if (verbose) {
    cat(">> preparing promoter regions...\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }

  if (verbose) {
    cat(">> preparing tag matrix...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }

  if(missing(windows)){
    windows <- getBioRegion(TxDb=TxDb,
                            upstream=upstream,
                            downstream=downstream,
                            by=by,
                            type=type)
  }


  if (listFlag) {
    tagMatrix <- lapply(peak, getTagMatrix,
                        weightCol=weightCol,
                        windows = windows,
                        upstream=upstream,
                        downstream=downstream,
                        TxDb = TxDb,
                        nbin = nbin,
                        verbose = verbose,
                        ignore_strand= ignore_strand)

    names(tagMatrix) <- names(peak)

  } else {
    tagMatrix <- getTagMatrix(peak,
                              weightCol=weightCol,
                              windows = windows,
                              TxDb = TxDb,
                              upstream=upstream,
                              downstream=downstream,
                              nbin = nbin,
                              verbose = verbose,
                              ignore_strand= ignore_strand)
  }

  if (verbose) {
    cat(">> generating figure...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }

  xlim <- NULL

  p <- peakHeatmap.internal2(tagMatrix = tagMatrix,
                             listFlag = listFlag,
                             palette = palette,
                             xlab = xlab,
                             ylab = ylab,
                             title = title,
                             nrow = nrow,
                             ncol = ncol)

  if (verbose) {
    cat(">> done...\t\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  invisible(tagMatrix)
  p
}

##' @importFrom aplot plot_list
##' @noRd
peakHeatmap.internal2 <- function(tagMatrix,
                                  listFlag,
                                  palette,
                                  xlab,
                                  ylab,
                                  title,
                                  nrow,
                                  ncol) {
  if ( is.null(xlab) || is.na(xlab))
    xlab <- ""
  if ( is.null(ylab) || is.na(ylab))
    ylab <- ""

  if (listFlag) {
    nc <- length(tagMatrix)
    if ( is.null(palette) || is.na(palette) ) {
      palette <- getPalette(nc)
    } else if (length(palette) != nc) {
      palette <- rep(palette[1], nc)
    } else {
      palette <- palette
    }

    if (is.null(title) || is.na(title))
      title <- names(tagMatrix)
    if (length(xlab) != nc) {
      xlab <- rep(xlab[1], nc)
    }
    if (length(ylab) != nc) {
      ylab <- rep(ylab[1], nc)
    }
    if (length(title) != nc) {
      title <- rep(title[1], nc)
    }

    tmp <- list()

    for (i in 1:nc) {

      p <- peakHeatmap.internal(tagMatrix = tagMatrix[[i]],
                                palette = palette[i],
                                xlab = xlab[i],
                                ylab = ylab[i],
                                title= title[i])

      p <- p + theme(plot.title = element_text(hjust = 0.5))

      tmp[[i]] <- p
    }

    if(is.null(nrow) && is.null(ncol))
      nrow <- 1

    p <- plot_list(gglist = tmp,
                   ncol = ncol,
                   nrow = nrow)
    return(p)

  } else {
    if (is.null(palette) || is.na(palette))
      palette <- "RdBu"
    if (is.null(title) || is.na(title))
      title <- ""
    peakHeatmap.internal(tagMatrix = tagMatrix,
                         palette = palette,
                         xlab = xlab,
                         ylab = ylab,
                         title = title)
  }
}


##' @import BiocGenerics
##' @importFrom yulab.utils mat2df
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 scale_fill_distiller
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 scale_x_continuous
##' @noRd
peakHeatmap.internal <- function(tagMatrix,
                                 palette="RdBu",
                                 xlab="",
                                 ylab="",
                                 title="") {

  upstream <- attr(tagMatrix, "upstream")
  downstream <- attr(tagMatrix, "downstream")
  binning_Flag <- attr(tagMatrix,"is.binning")
  type <- attr(tagMatrix,"type")

  body_Flag <- FALSE
  if(type == "body"){
    body_Flag <- TRUE
    label <- attr(tagMatrix,"label")
  }

  if(binning_Flag){
    nbin <- dim(tagMatrix)[2]
  }

  tagMatrix <- t(apply(tagMatrix, 1, function(x) x/max(x)))
  ii <- order(rowSums(tagMatrix))
  tagMatrix <- tagMatrix[ii,]

  colnames(tagMatrix) <- seq_len(dim(tagMatrix)[2])
  rownames(tagMatrix) <- seq_len(dim(tagMatrix)[1])

  tagMatrix <- mat2df(tagMatrix)
  colnames(tagMatrix) <- c("values","sample_ID","coordinate")

  sample_ID <- coordinate <- NULL

  p <- ggplot(tagMatrix, aes(x = coordinate,y = sample_ID)) +
    geom_tile(aes(fill = values)) +
    scale_fill_distiller(palette = palette)  +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          panel.grid=element_blank(),
          panel.background = element_blank()) +
    labs(x = xlab, y = ylab, title = title)

  if(body_Flag){

    if(inherits(upstream, 'rel')){

      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*(as.numeric(upstream)*100/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+25)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+50)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+75)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+100)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           nbin),
                                  labels=c(paste0("-",as.numeric(upstream)*100,"%"),
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0("+",as.numeric(downstream)*100,"%")))
    }

    if(is.null(upstream)){
      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*0.25),
                                           floor(nbin*0.5),
                                           floor(nbin*0.75),
                                           nbin),
                                  labels=c(label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2]))
    }

    if(!is.null(upstream) && !inherits(upstream, 'rel')){

      upstreamPer <- floor(upstream/1000)*0.1
      downstreamPer <- floor(downstream/1000)*0.1

      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.25)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.5)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.75)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))),
                                           nbin),
                                  labels=c(paste0("-",upstream,"bp"),
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0(downstream,"bp")))
    }

    p <- p + scale_y_continuous(expand = c(0,0))
    return(p)

  }

  if(binning_Flag){

    p <- p + scale_x_continuous(breaks = c(1,
                                           floor(nbin*(downstream*0.5/(downstream+upstream))),
                                           floor(nbin*(downstream/(downstream+upstream))),
                                           floor(nbin*((downstream + upstream*0.5)/(downstream+upstream))),
                                           nbin),
                                labels = c((-1*downstream),
                                           floor(-1*downstream*0.5),
                                           0,
                                           floor(upstream*0.5),
                                           upstream))
  }else{

    p <- p + scale_x_continuous(labels = function(x) x - upstream)

  }

  p <- p + scale_y_continuous(expand = c(0,0))

  p
}

##' Plot heatmap of peaks aligned to multiple sets of regions
##'
##' This function creates heatmaps showing peak distributions across multiple
##' sets of genomic regions (windows), allowing comparison of peak patterns
##' across different feature types or conditions.
##'
##' @description
##' The function extends \code{peakHeatmap()} to support multiple window sets.
##' It generates tag matrices for each window set and creates a faceted heatmap
##' visualization, with each facet showing peaks aligned to a different set of
##' regions. This is useful for comparing peak patterns across different feature
##' types (e.g., promoters vs enhancers) or conditions.
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Generates tag matrices for each window set (specified via \code{by}
##'         or custom \code{TxDb} GRanges)
##'   \item Creates heatmaps for each window set
##'   \item Arranges them in a grid with faceting
##'   \item Supports both absolute positioning and binning modes
##' }
##'
##' All windows must have the same \code{type} (start_site, end_site, or body)
##' and the same \code{upstream} and \code{downstream} parameters.
##'
##' @param peak peak file (BED format) or GRanges object, or a named list of
##'   peak files/GRanges objects for comparing multiple datasets
##' @param weightCol character, name of metadata column to use as weights. If NULL,
##'   all peaks have equal weight. Default is NULL
##' @param TxDb TxDb or EnsDb annotation object, or a named list of GRanges
##'   objects containing custom genomic features. If a list, names are used as
##'   window names
##' @param upstream numeric, upstream distance in base pairs. Must be the same
##'   for all windows. Default is 1000
##' @param downstream numeric, downstream distance in base pairs. Must be the
##'   same for all windows. Default is 1000
##' @param xlab character, x-axis label. Default is "" (empty)
##' @param ylab character, y-axis label. Default is "" (empty)
##' @param title character, plot title. If NULL, no title is shown. Default is NULL
##' @param palette character, color palette name from RColorBrewer. If NULL, uses
##'   default. Default is NULL
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param by character vector, feature types for each window. With TxDb: one
##'   of "gene", "transcript", "exon", "intron", "3UTR", "5UTR", "UTR". Length
##'   must match number of windows. Default is "gene"
##' @param type character, one of "start_site", "end_site", or "body". Must be
##'   the same for all windows. Default is "start_site"
##' @param nbin integer, number of bins for binning mode. If NULL, uses absolute
##'   positioning. Default is NULL
##' @param ignore_strand logical, whether to ignore strand information. Default is FALSE
##' @param windows_name character vector, names for each window (displayed as
##'   facet labels). If NULL and TxDb is a named list, uses names from TxDb.
##'   Default is NULL
##' @param nrow integer, number of rows for arranging heatmaps. Default is NULL (auto)
##' @param ncol integer, number of columns for arranging heatmaps. Default is NULL (auto)
##' @param facet_label_text_size numeric, size of facet label text. Default is 12
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 scale_fill_distiller
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 scale_x_continuous
##' @return figure
##' @export
peakHeatmap_multiple_Sets <- function(peak,
                                      weightCol=NULL,
                                      TxDb=NULL,
                                      upstream=1000,
                                      downstream=1000,
                                      xlab="",
                                      ylab="",
                                      title=NULL,
                                      palette=NULL,
                                      verbose=TRUE,
                                      by="gene",
                                      type="start_site",
                                      nbin = NULL,
                                      ignore_strand = FALSE,
                                      windows_name = NULL,
                                      ncol = NULL,
                                      nrow = NULL,
                                      facet_label_text_size = 12){
  listFlag <- FALSE
  if ( is(peak, "list") ) {
    listFlag <- TRUE
    if (is.null(names(peak)))
      stop("peak should be a peak file or a name list of peak files...")
  }

  if (verbose) {
    cat(">> preparing promoter regions...\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }


  ## check type
  if(length(type) != 1){
    stop("It should be only one type...")
  }

  if(is.null(windows_name) && !is.null(names(TxDb)))
    windows_name <- names(TxDb)

  ## make the window name
  if (is.null(windows_name)) {
    nn <- by
    warning("set the name automatically to ", paste(nn, collapse=' '))
    windows_name <- nn
  }else{
    if (length(windows_name) != length(by)) {
      stop("the length of the window name and the by should be equal...")
    }
  }

  if ( is(peak, "list") ) {
    tagMatrix <- lapply(peak, getTagMatrix2,
                        upstream=upstream,
                        downstream=downstream,
                        windows_name=windows_name,
                        type=type,
                        by=by,
                        TxDb=TxDb,
                        weightCol = weightCol,
                        nbin = nbin,
                        verbose = verbose,
                        ignore_strand= ignore_strand)
  } else {
    tagMatrix <- getTagMatrix2(peak=peak,
                               upstream=upstream,
                               downstream=downstream,
                               windows_name=windows_name,
                               type=type,
                               by=by,
                               TxDb=TxDb,
                               weightCol = weightCol,
                               nbin = nbin,
                               verbose = verbose,
                               ignore_strand= ignore_strand)
  }

  if(listFlag){

    nc <- length(tagMatrix)
    if ( is.null(palette) || is.na(palette) ) {
      palette <- getPalette(nc)
    } else if (length(palette) != nc) {
      palette <- rep(palette[1], nc)
    } else {
      palette <- palette
    }

    if (is.null(title) || is.na(title))
      title <- names(tagMatrix)
    if (length(xlab) != nc) {
      xlab <- rep(xlab[1], nc)
    }
    if (length(ylab) != nc) {
      ylab <- rep(ylab[1], nc)
    }
    if (length(title) != nc) {
      title <- rep(title[1], nc)
    }

    tmp <- list()

    for (i in 1:nc) {

      p <- peakHeatmap_multiple_Sets.internal(tagMatrix = tagMatrix[[i]],
                                              upstream=upstream,
                                              downstream=downstream,
                                              xlab=xlab[[i]],
                                              ylab=ylab[[i]],
                                              title=title[[i]],
                                              palette=palette[[i]],
                                              ncol = ncol,
                                              nrow = nrow,
                                              facet_label_text_size = facet_label_text_size)

      p <- p + theme(plot.title = element_text(hjust = 0.5))

      tmp[[i]] <- p
    }

    if(is.null(nrow) && is.null(ncol))
      nrow <- 1

    p <- plot_list(gglist = tmp,
                   ncol = ncol,
                   nrow = nrow)

  }else{

    if (is.null(palette) || is.na(palette))
      palette <- "RdBu"
    if (is.null(title) || is.na(title))
      title <- ""

    p <- peakHeatmap_multiple_Sets.internal(tagMatrix = tagMatrix,
                                            upstream=upstream,
                                            downstream=downstream,
                                            xlab=xlab,
                                            ylab=ylab,
                                            title=title,
                                            palette=palette,
                                            ncol = ncol,
                                            nrow = nrow,
                                            facet_label_text_size = facet_label_text_size)

  }

  return(p)

}


##' @importFrom yulab.utils mat2df
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 scale_fill_distiller
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 facet_grid
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 element_blank
##' @noRd
peakHeatmap_multiple_Sets.internal <- function(tagMatrix,
                                               upstream=1000,
                                               downstream=1000,
                                               xlab="",
                                               ylab="",
                                               title=NULL,
                                               palette=NULL,
                                               ncol = NULL,
                                               nrow = NULL,
                                               facet_label_text_size = 12){

  binning_Flag <- attr(tagMatrix[[1]],"is.binning")
  if(binning_Flag) nbin <- dim(tagMatrix[[1]])[2]

  type <- attr(tagMatrix,"type")
  body_Flag <- FALSE
  if(attr(tagMatrix[[1]],"type") == "body"){
    body_Flag <- TRUE
    label <- attr(tagMatrix,"label")
  }

  name_of_list <- as.list(names(tagMatrix))

  peak_list <- lapply(name_of_list,function(x){

    tagMatrix[[x]] <- t(apply(tagMatrix[[x]], 1, function(x) x/max(x)))
    ii <- order(rowSums(tagMatrix[[x]]))
    tagMatrix[[x]] <- tagMatrix[[x]][ii,]

    colnames(tagMatrix[[x]]) <- seq_len(dim(tagMatrix[[x]])[2])
    rownames(tagMatrix[[x]]) <- seq_len(dim(tagMatrix[[x]])[1])

    tagMatrix[[x]] <- mat2df(tagMatrix[[x]])
    colnames(tagMatrix[[x]]) <- c("values","sample_ID","coordinate")

    tagMatrix[[x]]$sample <- x
    return(tagMatrix[[x]])
  })

  peak_df <- list_to_dataframe(peak_list)

  sample_ID <- coordinate <- NULL

  p <- ggplot(peak_df, aes(x = coordinate,y = sample_ID)) +
    geom_tile(aes(fill = values)) +
    scale_fill_distiller(palette = palette)  +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          panel.grid=element_blank(),
          panel.background = element_blank()) +
    labs(x = xlab, y = ylab, title = title)

  if(body_Flag){

    if(inherits(upstream, 'rel')){

      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*(as.numeric(upstream)*100/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+25)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+50)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+75)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+100)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           nbin),
                                  labels=c(paste0("-",as.numeric(upstream)*100,"%"),
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0("+",as.numeric(downstream)*100,"%")))
    }

    if(is.null(upstream)){
      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*0.25),
                                           floor(nbin*0.5),
                                           floor(nbin*0.75),
                                           nbin),
                                  labels=c(label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2]))
    }

    if(!is.null(upstream) && !inherits(upstream, 'rel')){

      upstreamPer <- floor(upstream/1000)*0.1
      downstreamPer <- floor(downstream/1000)*0.1

      p <- p + scale_x_continuous(breaks=c(1,
                                           floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.25)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.5)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.75)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))),
                                           nbin),
                                  labels=c(paste0("-",upstream,"bp"),
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0(downstream,"bp")))
    }

    p <-  p + facet_grid(sample ~ .,switch = "y",space = "free_y",scales = "free_y") +
      theme(strip.text.y.left = element_text(color = "black",face = "bold",
                                             size = facet_label_text_size),
            strip.background = element_blank())

    return(p)

  }

  if(binning_Flag){

    p <- p + scale_x_continuous(breaks = c(1,
                                           floor(nbin*(downstream*0.5/(downstream+upstream))),
                                           floor(nbin*(downstream/(downstream+upstream))),
                                           floor(nbin*((downstream + upstream*0.5)/(downstream+upstream))),
                                           nbin),
                                labels = c((-1*downstream),
                                           floor(-1*downstream*0.5),
                                           0,
                                           floor(upstream*0.5),
                                           upstream))
  }else{

    p <- p + scale_x_continuous(breaks = c(1,
                                           floor(downstream*0.5),
                                           (downstream + 1),
                                           (downstream + 1 + floor(upstream * 0.5)),
                                           upstream+downstream+1),
                                labels = c((-1*downstream),
                                           floor(-1*downstream*0.5),
                                           0,
                                           floor(upstream*0.5),
                                           upstream))

  }

  p <-  p + facet_grid(sample ~ .,switch = "y",scales = "free_y",space = "free") +
    theme(strip.text.y.left = element_text(color = "black",face = "bold",
                                           size = facet_label_text_size),
          strip.background = element_blank()) +
    scale_y_continuous(expand = c(0,0))

  return(p)

}




##' Plot combined peak heatmap and average profile
##'
##' This function creates a combined visualization showing both a heatmap of
##' individual features and an average profile line plot in a single figure,
##' providing both detailed and summary views of peak distributions.
##'
##' @description
##' The function generates a two-panel figure:
##' \itemize{
##'   \item Top panel: Heatmap showing signal intensity for each feature (row)
##'         across genomic positions (columns)
##'   \item Bottom panel: Average profile line showing the mean signal across
##'         all features at each position
##' }
##'
##' This combination allows users to see both individual feature patterns and
##' the overall trend simultaneously.
##'
##' @param peak peak file (BED format) or GRanges object, or a named list of
##'   peak files/GRanges objects for comparing multiple datasets
##' @param weightCol character, name of metadata column to use as weights. If NULL,
##'   all peaks have equal weight. Default is NULL
##' @param TxDb TxDb or EnsDb annotation object, or a GRanges object containing
##'   custom genomic features
##' @param upstream numeric, upstream distance in base pairs. Default is 1000
##' @param downstream numeric, downstream distance in base pairs. Default is 1000
##' @param xlab character, x-axis label. Default is "" (empty)
##' @param ylab character, y-axis label. Default is "" (empty)
##' @param title character, plot title. If NULL, no title is shown. Default is NULL
##' @param palette character, color palette name from RColorBrewer for the heatmap.
##'   See \code{\link[ggplot2]{scale_colour_brewer}} for options. Default is NULL
##' @param verbose logical, whether to print progress messages. Default is TRUE
##' @param by character, feature type. With TxDb: one of "gene", "transcript",
##'   "exon", "intron", "3UTR", "5UTR", "UTR". Default is "gene"
##' @param type character, one of "start_site", "end_site", or "body". Default
##'   is "start_site"
##' @param nbin integer, number of bins for binning mode. If NULL, uses absolute
##'   positioning. Default is NULL
##' @param ignore_strand logical, whether to ignore strand information. Default is FALSE
##' @param windows_name character vector, names for windows (for multi-window plots).
##'   Default is NULL
##' @param nrow integer, number of rows for arranging multiple plots. Default is NULL
##' @param ncol integer, number of columns for arranging multiple plots. Default is NULL
##' @param facet_label_text_size numeric, size of facet label text. Default is 12
##' @param conf numeric, confidence level for confidence intervals in the profile
##'   plot (e.g., 0.95). If NA or missing, no confidence intervals are shown
##' @param facet character, one of "none", "row", or "column" for faceting the
##'   profile plot. Default is "none"
##' @param free_y logical, if TRUE, y-axis scales independently per facet in the
##'   profile plot. Default is TRUE
##' @param height_proportion numeric, proportion of height allocated to profile
##'   plot vs heatmap. Values between 0 and 1, where 0.5 means equal heights.
##'   Default is 0.5
##' @importFrom aplot insert_bottom
##' @importFrom aplot plot_list
##' @export
peak_Profile_Heatmap <- function(peak,
                                 weightCol=NULL,
                                 TxDb=NULL,
                                 upstream=1000,
                                 downstream=1000,
                                 xlab="",
                                 ylab="",
                                 title=NULL,
                                 palette=NULL,
                                 verbose=TRUE,
                                 by="gene",
                                 type="start_site",
                                 nbin = NULL,
                                 ignore_strand = FALSE,
                                 windows_name = NULL,
                                 ncol = NULL,
                                 nrow = NULL,
                                 facet_label_text_size = 12,
                                 conf,
                                 facet = "row",
                                 free_y = TRUE,
                                 height_proportion = 4){

  conf <- if(missingArg(conf)) NA else conf

  if(is(peak, "list")){

    nc <- length(peak)

    tmp <- list()

    if ( is.null(names(peak)) ) {
      nn <- paste0("peak", seq_along(peak))
      warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
      names(peak) <- nn
      ## stop("tagMatrix should be a named list...")
    }

    if(is.null(palette)) palette <- getPalette(nc)

    if(is.null(title)) title_of_plot <- names(peak)

    for (i in 1:nc) {
      peak_profile <- plotPeakProf(peak = peak[[i]],
                                   upstream = upstream,
                                   downstream = downstream,
                                   conf = conf,
                                   by = by,
                                   type = type,
                                   windows_name = windows_name,
                                   weightCol = weightCol,
                                   TxDb = TxDb,
                                   xlab = xlab,
                                   ylab = ylab,
                                   facet = facet,
                                   free_y = free_y,
                                   verbose = verbose,
                                   nbin = nbin,
                                   ignore_strand = ignore_strand)

      peak_profile <- peak_profile + labs(title = title_of_plot[i]) +
        theme(plot.title = element_text(hjust = 0.5))

      if(length(by) != 1){
        peak_heatmap <- peakHeatmap_multiple_Sets(peak = peak[[i]],
                                                  weightCol=weightCol,
                                                  TxDb=TxDb,
                                                  upstream=upstream,
                                                  downstream=downstream,
                                                  xlab=xlab,
                                                  ylab=ylab,
                                                  title=title,
                                                  palette=palette[[i]],
                                                  verbose=verbose,
                                                  by=by,
                                                  type=type,
                                                  nbin = nbin,
                                                  ignore_strand = ignore_strand,
                                                  windows_name = windows_name,
                                                  ncol = ncol,
                                                  nrow = nrow,
                                                  facet_label_text_size = facet_label_text_size)
      }else{

        peak_heatmap <- peakHeatmap(peak[[i]],
                                    weightCol=weightCol,
                                    TxDb=TxDb,
                                    upstream=upstream,
                                    downstream=downstream,
                                    xlab=xlab,
                                    ylab=ylab,
                                    title=title,
                                    palette=palette[[i]],
                                    verbose=verbose,
                                    by=by,
                                    type=type,
                                    nbin = nbin,
                                    ignore_strand = ignore_strand,
                                    ncol = ncol,
                                    nrow = nrow)

      }

      p <- peak_profile %>%
        insert_bottom(peak_heatmap,height = height_proportion)

      tmp[[i]] <- p
    }

    if (is.null(ncol) && is.null(nrow))
      nrow <- 1

    p <- plot_list(gglist = tmp,
                   ncol = ncol,
                   nrow = nrow)

    return(p)

  }

  peak_profile <- plotPeakProf(peak = peak,
                               upstream = upstream,
                               downstream = downstream,
                               conf = conf,
                               by = by,
                               type = type,
                               windows_name = windows_name,
                               weightCol = weightCol,
                               TxDb = TxDb,
                               xlab = xlab,
                               ylab = ylab,
                               facet = facet,
                               free_y = free_y,
                               verbose = verbose,
                               nbin = nbin,
                               ignore_strand = ignore_strand)


  if(length(by) != 1){
    peak_heatmap <- peakHeatmap_multiple_Sets(peak = peak,
                                              weightCol=weightCol,
                                              TxDb=TxDb,
                                              upstream=upstream,
                                              downstream=downstream,
                                              xlab=xlab,
                                              ylab=ylab,
                                              title=title,
                                              palette=palette,
                                              verbose=verbose,
                                              by=by,
                                              type=type,
                                              nbin = nbin,
                                              ignore_strand = ignore_strand,
                                              windows_name = windows_name,
                                              ncol = ncol,
                                              nrow = nrow,
                                              facet_label_text_size = facet_label_text_size)
  }else{

    peak_heatmap <- peakHeatmap(peak = peak,
                                weightCol=weightCol,
                                TxDb=TxDb,
                                upstream=upstream,
                                downstream=downstream,
                                xlab=xlab,
                                ylab=ylab,
                                title=title,
                                palette=palette,
                                verbose=verbose,
                                by=by,
                                type=type,
                                nbin = nbin,
                                ignore_strand = ignore_strand,
                                ncol = ncol,
                                nrow = nrow)

  }

  p <- peak_profile %>%
    insert_bottom(peak_heatmap,height = height_proportion)


  return(p)
}