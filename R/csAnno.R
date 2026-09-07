##' Class "csAnno"
##'
##' This class represents the output of ChIPseeker peak annotation. It contains
##' annotated peaks along with their genomic annotations, nearest gene information,
##' and annotation statistics.
##'
##' @description
##' The \code{csAnno} class is the main output object from \code{annotatePeak()}.
##' It stores annotated peaks as a GRanges object along with metadata about the
##' annotation process and summary statistics. This object can be used for
##' visualization, further analysis, and conversion to other formats.
##'
##' @details
##' The class contains the following slots:
##' \itemize{
##'   \item \code{anno}: GRanges object containing all annotated peaks with
##'         metadata columns (annotation, geneId, distanceToTSS, etc.)
##'   \item \code{tssRegion}: numeric vector of length 2 specifying the TSS
##'         region used for promoter annotation (e.g., c(-3000, 3000))
##'   \item \code{level}: character, annotation level - either "transcript" or
##'         "gene"
##'   \item \code{hasGenomicAnnotation}: logical, whether genomic feature
##'         annotations (Promoter, Exon, UTR, etc.) were assigned
##'   \item \code{detailGenomicAnnotation}: data.frame with logical columns
##'         indicating which genomic features each peak overlaps
##'   \item \code{annoStat}: data.frame containing annotation statistics showing
##'         the percentage of peaks in each genomic feature category
##'   \item \code{peakNum}: numeric, total number of peaks that were annotated
##' }
##'
##' @name csAnno-class
##' @aliases csAnno-class show,csAnno-method vennpie,csAnno-method plotDistToTSS,csAnno-method plotAnnoBar,csAnno-method plotAnnoPie,csAnno-method upsetplot,csAnno-method subset,csAnno-method
##'
##' @docType class
##' @slot anno GRanges object containing annotated peaks with all metadata
##' @slot tssRegion numeric vector of length 2, TSS region for promoter annotation
##' @slot level character, annotation level ("transcript" or "gene")
##' @slot hasGenomicAnnotation logical, whether genomic annotations were assigned
##' @slot detailGenomicAnnotation data.frame, detailed genomic annotation matrix
##' @slot annoStat data.frame, annotation statistics summary
##' @slot peakNum numeric, total number of annotated peaks
##' @exportClass csAnno
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @seealso \code{\link{annotatePeak}}
##' @keywords classes
setClass("csAnno",
         representation=representation(
             anno = "GRanges",
             tssRegion = "numeric",
             level = "character",
             hasGenomicAnnotation = "logical",
             detailGenomicAnnotation="data.frame",
             annoStat="data.frame",
             peakNum="numeric"
             ))


##' Convert csAnno object to GRanges
##'
##' This function extracts the GRanges object containing annotated peaks from a
##' \code{csAnno} object.
##'
##' @description
##' The function returns the \code{anno} slot from the \code{csAnno} object,
##' which contains all annotated peaks as a GRanges object with metadata columns
##' (annotation, geneId, distanceToTSS, etc.).
##'
##' @param x \code{csAnno} object to convert
##' @return GRanges object containing all annotated peaks with their metadata
##'   columns. This is the same as accessing \code{x@anno} directly
##' @export
##' @examples
##' \dontrun{
##' peak.gr <- as.GRanges(peakAnno)
##' ## Equivalent to: peak.gr <- peakAnno@anno
##' }
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
as.GRanges <- function(x) {
    if (!is(x, "csAnno"))
        stop("not supported...")
    return(x@anno)
}

##' Get annotation statistics from csAnno object
##'
##' This function extracts the annotation statistics data frame from a
##' \code{csAnno} object.
##'
##' @description
##' The function returns the \code{annoStat} slot, which contains a summary of
##' how peaks are distributed across different genomic feature categories
##' (Promoter, Exon, Intron, UTR, Intergenic, etc.) with percentages.
##'
##' @param x \code{csAnno} object
##' @return A data.frame with columns:
##'   \itemize{
##'     \item \code{Feature}: genomic feature category (Promoter, 5' UTR, 3' UTR,
##'           Exon, Intron, Downstream, Intergenic, etc.)
##'     \item \code{Frequency}: percentage of peaks in each category
##'   }
##'   Only available if \code{hasGenomicAnnotation=TRUE}
##' @export
##' @examples
##' \dontrun{
##' stats <- getAnnoStat(peakAnno)
##' ## View annotation distribution
##' print(stats)
##' }
##' @seealso \code{\link{plotAnnoBar}} for visualizing these statistics
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
getAnnoStat <- function(x) {
    if (!is(x, "csAnno"))
        stop("not supported...")
    return(x@annoStat)
}



##' Combine multiple csAnno objects
##'
##' This function combines two or more \code{csAnno} objects into a single
##' object, merging their peaks and aggregating annotation statistics.
##'
##' @description
##' The function merges multiple \code{csAnno} objects that were created with
##' the same annotation parameters (TSS region, level, genomic annotation status).
##' It combines all peaks, merges annotation statistics, and creates a unified
##' summary. This is useful for comparing or combining results from multiple
##' ChIP-seq experiments.
##'
##' @details
##' The function performs the following operations:
##' \enumerate{
##'   \item Validates that all objects are \code{csAnno} instances
##'   \item Checks that all objects have the same \code{tssRegion}, \code{level},
##'         and \code{hasGenomicAnnotation} values
##'   \item Combines all GRanges objects from the \code{anno} slots
##'   \item Merges \code{detailGenomicAnnotation} data frames by row binding
##'   \item Aggregates annotation statistics by merging and summing frequencies
##'         across all objects
##'   \item Sums the total peak counts
##' }
##'
##' The resulting object contains all peaks from all input objects with combined
##' annotation statistics showing the overall distribution across genomic features.
##'
##' @param x \code{csAnno} object (first object to combine)
##' @param ... additional \code{csAnno} objects to combine. At least one additional
##'   object must be provided
##' @return A new \code{csAnno} object containing:
##'   \itemize{
##'     \item All peaks from all input objects combined
##'     \item Merged annotation statistics with aggregated frequencies
##'     \item Combined detail genomic annotations
##'     \item Sum of peak counts from all objects
##'   }
##'   The \code{tssRegion}, \code{level}, and \code{hasGenomicAnnotation} are
##'   inherited from the first object (all must be identical)
##' @export
##' @examples
##' \dontrun{
##' ## Combine two peak annotation results
##' combined <- combine_csAnno(peakAnno1, peakAnno2, peakAnno3)
##'
##' ## Visualize combined statistics
##' plotAnnoBar(combined)
##' }
##' @seealso \code{\link{annotatePeak}} for creating csAnno objects
##' @references https://github.com/YuLab-SMU/ChIPseeker/issues/157
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
combine_csAnno <- function(x, ...){
    z <- list(x, ...)

    if(sum(vapply(z, function(x) !is(x, "csAnno"), FUN.VALUE = logical(1))) != 0){
        stop("not supported...")
    }

    if(length(z)<2){
        stop("need two or more csAnno object...")
    }


    if(sum(!duplicated(lapply(z, function(x) x@tssRegion[1]))) != 1
       && sum(!duplicated(lapply(z, function(x) x@tssRegion[2]))) != 1){
        stop("the tss regions of different csAnno objects should be the same...")
    }

    if(sum(!duplicated(lapply(z, function(x) x@level))) != 1){
        stop("the level of different csAnno object should be the same...")
    }

    if(sum(!duplicated(lapply(z, function(x) x@hasGenomicAnnotation))) != 1){
        stop("the status of GenomicAnnotation should be the same...")
    }

    combine_tssRegion <- x@tssRegion
    combine_level <- x@level
    combine_hasGenomicAnnotation <- x@hasGenomicAnnotation

    combine_anno <- x@anno
    for(i in 2:length(z)){
        combine_anno <- c(combine_anno,z[[i]]@anno)
    }

    combine_detailGenomicAnnotation <- lapply(z, function(x) x@detailGenomicAnnotation)
    combine_detailGenomicAnnotation <- do.call("rbind",combine_detailGenomicAnnotation)

    combine_peakNum <- x@peakNum
    for(i in 2:length(z)){
        combine_peakNum <- combine_peakNum+z[[i]]@peakNum
    }

    feature <- x@annoStat$Feature
    for(i in 2:length(z)){
        if(length(feature)<length(z[[i]]@annoStat$Feature)){
            feature_levels <- levels(z[[i]]@annoStat$Feature)
            feature <- c(as.vector(feature),as.vector(z[[i]]@annoStat$Feature))
            feature <- feature[!duplicated(feature)]
            feature <- factor(feature,
                              levels = feature_levels)
            feature <- sort(feature)
        }else{
            feature_levels <- levels(feature)
            feature <- c(as.vector(feature),as.vector(z[[i]]@annoStat$Feature))
            feature <- feature[!duplicated(feature)]
            feature <- factor(feature,
                              levels = feature_levels)
            feature <- sort(feature)
        }
    }

    combine_annoStat <- data.frame(Feature=feature)

    for(i in 1:length(z)){
        combine_annoStat <- merge(combine_annoStat, z[[i]]@annoStat,
                                  by = "Feature", all = T, sort = F)
        combine_annoStat[is.na(combine_annoStat)] <- 0
        combine_annoStat <- combine_annoStat[order(combine_annoStat$Feature),]
    }

    total <- (ncol(combine_annoStat)-1)*100
    combine_annoStat$sum <- rowSums(combine_annoStat[, 2:ncol(combine_annoStat)])


    for (i in 1:length(combine_annoStat$sum)) {
        combine_annoStat$result[i] <- (combine_annoStat$sum[i]/total)*100
    }

    annoStat_result <- data.frame(Feature=combine_annoStat[,1],Frequency=combine_annoStat[,ncol(combine_annoStat)])

    res <- new("csAnno",
               anno = combine_anno,
               tssRegion = combine_tssRegion,
               level = combine_level,
               hasGenomicAnnotation = combine_hasGenomicAnnotation,
               detailGenomicAnnotation = combine_detailGenomicAnnotation,
               annoStat = annoStat_result,
               peakNum = combine_peakNum
    )

    return(res)
}

##' Venn pie plot method for csAnno objects
##'
##' This function creates a Venn pie plot showing the overlap between different
##' genomic annotation categories.
##'
##' @description
##' The function generates a specialized visualization that combines elements of
##' Venn diagrams and pie charts to show how peaks are distributed across genomic
##' features and how these features overlap. This is useful for understanding
##' the complexity of peak annotations.
##'
##' @param x \code{csAnno} object containing annotated peaks
##' @param r numeric, initial radius for the plot. Default is 0.2
##' @param cex numeric, character expansion factor for adjusting legend text size.
##'   Default is 1.2
##' @param ... additional parameters passed to the underlying plotting function
##' @return A plot object (typically displayed directly, not returned)
##' @usage vennpie(x, r = 0.2, cex=1.2, ...)
##' @exportMethod vennpie
##' @examples
##' \dontrun{
##' vennpie(peakAnno)
##' vennpie(peakAnno, r=0.3, cex=1.5)
##' }
##' @seealso \code{\link{plotAnnoPie}} for standard pie charts,
##'   \code{\link{upsetplot}} for upset plots
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @noRd
setMethod("vennpie", signature(x="csAnno"),
          function(x,
                   r = 0.2,
                   cex = 1.2,
                   ...) {
            vennpie.csAnno(x, r, cex, ...)
          }
          )


##' Upset plot method for csAnno objects
##'
##' This function creates an upset plot (also called UpSet plot) showing the
##' intersections between different genomic annotation categories.
##'
##' @description
##' The function generates an upset plot that visualizes how peaks overlap across
##' different genomic feature categories (Promoter, Exon, Intron, UTR, etc.).
##' This is particularly useful for understanding complex annotation patterns
##' where peaks may overlap multiple features simultaneously.
##'
##' @details
##' Upset plots are an alternative to Venn diagrams that can handle more than
##' 3-4 sets effectively. They show:
##' \itemize{
##'   \item Set intersections (combinations of features) as bars
##'   \item Individual sets (features) as dots connected by lines
##'   \item The size of each intersection
##' }
##'
##' @param x \code{csAnno} object containing annotated peaks with genomic
##'   annotations (requires \code{hasGenomicAnnotation=TRUE})
##' @param ... additional parameters passed to \code{enrichplot::upsetplot()}
##' @return A plot object (typically displayed directly, not returned)
##' @usage upsetplot(x, ...)
##' @importFrom enrichplot upsetplot
##' @exportMethod upsetplot
##' @examples
##' \dontrun{
##' upsetplot(peakAnno)
##' }
##' @seealso \code{\link{vennpie}} for Venn pie plots,
##'   \code{\link{plotAnnoBar}} for bar plots
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @noRd
setMethod("upsetplot", signature(x="csAnno"),
          function(x, ...) {
              upsetplot.csAnno(x, ...)
          }
          )

##' Convert csAnno object to data.frame
##'
##' This function converts a \code{csAnno} object to a data.frame by extracting
##' and converting the GRanges object containing annotated peaks.
##'
##' @description
##' The function converts the \code{anno} slot (GRanges object) to a data.frame,
##' preserving all metadata columns. This is useful for exporting data,
##' performing data.frame operations, or using with functions that require
##' data.frame input.
##'
##' @param x \code{csAnno} object to convert
##' @param row.names character vector or NULL, row names for the resulting
##'   data.frame. If NULL, default row names are used. Default is NULL
##' @param optional logical, kept for compatibility with generic method but
##'   should be omitted. Default is FALSE
##' @param ... additional parameters (currently unused)
##' @return A data.frame containing all peak information with columns:
##'   \itemize{
##'     \item Standard GRanges columns: seqnames, start, end, width, strand
##'     \item All metadata columns: annotation, geneId, distanceToTSS, geneChr,
##'           geneStart, geneEnd, etc.
##'   }
##' @method as.data.frame csAnno
##' @export
##' @examples
##' \dontrun{
##' peak.df <- as.data.frame(peakAnno)
##' ## Export to CSV
##' write.csv(peak.df, "annotated_peaks.csv")
##' }
##' @seealso \code{\link{as.GRanges}} for converting to GRanges object
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
as.data.frame.csAnno <- function(x, row.names=NULL, optional=FALSE, ...) {
    y <- as.GRanges(x)
    if (!(is.null(row.names) || is.character(row.names)))
        stop("'row.names' must be NULL or a character vector")
    df <- as.data.frame(y)
    rownames(df) <- row.names
    return(df)
}

##' Show method for csAnno objects
##'
##' This function displays a summary of a \code{csAnno} object when it is printed
##' or shown in the console.
##'
##' @description
##' The function prints a concise summary including the number of annotated peaks,
##' and if genomic annotations are available, displays the annotation statistics
##' showing the distribution of peaks across genomic features.
##'
##' @param object \code{csAnno} object to display
##' @return No return value, called for side effects (printing to console)
##' @importFrom methods show
##' @exportMethod show
##' @usage show(object)
##' @examples
##' \dontrun{
##' peakAnno  # Automatically calls show()
##' show(peakAnno)  # Explicit call
##' }
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @noRd
setMethod("show", signature(object="csAnno"),
          function(object) {
              cat("Annotated peaks generated by ChIPseeker\n")
              cat(paste(length(object@anno), object@peakNum, sep="/"),
                  " peaks were annotated\n")
              if (object@hasGenomicAnnotation) {
                  cat("Genomic Annotation Summary:\n")
                  print(object@annoStat)
              }
          }
          )

##' Plot annotation bar chart for list of csAnno objects
##'
##' This function creates a bar plot comparing annotation distributions across
##' multiple \code{csAnno} objects.
##'
##' @description
##' The function generates a grouped bar chart showing the percentage of peaks
##' in each genomic feature category for multiple ChIP-seq experiments, allowing
##' easy comparison of annotation patterns across datasets.
##'
##' @param x list of \code{csAnno} objects to compare. If unnamed, automatic
##'   names will be assigned
##' @param xlab character, label for x-axis. Default is "" (empty)
##' @param ylab character, label for y-axis. Default is "Percentage(\%)"
##' @param title character, plot title. Default is "Feature Distribution"
##' @param ... additional parameters passed to the underlying plotting function
##' @return A ggplot2 bar plot object showing annotation distributions for all
##'   datasets side by side
##' @name plotAnnoBar
##' @docType methods
##' @rdname plotAnnoBar-methods
##' @aliases plotAnnoBar,list-method
##' @exportMethod plotAnnoBar
##' @examples
##' \dontrun{
##' ## Compare multiple experiments
##' peakList <- list(experiment1=peakAnno1, experiment2=peakAnno2)
##' plotAnnoBar(peakList)
##' }
##' @seealso \code{\link{plotAnnoBar,csAnno-method}} for single object plotting
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @noRd
setMethod("plotAnnoBar", signature(x="list"),
          function(x,
                   xlab="",
                   ylab='Percentage(%)',
                   title="Feature Distribution",
                   ...) {
              if (is.null(names(x))) {
                  nn <- paste0("Peak", seq_along(x))
                  warning("input is not a named list, set the name automatically to ", paste(nn, collapse = " "))
                  names(x) <- nn
                  ## stop("input object should be a named list...")
              }
              anno <- lapply(x, getAnnoStat)
              ## anno.df <- ldply(anno)
              anno.df <- list_to_dataframe(anno)
              categoryColumn <- ".id"
              plotAnnoBar.data.frame(anno.df, xlab, ylab, title, categoryColumn)
          })

##' Plot annotation bar chart for csAnno object
##'
##' This function creates a bar plot showing the distribution of peaks across
##' different genomic feature categories.
##'
##' @description
##' The function generates a bar chart displaying the percentage of peaks in
##' each genomic annotation category (Promoter, 5' UTR, 3' UTR, Exon, Intron,
##' Downstream, Intergenic, etc.). This provides a quick visual summary of
##' where peaks are located in the genome.
##'
##' @param x \code{csAnno} object containing annotated peaks with genomic
##'   annotations (requires \code{hasGenomicAnnotation=TRUE})
##' @param xlab character, label for x-axis. Default is "" (empty)
##' @param ylab character, label for y-axis. Default is "Percentage(\%)"
##' @param title character, plot title. Default is "Feature Distribution"
##' @param ... additional parameters passed to the underlying plotting function
##' @return A ggplot2 bar plot object showing the percentage distribution of
##'   peaks across genomic features
##' @exportMethod plotAnnoBar
##' @usage plotAnnoBar(x, xlab="", ylab='Percentage(\%)',title="Feature Distribution", ...)
##' @examples
##' \dontrun{
##' plotAnnoBar(peakAnno)
##' plotAnnoBar(peakAnno, title="My ChIP-seq Experiment")
##' }
##' @seealso \code{\link{plotAnnoPie}} for pie chart visualization,
##'   \code{\link{getAnnoStat}} to get the underlying statistics
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @noRd
setMethod("plotAnnoBar", signature(x="csAnno"),
          function(x,
                   xlab="",
                   ylab="Percentage(%)",
                   title="Feature Distribution",
                   ...) {
              anno.df <- getAnnoStat(x)
              categoryColumn <- 1
              plotAnnoBar.data.frame(anno.df, xlab, ylab, title, categoryColumn)
          })



##' Plot annotation pie chart for csAnno object
##'
##' This function creates a pie chart showing the distribution of peaks across
##' different genomic feature categories.
##'
##' @description
##' The function generates a pie chart (2D or 3D) displaying the percentage of
##' peaks in each genomic annotation category. Each slice represents a different
##' genomic feature (Promoter, Exon, Intron, UTR, etc.), with the size
##' proportional to the percentage of peaks in that category.
##'
##' @param x \code{csAnno} object containing annotated peaks with genomic
##'   annotations (requires \code{hasGenomicAnnotation=TRUE})
##' @param ndigit integer, number of decimal places to display in percentage
##'   labels. Default is 2
##' @param cex numeric, character expansion factor for adjusting label text size.
##'   Default is 0.9
##' @param col character vector or NA, colors for pie slices. If NA, default
##'   colors are used. Default is NA
##' @param legend.position character, position of the legend. Options include
##'   "rightside", "topright", "bottomright", etc. Default is "rightside"
##' @param pie3D logical, whether to create a 3D pie chart. If FALSE, creates
##'   a standard 2D pie chart. Default is FALSE
##' @param radius numeric, radius of the pie chart. Values between 0 and 1.
##'   Default is 0.8
##' @param ... additional parameters passed to the underlying plotting function
##' @return A plot object (typically displayed directly, not returned)
##' @exportMethod plotAnnoPie
##' @usage plotAnnoPie(x,ndigit=2,cex=0.9,col=NA,legend.position="rightside",pie3D=FALSE,radius=0.8,...)
##' @examples
##' \dontrun{
##' plotAnnoPie(peakAnno)
##' plotAnnoPie(peakAnno, pie3D=TRUE, radius=0.9)
##' plotAnnoPie(peakAnno, col=c("red", "blue", "green"), legend.position="topright")
##' }
##' @seealso \code{\link{plotAnnoBar}} for bar chart visualization
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @noRd
setMethod("plotAnnoPie", signature(x="csAnno"),
          function(x,
                   ndigit=2,
                   cex=0.9,
                   col=NA,
                   legend.position="rightside",
                   pie3D=FALSE,
                   radius=0.8,
                   ...){
              plotAnnoPie.csAnno(x, ndigit, cex, col, legend.position, pie3D, radius, ...)
          })



##' Plot distance to TSS for list of csAnno objects
##'
##' This function creates a diverging bar plot comparing the distribution of
##' peak distances to TSS across multiple \code{csAnno} objects.
##'
##' @description
##' The function generates a grouped diverging bar chart showing the percentage
##' of peaks at different distances from TSS for multiple ChIP-seq experiments,
##' allowing comparison of binding site distributions across datasets.
##'
##' @param x list of \code{csAnno} objects to compare. If unnamed, automatic
##'   names will be assigned
##' @param distanceColumn character, name of the column containing distances
##'   to TSS. Default is "distanceToTSS"
##' @param xlab character, label for x-axis. Default is "" (empty)
##' @param ylab character, label for y-axis. Default is "Binding sites (\%) (5'->3')"
##' @param title character, plot title. Default is "Distribution of transcription
##'   factor-binding loci relative to TSS"
##' @param distanceBreaks numeric vector, breakpoints for distance categories
##'   in base pairs. Default is c(0, 1000, 3000, 5000, 10000, 100000)
##' @param palette character, color palette name from RColorBrewer or NULL for
##'   default colors. Default is NULL
##' @param ... additional parameters passed to the underlying plotting function
##' @return A ggplot2 diverging bar plot object showing distance distributions
##'   for all datasets
##' @name plotDistToTSS
##' @docType methods
##' @rdname plotDistToTSS-methods
##' @aliases plotDistToTSS,list-method
##' @exportMethod plotDistToTSS
##' @examples
##' \dontrun{
##' ## Compare multiple experiments
##' peakList <- list(exp1=peakAnno1, exp2=peakAnno2)
##' plotDistToTSS(peakList)
##' }
##' @seealso \code{\link{plotDistToTSS,csAnno-method}} for single object plotting
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @noRd
setMethod("plotDistToTSS", signature(x="list"),
          function(x, distanceColumn="distanceToTSS",
                                     xlab="", ylab="Binding sites (%) (5'->3')",
                                     title="Distribution of transcription factor-binding loci relative to TSS",
                     distanceBreaks=c(0, 1000, 3000, 5000, 10000, 100000),
                     palette = NULL, ...) {
              if (is.null(names(x))) {
                  nn <- paste0("Peak", seq_along(x))
                  warning("input is not a named list, set the name automatically to ", paste(nn, collapse = " "))
                  names(x) <- nn
                  ## stop("input object should be a named list...")
              }

              peakAnno <- lapply(x, as.data.frame)
              ## peakDist <- ldply(peakAnno)
              peakDist <- list_to_dataframe(peakAnno)
              categoryColumn <- ".id"
              plotDistToTSS.data.frame(peakDist, distanceColumn = distanceColumn,
                                       distanceBreaks = distanceBreaks, palette = palette,
                                       xlab = xlab, ylab = ylab, title = title, categoryColumn = categoryColumn)
          })


##' Plot distance to TSS for csAnno object
##'
##' This function creates a diverging bar plot showing the distribution of peaks
##' relative to transcription start sites (TSS).
##'
##' @description
##' The function generates a diverging bar chart that visualizes where peaks are
##' located relative to gene TSSs. Peaks are categorized into distance bins
##' (e.g., 0-1kb, 1-3kb, etc.), and the plot shows:
##' \itemize{
##'   \item Upstream peaks (5' of TSS) on the left side (negative bars)
##'   \item Downstream peaks (3' of TSS) on the right side (positive bars)
##'   \item TSS marked at the center (0)
##'   \item Distance categories color-coded
##' }
##'
##' This visualization helps understand whether binding sites are enriched near
##' promoters, evenly distributed, or show other patterns.
##'
##' @param x \code{csAnno} object containing annotated peaks with distanceToTSS
##'   information
##' @param distanceColumn character, name of the column containing distances
##'   to TSS. Default is "distanceToTSS"
##' @param xlab character, label for x-axis. Default is "" (empty)
##' @param ylab character, label for y-axis. Default is "Binding sites (\%) (5'->3')"
##' @param title character, plot title. Default is "Distribution of transcription
##'   factor-binding loci relative to TSS"
##' @param distanceBreaks numeric vector, breakpoints for distance categories
##'   in base pairs. The function automatically adds 0 and Inf if not present.
##'   Default is c(0, 1000, 3000, 5000, 10000, 100000)
##' @param palette character, color palette name from RColorBrewer (e.g., "Set1",
##'   "Dark2") or NULL for default colors. Run
##'   \code{RColorBrewer::display.brewer.all()} to see available palettes.
##'   Default is NULL
##' @param ... additional parameters passed to the underlying plotting function
##' @return A ggplot2 diverging bar plot object showing the distribution of
##'   peaks relative to TSS
##' @exportMethod plotDistToTSS
##' @usage plotDistToTSS(x,distanceColumn="distanceToTSS", xlab="",
##' ylab="Binding sites (\%) (5'->3')",
##' title="Distribution of transcription factor-binding loci relative to TSS",...)
##' @examples
##' \dontrun{
##' plotDistToTSS(peakAnno)
##' plotDistToTSS(peakAnno, distanceBreaks=c(0, 2000, 5000, 10000, 50000))
##' plotDistToTSS(peakAnno, palette="Set2")
##' }
##' @seealso \code{\link{annotatePeak}} for creating csAnno objects with
##'   distanceToTSS information
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @noRd
setMethod("plotDistToTSS", signature(x="csAnno"),
          function(x, distanceColumn="distanceToTSS",
                                     xlab="", ylab="Binding sites (%) (5'->3')",
                                     title="Distribution of transcription factor-binding loci relative to TSS",
                                     distanceBreaks=c(0, 1000, 3000, 5000, 10000, 100000),
                                     palette = NULL,...) {
              peakDist <- as.data.frame(x)
              categoryColumn <- 1
              plotDistToTSS.data.frame(peakDist, distanceColumn = distanceColumn, distanceBreaks = distanceBreaks, palette = palette,
                                       xlab = xlab, ylab = ylab, title = title, categoryColumn = categoryColumn)
          })

