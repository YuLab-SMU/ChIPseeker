##' Create UpSet plot for genomic annotation combinations
##'
##' This function visualizes the combinations of genomic annotations (e.g.,
##' Promoter, Exon, Intron, UTR) that peaks can have, using an UpSet plot to
##' show set intersections and their frequencies.
##'
##' @description
##' The function creates an UpSet plot (an alternative to Venn diagrams) that
##' shows how peaks are distributed across different combinations of genomic
##' annotations. It displays:
##' \itemize{
##'   \item The frequency of each annotation combination (intersection)
##'   \item Which annotations are present in each combination
##'   \item Optionally, a vennpie plot overlay showing the overall annotation
##'         distribution
##' }
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Extracts the \code{detailGenomicAnnotation} data frame from the
##'         \code{csAnno} object, which contains logical columns indicating
##'         which annotations each peak has
##'   \item Converts the annotation matrix into a list format where each peak
##'         is represented by the set of annotations it has
##'   \item Creates an UpSet plot using \code{ggupset} showing:
##'     \itemize{
##'       \item Left panel: bar chart of intersection sizes (frequencies)
##'       \item Top panel: bar chart showing individual set sizes
##'       \item Bottom panel: matrix showing which sets are included in each
##'             intersection
##'     }
##'   \item Optionally overlays a vennpie plot if \code{vennpie=TRUE}
##' }
##'
##' UpSet plots are particularly useful when there are many annotation categories
##' and peaks can have multiple annotations simultaneously (e.g., a peak can
##' be both in a Promoter and overlap an Exon).
##'
##' @param x \code{csAnno} object containing annotated peaks with
##'   \code{detailGenomicAnnotation} slot
##' @param order_by character, how to order the intersections in the plot.
##'   Options: "freq" (by frequency, default) or "degree" (by number of sets
##'   in the intersection). Default is "freq"
##' @param vennpie logical, whether to overlay a vennpie plot showing the
##'   overall annotation distribution. If TRUE, the vennpie is embedded as
##'   a subplot. Default is FALSE
##' @param vp list, viewport parameters for positioning the vennpie subplot
##'   when \code{vennpie=TRUE}. Contains:
##'   \itemize{
##'     \item \code{x}: x position (0-1, left to right). Default is 0.6
##'     \item \code{y}: y position (0-1, bottom to top). Default is 0.7
##'     \item \code{width}: width of subplot (0-1). Default is 0.8
##'     \item \code{height}: height of subplot (0-1). Default is 0.8
##'   }
##' @return A ggplot2 object (or combined plot if vennpie=TRUE) showing the
##'   UpSet visualization. The plot can be further customized using ggplot2
##'   functions or printed directly
##' @importFrom ggplot2 coord_fixed
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_minimal
##' @seealso \code{\link{vennpie}} for the vennpie plot function,
##'   \code{\link{plotAnnoBar}} for bar plots of annotation distribution
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, TxDb=txdb)
##'
##' ## Basic UpSet plot
##' upsetplot(peakAnno)
##'
##' ## UpSet plot with vennpie overlay
##' upsetplot(peakAnno, vennpie=TRUE)
##'
##' ## Order by degree instead of frequency
##' upsetplot(peakAnno, order_by="degree")
##' }
##' @author Guangchuang Yu
upsetplot.csAnno <- function(x, order_by = "freq", vennpie=FALSE, vp = list(x=.6, y=.7, width=.8, height=.8)) {
    y <- x@detailGenomicAnnotation
    nn <- names(y)
    y <- as.matrix(y)

    res <- tibble::tibble(anno = lapply(1:nrow(y), function(i) nn[y[i,]]))
    g <- ggplot(res, aes_(x = ~anno)) + geom_bar() +
        xlab(NULL) + ylab(NULL) + theme_minimal() +
        ggupset::scale_x_upset(n_intersections = 20, order_by = order_by)

    if (!vennpie) return(g)

    f <- function() vennpie(x, cex = .9)

    p <- ggplotify::as.ggplot(f) + coord_fixed()

    ggplotify::as.ggplot(g) +
        ggimage::geom_subview(subview = p, x = vp$x, y = vp$y, width = vp$width, height = vp$height)


    ## y[y] <- 1
    ## y <- as.data.frame(y)
    ## ## cn <- colnames(y)
    ## ## cn[cn == "fiveUTR"] <- "5 UTR"
    ## ## cn[cn == "threeUTR"] <- "3 UTR"
    ## ## colnames(y) <- cn

    ## if (is.null(sets)) {
    ##     sets <- c("distal_intergenic", "downstream",
    ##               "threeUTR", "fiveUTR", "Intron",
    ##               "Exon", "Promoter")
    ##     if (vennpie && is.null(sets.bar.color)) {
    ##         sets.bar.color <- c("#d95f0e", "#fee0d2", "#98D277",
    ##                             "#6F9E4C", "#fc9272", "#9ecae1", "#ffeda0")
    ##     }
    ## }

    ## if (is.null(sets.bar.color)) {
    ##     sets.bar.color <- "black"
    ## }

    ## if (vennpie) {
    ##     plot.new()
    ##     # grid.rect(gp = gpar(fill="white"))
    ##     upset(y, sets=sets, sets.bar.color=sets.bar.color,
    ##           order.by = order.by, ...)
    ##     pushViewport(vp)
    ##     ##par(plt=gridPLT(), new=TRUE)
    ##     vennpie(x)
    ##     popViewport()
    ## } else {
    ##     upset(y, sets=sets,sets.bar.color=sets.bar.color,
    ##           order.by = order.by, ...)
    ## }
}
