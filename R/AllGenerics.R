##' Create Venn pie chart for genomic annotation distribution
##'
##' This generic function creates a nested pie chart (venn diagram style)
##' visualization showing the distribution of peaks across different genomic
##' annotation categories, with hierarchical nesting to show overlapping annotations.
##'
##' @description
##' The function visualizes the genomic annotation of peaks using a nested pie
##' chart approach. It displays three concentric pie charts showing:
##' \itemize{
##'   \item \strong{Outer ring}: Genic vs Intergenic distribution
##'   \item \strong{Middle ring}: Breakdown of Genic (Intron, Exon) and Intergenic
##'         (Upstream, Downstream, Distal Intergenic) categories
##'   \item \strong{Inner ring}: Further breakdown showing Exon and Downstream
##'         regions in detail
##' }
##'
##' This visualization helps understand the overlap between different annotation
##' categories, as peaks can be annotated with multiple features simultaneously.
##'
##' @param x \code{csAnno} object containing annotated peaks with
##'   \code{detailGenomicAnnotation} slot populated
##' @param r numeric, initial radius for the base pie chart. Controls the overall
##'   size of the plot. The nested pies use multiples of this radius (2*r, 3*r, 4*r).
##'   Default is 0.2
##' @param cex numeric, character expansion factor for the legend text. Larger
##'   values make the legend text bigger. Default is 1.2
##' @param ... additional parameters passed to the plotting functions
##' @return A plot showing nested pie charts with genomic annotation distribution
##' @docType methods
##' @name vennpie
##' @rdname vennpie-methods
##' @seealso \code{\link{plotAnnoBar}} \code{\link{plotAnnoPie}} for alternative
##'   visualization methods
##' @export
setGeneric("vennpie",
  function(x, r = 0.2, cex = 1.2, ...)
  standardGeneric("vennpie")
)


##' Plot feature distribution based on distances to TSS
##'
##' This generic function creates a bar plot showing the distribution of peaks
##' relative to transcription start sites (TSS), categorizing peaks by their
##' distance from the nearest gene's TSS.
##'
##' @description
##' The function visualizes where ChIP-seq peaks are located relative to gene
##' TSSs. It categorizes peaks into distance bins (e.g., 0-1kb, 1-3kb, etc.)
##' and displays the percentage of peaks in each category, separately for
##' upstream (5') and downstream (3') regions. The plot uses a diverging bar
##' chart with TSS at the center (0).
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Categorizes peaks into distance bins based on \code{distanceBreaks}
##'   \item Separates upstream (negative distance) and downstream (positive
##'         distance) peaks
##'   \item Calculates percentages for each distance category
##'   \item Creates a diverging bar plot with upstream peaks on the left and
##'         downstream peaks on the right, with TSS at the center
##' }
##'
##' @param x \code{csAnno} object or list of \code{csAnno} objects containing
##'   annotated peaks with distance information
##' @param distanceColumn character, name of the column containing distances
##'   from peaks to nearest gene TSS. Default is "distanceToTSS"
##' @param xlab character, label for x-axis. Default is "" (empty)
##' @param ylab character, label for y-axis. Default is "Binding sites (\%) (5'->3')"
##' @param title character, plot title. Default is "Distribution of transcription
##'   factor-binding loci relative to TSS"
##' @param ... additional parameters passed to the implementation methods, such as:
##'   \itemize{
##'     \item \code{distanceBreaks}: numeric vector of breakpoints for distance
##'           categories (default: c(0, 1000, 3000, 5000, 10000, 100000))
##'     \item \code{palette}: color palette name from RColorBrewer
##'     \item \code{categoryColumn}: column for grouping multiple datasets
##'   }
##' @return A ggplot2 bar plot object showing the distribution of peaks relative
##'   to TSS. The plot uses a diverging bar chart format with upstream (5') on
##'   the left and downstream (3') on the right, with TSS at the center
##' @docType methods
##' @name plotDistToTSS
##' @rdname plotDistToTSS-methods
##' @aliases plotDistToTSS,list-method
##' @seealso \code{\link{annotatePeak}} for peak annotation,
##'   \code{\link{plotAnnoBar}} for annotation distribution plots
##' @export
setGeneric("plotDistToTSS",
  function(x,
    distanceColumn="distanceToTSS",
    xlab="", ylab="Binding sites (%) (5'->3')",
    title="Distribution of transcription factor-binding loci relative to TSS",
    ...)
  standardGeneric("plotDistToTSS")
)

##' Plot bar chart of genomic annotation distribution
##'
##' This generic function creates a bar plot showing the distribution of peaks
##' across different genomic annotation categories (Promoter, Exon, Intron, etc.).
##'
##' @description
##' The function visualizes the percentage of peaks falling into each genomic
##' annotation category. It creates a horizontal bar chart where each bar
##' represents a category, and the bar length represents the percentage of peaks
##' in that category. The function supports both single \code{csAnno} objects
##' and lists of \code{csAnno} objects for comparing multiple datasets.
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Extracts annotation statistics from the \code{csAnno} object(s)
##'   \item Calculates percentages for each annotation category
##'   \item Creates a horizontal bar chart with categories on the y-axis and
##'         percentages on the x-axis
##'   \item Color-codes each annotation category
##' }
##'
##' Common annotation categories include: Promoter, 5' UTR, 3' UTR, Exon,
##' Intron, Downstream, and Intergenic.
##'
##' @param x \code{csAnno} object or list of \code{csAnno} objects containing
##'   annotated peaks with genomic annotation information
##' @param xlab character, label for x-axis. Default is "" (empty)
##' @param ylab character, label for y-axis. Default is "Percentage(\%)"
##' @param title character, plot title. Default is "Feature Distribution"
##' @param ... additional parameters passed to the implementation methods
##' @return A ggplot2 bar plot object showing the distribution of peaks across
##'   genomic annotation categories. For lists of \code{csAnno} objects, the plot
##'   shows side-by-side bars for comparison
##' @docType methods
##' @name plotAnnoBar
##' @rdname plotAnnoBar-methods
##' @aliases plotAnnoBar,list-method
##' @seealso \code{\link{plotAnnoPie}} for pie chart visualization,
##'   \code{\link{annotatePeak}} for peak annotation
##' @export
setGeneric("plotAnnoBar",
  function(x,
    xlab="",
    ylab="Percentage(%)",
    title="Feature Distribution",
    ...)
  standardGeneric("plotAnnoBar")
)


##' Plot pie chart of genomic annotation distribution
##'
##' This generic function creates a pie chart showing the distribution of peaks
##' across different genomic annotation categories.
##'
##' @description
##' The function visualizes the percentage of peaks falling into each genomic
##' annotation category using a pie chart. Each slice of the pie represents an
##' annotation category, with the slice size proportional to the percentage of
##' peaks in that category. Labels show the category name and percentage.
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Extracts annotation statistics from the \code{csAnno} object
##'   \item Calculates percentages for each annotation category
##'   \item Creates a pie chart with color-coded slices
##'   \item Adds labels showing category names and percentages (rounded to
##'         \code{ndigit} decimal places)
##'   \item Optionally creates a 3D pie chart for enhanced visualization
##' }
##'
##' Common annotation categories include: Promoter, 5' UTR, 3' UTR, Exon,
##' Intron, Downstream, and Intergenic.
##'
##' @param x \code{csAnno} object containing annotated peaks with genomic
##'   annotation information
##' @param ndigit integer, number of decimal places to round percentages in
##'   labels. Default is 2
##' @param cex numeric, character expansion factor for labels. Larger values
##'   make labels bigger. Default is 0.9
##' @param col character vector, colors for pie slices. If NA (default),
##'   automatic colors are used
##' @param legend.position character, position of the legend. Options include
##'   "rightside", "leftside", "top", "bottom", or specific coordinates.
##'   Default is "rightside"
##' @param pie3D logical, whether to create a 3D pie chart. If TRUE, uses
##'   \code{pie3D()} from the \code{plotrix} package. Default is FALSE
##' @param radius numeric, radius of the pie chart. Values between 0 and 1.
##'   Default is 0.8
##' @param ... additional parameters passed to the pie chart plotting functions
##' @return A pie chart plot showing the distribution of peaks across genomic
##'   annotation categories. The plot includes labels with category names and
##'   percentages, and a legend
##' @docType methods
##' @name plotAnnoPie
##' @rdname plotAnnoPie-methods
##' @seealso \code{\link{plotAnnoBar}} for bar chart visualization,
##'   \code{\link{vennpie}} for nested pie chart visualization,
##'   \code{\link{annotatePeak}} for peak annotation
##' @export
setGeneric("plotAnnoPie",
  function(x,
    ndigit=2,
    cex=0.9,
    col=NA,
    legend.position="rightside",
    pie3D=FALSE,
    radius=0.8,
    ...)
  standardGeneric("plotAnnoPie")
)
