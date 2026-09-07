##' Create Venn pie chart for genomic annotation distribution
##'
##' This function creates a nested pie chart (venn diagram style) visualization
##' showing the distribution of peaks across different genomic annotation
##' categories, with hierarchical nesting to show overlapping annotations.
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
##' categories, as peaks can be annotated with multiple features simultaneously
##' (e.g., a peak can be both in an Exon and in a Promoter region).
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Extracts detailed genomic annotation from the \code{csAnno} object
##'   \item Calculates counts for each annotation category:
##'     \itemize{
##'       \item Genic: peaks overlapping exons or introns
##'       \item Intergenic: peaks not in genic regions
##'       \item Exon: peaks overlapping exons
##'       \item Intron: peaks overlapping introns
##'       \item Upstream: peaks in promoter regions with negative distance to TSS
##'       \item Downstream: peaks downstream of genes
##'       \item Distal Intergenic: intergenic peaks far from genes
##'     }
##'   \item Creates three nested floating pie charts using \code{floating.pie()}
##'         from the \code{plotrix} package
##'   \item Adds a legend showing all annotation categories with their colors
##' }
##'
##' The function uses pseudo-counts (+1) for each category to ensure proper
##' visualization even when some categories have zero counts, preventing color
##' mismatches in the pie charts.
##'
##' @param x \code{csAnno} object containing annotated peaks with
##'   \code{detailGenomicAnnotation} slot populated
##' @param r numeric, initial radius for the base pie chart. Controls the overall
##'   size of the plot. The nested pies use multiples of this radius (2*r, 3*r, 4*r).
##'   Default is 0.2
##' @param cex numeric, character expansion factor for the legend text. Larger
##'   values make the legend text bigger. Default is 1.2
##' @param col named character vector, custom colors for annotation categories.
##'   Names should match category names: "Genic", "Intergenic", "Intron", "Exon",
##'   "Upstream", "Downstream", "Distal_Intergenic". Colors not specified will
##'   use default colors. Default is NULL (uses default color scheme)
##' @return No return value. Creates a plot showing nested pie charts of genomic
##'   annotation distribution
##' @importFrom plotrix floating.pie
##' @seealso \code{\link{plotAnnoPie}} for a simple pie chart,
##'   \code{\link{plotAnnoBar}} for a bar chart visualization
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, TxDb=txdb)
##'
##' ## Create venn pie chart
##' vennpie(peakAnno)
##'
##' ## Customize colors
##' vennpie(peakAnno, col=c(Exon="red", Intron="blue"))
##'
##' ## Adjust size and legend
##' vennpie(peakAnno, r=0.3, cex=1.5)
##' }
##' @author G Yu
vennpie.csAnno <- function(x,
                           r = 0.2,
                           cex = 1.2,
                           col = NULL) {
    detailGenomicAnnotation <- x@detailGenomicAnnotation

    distance <- as.data.frame(x)$distanceToTSS
    total <- nrow(detailGenomicAnnotation)
    Genic <- sum(detailGenomicAnnotation$genic)

    Intergenic <- total-Genic
    Distal_Intergenic <- sum(detailGenomicAnnotation$distal_intergenic)
    Intron <- sum(detailGenomicAnnotation$Intron)
    Exon <- sum(detailGenomicAnnotation$Exon)
    Upstream <- sum(detailGenomicAnnotation$Promoter & distance < 0)

    ## fiveUTR <- sum(detailGenomicAnnotation$fiveUTR)
    ## threeUTR <- sum(detailGenomicAnnotation$threeUTR)
    Downstream <- sum(detailGenomicAnnotation$downstream)

    ## fiveUTR='#e5f5e0',threeUTR='#a1d99b',
    cols <- c(NO='white', Genic='#3182bd', Intergenic='#fec44f',
              Intron='#fc9272', Exon='#9ecae1', Upstream='#ffeda0',
              Downstream='#fee0d2', Distal_Intergenic='#d95f0e')

    cols[names(col)] <- col


    ##par(mai = c(0,0,0,0))
    ##layout(matrix(c(1,2), ncol=2), widths=c(0.7,0.3))
    pie(1, radius=r, init.angle=90, col="white", border=NA, labels='')

    ## https://www.biostars.org/p/326456/
    ## if count is 0, floating pie will ignore it
    ## and the color will mismatch with the category
    ## fixed by adding pseudo-count +1
    floating.pie(0,0, c(Exon,
                        Genic-Exon,
                        Distal_Intergenic,
                        Downstream,
                        Intergenic-Distal_Intergenic-Downstream
                        ) + 1,
                 radius=4*r,
                 startpos=pi/2,
                 col=cols[c("Exon", "NO", "NO", "Downstream", "NO")],
                 border=NA)

    floating.pie(0,0, c(Genic-Intron,
                        Intron,
                        Distal_Intergenic,
                        Intergenic-Upstream-Distal_Intergenic,
                        Upstream) +1 ,
                 radius=3*r,
                 startpos=pi/2,
                 col=cols[c("NO", "Intron", "Distal_Intergenic",
                     "NO", "Upstream")],
                 border=NA)

    floating.pie(0, 0, c(Genic, Intergenic) +1,
                 radius=2*r,
                 startpos=pi/2,
                 col=cols[c("Genic", "Intergenic")],
                 border=NA)
    ##plot.new()
    ##legend(center), legend=names(cols)[-1], fill=cols[-1], bty="n")
    legend(3*r, 3*r, legend=sub("_", " ", names(cols)[-1]),
           fill=cols[-1], bty="n", cex=cex)
}
