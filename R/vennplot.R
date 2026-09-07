##' Plot Venn diagram for overlapping sets
##'
##' This function creates Venn diagrams to visualize overlaps between multiple
##' sets of objects (e.g., peaks, genes, genomic regions). It supports three
##' different plotting methods with different features and customization options.
##'
##' @description
##' The function calculates all possible intersections between sets using the
##' \code{overlap()} function, then visualizes the results as a Venn diagram.
##' It supports three plotting backends, each with different features:
##' \itemize{
##'   \item \code{gplots}: Default method, produces simple black-and-white Venn
##'         diagrams using the \code{gplots} package
##'   \item \code{ggVennDiagram}: Produces colorful, customizable Venn diagrams
##'         using ggplot2, allowing full control over colors, labels, and styling
##'   \item \code{Vennerable}: Alternative method using the Vennerable package
##'         (requires installation from R-Forge)
##' }
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Validates input: ensures the list is named (auto-generates names if
##'         missing, e.g., "Set1", "Set2", ...)
##'   \item Calculates all possible intersections using \code{overlap()}, which
##'         handles both vectors and GRanges objects
##'   \item Creates the Venn diagram using the specified method:
##'     \itemize{
##'       \item \code{gplots}: Uses \code{plot.venn()} with a matrix format
##'       \item \code{ggVennDiagram}: Passes sets directly to
##'             \code{ggVennDiagram::ggVennDiagram()} with additional parameters
##'       \item \code{Vennerable}: Creates a Venn object and plots it
##'     }
##' }
##'
##' For GRanges objects, the function uses set intersection to find overlapping
##' genomic regions. For vectors, it uses standard set intersection.
##'
##' @param Sets named list of objects to compare. Each element can be:
##'   \itemize{
##'     \item A vector (e.g., gene IDs, peak names)
##'     \item A GRanges object (genomic regions)
##'   }
##'   If the list is not named, names will be auto-generated as "Set1", "Set2", etc.
##' @param by character, plotting method to use. One of:
##'   \itemize{
##'     \item "gplots" (default): Simple black-and-white Venn diagrams
##'     \item "ggVennDiagram": Colorful, customizable ggplot2-based diagrams
##'     \item "Vennerable": Alternative method (requires Vennerable package)
##'   }
##' @param ... additional parameters passed to \code{ggVennDiagram::ggVennDiagram()}
##'   when \code{by="ggVennDiagram"}. Common parameters include:
##'   \itemize{
##'     \item \code{show_percentage}: logical, whether to show percentages
##'     \item \code{label_alpha}: numeric, transparency of labels
##'     \item \code{edge_size}: numeric, size of circle edges
##'   }
##'   See \code{\link[ggVennDiagram]{ggVennDiagram}} for all available options
##' @return A Venn diagram plot. The return type depends on the \code{by} parameter:
##'   \itemize{
##'     \item \code{gplots}: Base R plot (invisible return)
##'     \item \code{ggVennDiagram}: ggplot2 object that can be further customized
##'     \item \code{Vennerable}: Grid plot (invisible return)
##'   }
##' @importFrom gplots plot.venn
## @importFrom ggVennDiagram ggVennDiagram
## @importFrom Vennerable Venn
## @importFrom grid grid.newpage
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfiles <- getSampleFiles()
##' peakAnnoList <- lapply(peakfiles, annotatePeak, TxDb=txdb)
##' names(peakAnnoList) <- names(peakfiles)
##'
##' ## Compare genes from different peak files
##' genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
##' vennplot(genes)
##'
##' ## Use ggVennDiagram for customizable plot
##' vennplot(genes, by="ggVennDiagram", show_percentage=TRUE)
##'
##' ## Compare GRanges objects directly
##' peaks <- lapply(peakfiles, readPeakFile)
##' names(peaks) <- names(peakfiles)
##' vennplot(peaks)
##' }
##' @seealso \code{\link{overlap}} for calculating set intersections,
##'   \code{\link{vennplot.peakfile}} for a convenience function for peak files
##' @export
##' @author G Yu
vennplot <- function(Sets, by="gplots",...) {
    if (is.null(names(Sets))) {
        nn <- paste0("Set", seq_along(Sets))
        warning("input is not a named list, set the name automatically to ", paste(nn, collapse = " "))
        names(Sets) <- nn
        ## stop("input object should be a named list...")
    }

    overlapDF <- overlap(Sets)
    if (by == "Vennerable") {
        ## setRepositories(ind=7)
        ## install.package("Vennerable")
        ## OR
        ## install.packages("Vennerable", repos="http://R-Forge.R-project.org")
        pkg <- "Vennerable"
        require(pkg, character.only=TRUE)
        Venn <- eval(parse(text="Venn"))
        v <- Venn(SetNames=names(Sets), Weight=overlapDF$Weight)
        plotVenn <- eval(parse(text="Vennerable:::plotVenn"))
        plotVenn(v)
    } else if (by == "gplots") {
        n <- ncol(overlapDF)
        colnames(overlapDF)[n] <- "num"
        overlapDF <- overlapDF[, c(n, 1:(n-1))]
        rownames(overlapDF)=apply(overlapDF, 1, function(i) paste(i[-1], sep="", collapse=""))
        vennCount <- as.matrix(overlapDF)
        class(vennCount) <- "venn"
        plot.venn(vennCount)
    } else if(by == "ggVennDiagram"){
	    ggVennDiagram::ggVennDiagram(Sets, ...)
    } else {
        stop("not supported...")
    }
}

##' Plot Venn diagram for peak files
##'
##' This is a convenience function that reads peak files and creates a Venn
##' diagram showing the overlap of genomic regions between different peak sets.
##'
##' @description
##' The function simplifies the workflow of comparing multiple peak files by
##' automatically reading the files and setting up the Venn diagram. It loads
##' each peak file using \code{readPeakFile()}, optionally labels them, and then
##' calls \code{vennplot()} to create the visualization.
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Reads each peak file using \code{readPeakFile()} (converts to GRanges)
##'   \item Generates labels from filenames if not provided (removes file extensions)
##'   \item Names the peak list with the labels
##'   \item Calls \code{vennplot()} to create the Venn diagram
##' }
##'
##' This is particularly useful for quickly comparing multiple ChIP-seq peak
##' files to see how many peaks are shared between experiments.
##'
##' @param files character vector of peak file paths. Supported formats include
##'   BED, narrowPeak, broadPeak, and other formats supported by
##'   \code{readPeakFile()}
##' @param labels character vector of labels for each peak file. If NULL, labels
##'   are automatically generated from filenames by removing file extensions.
##'   Default is NULL
##' @return A Venn diagram plot showing the overlap of peaks between different
##'   files. The plot is created using the default method (\code{by="gplots"}).
##'   See \code{\link{vennplot}} for details on the return type
##' @export
##' @examples
##' \dontrun{
##' ## Compare multiple peak files
##' peakfiles <- c("sample1.bed", "sample2.bed", "sample3.bed")
##' vennplot.peakfile(peakfiles)
##'
##' ## With custom labels
##' vennplot.peakfile(peakfiles, labels=c("Treatment", "Control", "Input"))
##'
##' ## Using getSampleFiles()
##' files <- getSampleFiles()
##' vennplot.peakfile(files)
##' }
##' @seealso \code{\link{vennplot}} for the main plotting function,
##'   \code{\link{readPeakFile}} for reading peak files
##' @author G Yu
vennplot.peakfile <- function(files, labels=NULL) {
    peak.Sets <- lapply(files, readPeakFile)
    if (is.null(labels)) {
        ## remove .xls or .bed of the file names as labels
        labels <- sub("\\.\\w+$", "", files)
    }
    names(peak.Sets) <- labels
    vennplot(peak.Sets)
}


