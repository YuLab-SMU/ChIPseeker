##' plot the overlap of a list of object
##'
##'
##' There are two ways to plot, which users can specify through `by`.
##' 
##' The first way is to use `gplots` packages, by setting `by = gplots`. This method
##' is default method. The venn plot produced through this way has no color.
##' 
##' The second way is to use `ggVennDiagram` packages, by setting `by = ggVennDiagram`. 
##' The venn plot produced through this way has colors which can be defined by users using
##' ggplot2 grammar e.g.(scale_fill_distiller()). And users can specify any details, like digital number,
##' text size and showing percentage or not, by inputting `...` extra parameters.
##' 
##' @title vennplot
##' @param Sets a list of object, can be vector or GRanges object
##' @param by one of gplots, ggVennDiagram or Vennerable
##' @param ... extra parameters using ggVennDiagram. Details see \link[ggVennDiagram]{ggVennDiagram}
##' @return venn plot that summarize the overlap of peaks
##' from different experiments or gene annotation from
##' different peak files.
##' @importFrom gplots plot.venn
## @importFrom ggVennDiagram ggVennDiagram
## @importFrom Vennerable Venn
## @importFrom grid grid.newpage
## @importFrom RColorBrewer brewer.pal
##' @examples
##' ## example not run
##' ## require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' ## txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' ## peakfiles <- getSampleFiles()
##' ## peakAnnoList <- lapply(peakfiles, annotatePeak)
##' ## names(peakAnnoList) <- names(peakfiles)
##' ## genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
##' ## vennplot(genes)
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
    }
    else {
        stop("not supported...")
    }
}

##' vennplot for peak files
##'
##'
##' @title vennplot.peakfile
##' @param files peak files
##' @param labels labels for peak files
##' @return figure
##' @export
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


