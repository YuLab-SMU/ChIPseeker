##' plot the overlap of a list of object
##'
##' 
##' @title vennplot 
##' @param Sets a list of object, can be vector or GRanges object
##' @param by one of gplots or Vennerable
##' @return venn plot that summarize the overlap of peaks
##' from different experiments or gene annotation from
##' different peak files.
##' @importFrom gplots plot.venn
## @importFrom Vennerable Venn
## @importFrom grid grid.newpage
## @importFrom RColorBrewer brewer.pal
##' @examples
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, TranscriptDb=txdb)
##' peakAnnoList <- lapply(1:3, function(i) peakAnno[sample(1:length(peakAnno), 100),])
##' names(peakAnnoList) <- paste("peak", 1:3, sep="_")
##' genes= lapply(peakAnnoList, function(i) unlist(i$geneId))
##' vennplot(genes)
##' @export
##' @author G Yu
vennplot <- function(Sets, by="gplots") {
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
    } else {
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


