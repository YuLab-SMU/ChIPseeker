##' plot feature distribution based on their chromosome region
##'
##' plot chromosome region features
##' @title plotAnnoBar.data.frame
##' @rdname plotAnnoBar
##' @param anno.df annotation stats
##' @param xlab xlab
##' @param ylab ylab
##' @param title plot title
##' @param categoryColumn category column
##' @return bar plot that summarize genomic features of peaks
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 guide_legend
##' @seealso \code{\link{annotatePeak}} \code{\link{plotAnnoPie}}
##' @author Guangchuang Yu \url{https://yulab-smu.top}
plotAnnoBar.data.frame <- function(anno.df,
                                   xlab="",
                                   ylab="Percentage(%)",
                                   title="Feature Distribution",
                                   categoryColumn) {


    anno.df$Feature <- factor(anno.df$Feature, levels = rev(levels(anno.df$Feature)))

    p <- ggplot(anno.df, aes_string(x = categoryColumn,
                                    fill = "Feature",
                                    y = "Frequency"))

    p <- p + geom_bar(stat="identity") + coord_flip() + theme_bw()
    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

    if (categoryColumn == 1) {
        p <- p + scale_x_continuous(breaks=NULL)
        p <- p+scale_fill_manual(values=rev(getCols(nrow(anno.df))), guide=guide_legend(reverse=TRUE))
    } else {
        p <- p+scale_fill_manual(values=rev(getCols(length(unique(anno.df$Feature)))), guide=guide_legend(reverse=TRUE))
    }

    return(p)
}

##' pieplot from peak genomic annotation
##'
##'
##' @title plotAnnoPie
##' @rdname plotAnnoPie
##' @param x csAnno object
##' @param ndigit number of digit to round
##' @param cex label cex
##' @param col color
##' @param legend.position topright or other.
##' @param pie3D plot in 3D or not
##' @param radius radius of Pie
##' @param ... extra parameter
##' @return pie plot of peak genomic feature annotation
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="chipseeker")
##' peakAnno <- annotatePeak(peakfile, TxDb=txdb)
##' plotAnnoPie(peakAnno)
##' }
##' @seealso \code{\link{annotatePeak}} \code{\link{plotAnnoBar}}
##' @export
##' @author Guangchuang Yu \url{https://yulab-smu.top}
plotAnnoPie.csAnno <- function(x,
                        ndigit=2,
                        cex=0.8,
                        col=NA,
                        legend.position="rightside",
                        pie3D=FALSE,
                        radius=0.8,
                        ...){

    anno.df <- getAnnoStat(x)
    if (is.na(col)) {
        col <- getCols(nrow(anno.df))
    }

    if (pie3D)
        annoPie3D(anno.df, ndigit=ndigit, cex=cex, col=col, ...)

    annoPie(anno.df, ndigit=ndigit, cex=cex, col=col, legend.position=legend.position, radius=radius, ...)
 }

##' @importFrom RColorBrewer brewer.pal
##' @importFrom grDevices colorRampPalette
##' @importFrom graphics par
##' @importFrom graphics layout
##' @importFrom graphics pie
##' @importFrom graphics legend
##' @importFrom graphics plot.new
annoPie <- function(anno.df, ndigit=2, cex=0.8, col=NA, legend.position, radius=0.8, ...) {
    if ( ! all(c("Feature", "Frequency") %in% colnames(anno.df))) {
        stop("check your input...")
    }

    if (legend.position == "rightside") {
        labels=paste(anno.df$Feature, " (",
            round(anno.df$Frequency/sum(anno.df$Frequency)*100, ndigit),
            "%)", sep="")

        par(mai = c(0,0,0,0))
        layout(matrix(c(1,2), ncol=2), widths=c(0.6,0.4))
        pie(anno.df$Frequency, labels=NA, cex=cex, col=col, ...)
        plot.new()
        legend("center", legend = labels, fill = col, bty = "n", cex = cex)
    } else {
        par(mai = c(0,0,0,0))
        pie(anno.df$Frequency,
            ##     ## labels=paste(round(anno.df$Frequency/sum(anno.df$Frequency)*100, 2), "%", sep=""),
            labels=paste(anno.df$Feature, " (",
                round(anno.df$Frequency/sum(anno.df$Frequency)*100, ndigit),
                "%)", sep=""),
            cex=cex,
            col=col,
            radius=radius,
            ...
            )
    }
}

## @param ndigit ndigit
## @param radius the radius of the pie
## @param explode the amount to "explode" the pie
## @param labelcex label font size
## @importFrom plotrix pie3D
annoPie3D <- function(anno.df,
                      ndigit=2,
                      cex=1,
                      ...){

    ## anno.df <- getGenomicAnnoStat(peakAnno)

    pkg <- "plotrix"
    require(pkg, character.only=TRUE)
    pie3D <- eval(parse(text="pie3D"))
    pie3D(anno.df$Frequency,
          labels=paste(
              anno.df$Feature,
              "(",
              paste(round(anno.df$Frequency, ndigit), "%", sep=""),
              ")",
              sep=""),
          labelcex=cex,
          col=col,
          ...)
}

getGenomicAnnoStat <- function(peakAnno) {
    if ( class(peakAnno) == "GRanges" )
        peakAnno <- as.data.frame(peakAnno)
    anno <- peakAnno$annotation
    ## anno <- sub(" \\(.+", "", anno)

    e1 <- getOption("ChIPseeker.ignore_1st_exon")
    i1 <- getOption("ChIPseeker.ignore_1st_intron")
    ids <- getOption("ChIPseeker.ignore_downstream")

    if (is.null(e1) || !e1) {
        e1lab <- "1st Exon"
        anno[grep("exon 1 of", anno)] <- e1lab
        exonlab <- "Other Exon"
    } else {
        e1lab <- NULL
        exonlab <- "Exon"
    }

    if (is.null(i1) || !i1) {
        i1lab <- "1st Intron"
        anno[grep("intron 1 of", anno)] <- i1lab
        intronlab <- "Other Intron"
    } else {
        i1lab <- NULL
        intronlab <- "Intron"
    }

    anno[grep("Exon \\(", anno)] <- exonlab
    anno[grep("Intron \\(", anno)] <- intronlab

    if (is.null(ids) || !ids) {
        dsd <- getOption("ChIPseeker.downstreamDistance")
        if (is.null(dsd))
            dsd <- 3000 ## downstream 3k by default
        if (dsd > 1000) {
            dsd <- round(dsd/1000, 1)
            dsd <- paste0(dsd, "kb")
        }
        dslab <- paste0("Downstream (<=", dsd, ")")

        anno[grep("Downstream", anno)] <- dslab
        iglab <- "Distal Intergenic"
    } else {
        dslab <- NULL
        iglab <- "Intergenic"
        anno[grep("Downstream", anno)] <- iglab
    }
    anno[grep("^Distal", anno)] <- iglab

    lvs <- c(
        "5' UTR",
        "3' UTR",
        e1lab,
        exonlab,
        i1lab,
        intronlab,
        dslab,
        iglab
    )

    promoter <- unique(anno[grep("Promoter", anno)])
    ip <- getOption("ChIPseeker.ignore_promoter_subcategory")
    if ((is.null(ip) || !ip) && (length(promoter) > 0)) {
        plab <- sort(as.character(promoter))
    } else {
        plab <- "Promoter"
        anno[grep("^Promoter", anno)] <- plab
    }
    lvs <- c(plab, lvs)

    ## count frequency
    anno.table <- table(anno)

    ## calculate ratio
    anno.ratio <- anno.table/ sum(anno.table) * 100
    anno.df <- as.data.frame(anno.ratio)
    colnames(anno.df) <- c("Feature", "Frequency")

    anno.df$Feature <- factor(anno.df$Feature, levels=lvs[lvs %in% anno.df$Feature])
    anno.df <- anno.df[order(anno.df$Feature),]
    return(anno.df)
}




