##' plot feature distribution based on their chromosome region
##'
##' plot chromosome region features
##' @title plotAnnoBar.data.frame
##' @param anno.df annotation stats
##' @param xlab xlab
##' @param ylab ylab
##' @param title plot title
##' @param categoryColumn category column
##' @return bar plot that summarize genomic features of peaks
##' @importFrom plyr ldply
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
##' @seealso \code{\link{annotatePeak}} \code{\link{plotAnnoPie}}
##' @author Guangchuang Yu \url{http://ygc.name}
plotAnnoBar.data.frame <- function(anno.df,
                                   xlab="",
                                   ylab="Percentage(%)",
                                   title="Feature Distribution",
                                   categoryColumn) {

    
    p <- ggplot(anno.df, aes_string(x = categoryColumn,
                                    fill = "Feature",
                                    y = "Frequency"))
    
    p <- p + geom_bar(stat="identity") + coord_flip() + theme_bw()
    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

    if (categoryColumn == 1) {
        p <- p + scale_x_continuous(breaks=NULL)
        p <- p+scale_fill_manual(values=getCols(nrow(anno.df)))
    } else {
        p <- p+scale_fill_manual(values=getCols(length(unique(anno.df$Feature))))
    }
    
    return(p)
}

##' pieplot from peak genomic annotation
##'
##' 
##' @title plotAnnoPie
##' @param x csAnno object
##' @param ndigit number of digit to round
##' @param cex label cex
##' @param col color
##' @param legend.position topright or other.
##' @param pie3D plot in 3D or not
##' @param ... extra parameter
##' @return pie plot of peak genomic feature annotation
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg38.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="chipseeker")
##' peakAnno <- annotatePeak(peakfile, TxDb=txdb)
##' plotAnnoPie(peakAnno)
##' }
##' @seealso \code{\link{annotatePeak}} \code{\link{plotAnnoBar}}
##' @export
##' @author G Yu
plotAnnoPie.csAnno <- function(x,
                        ndigit=2,
                        cex=0.9,
                        col=NA,
                        legend.position="rightside",
                        pie3D=FALSE,
                        ...){
    
    anno.df <- getAnnoStat(x)
    if (is.na(col)) {
        col <- getCols(nrow(anno.df))
    }
    
    if (pie3D)
        annoPie3D(anno.df, ndigit=ndigit, cex=cex, col=col, ...)
    
    annoPie(anno.df, ndigit=ndigit, cex=cex, col=col, legend.position=legend.position, ...)
 }

##' @importFrom RColorBrewer brewer.pal
##' @importFrom grDevices colorRampPalette
annoPie <- function(anno.df, ndigit=2, cex=0.9, col=NA, legend.position, ...) {
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
        legend("center", legend = labels, fill=col, bty="n")
    } else {
        pie(anno.df$Frequency,
            ##     ## labels=paste(round(anno.df$Frequency/sum(anno.df$Frequency)*100, 2), "%", sep=""),
            labels=paste(anno.df$Feature, " (",
                round(anno.df$Frequency/sum(anno.df$Frequency)*100, ndigit),
                "%)", sep=""),
            cex=cex,
            col=col,
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

    anno[grep("exon 1 of", anno)] <- "1st Exon"
    anno[grep("intron 1 of", anno)] <- "1st Intron"
    anno[grep("Exon \\(", anno)] <- "Other Exon"
    anno[grep("Intron \\(", anno)] <- "Other Intron"
    anno[grep("Downstream", anno)] <- "Downstream (<=3kb)"
    
    ## count frequency
    anno.table <- table(anno)
    
    ## calculate ratio
    anno.ratio <- anno.table/ sum(anno.table) * 100
    anno.df <- as.data.frame(anno.ratio)
    colnames(anno.df) <- c("Feature", "Frequency")
    lvs <- c("Promoter (<=1kb)",
             "Promoter (1-2kb)",
             "Promoter (2-3kb)",
             "5' UTR",
             "3' UTR",
             "1st Exon",
             "Other Exon",
             "1st Intron",
             "Other Intron",
             "Downstream (<=3kb)",
             "Distal Intergenic")
    anno.df$Feature <- factor(anno.df$Feature, levels=lvs[lvs %in% anno.df$Feature])
    anno.df <- anno.df[order(anno.df$Feature),]
    return(anno.df)
}




