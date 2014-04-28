##' plot feature distribution based on their chromosome region
##'
##' plot chromosome region features
##' @title plotAnnoBar
##' @param peakAnno peakAnno in data.frame
##' @param title plot title
##' @param xlab xlab
##' @param ylab ylab
##' @return bar plot that summarize genomic features of peaks
##' @importFrom plyr ldply
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##' @examples
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, TranscriptDb=txdb)
##' plotAnnoBar(peakAnno)
##' @seealso \code{\link{annotatePeak}} \code{\link{plotAnnoPie}}
##' @export
##' @author G Yu
plotAnnoBar <- function(peakAnno,
                        title="Feature Distribution",
                        xlab="",
                        ylab="Percentage(%)") {
    
    if ( class(peakAnno) == "data.frame" || class(peakAnno) == "GRanges" ) {
        anno.df <- getGenomicAnnoStat(peakAnno)
        categoryColumn <- 1
    } else if ( class(peakAnno) == "list" ) {
        anno <- lapply(peakAnno, getGenomicAnnoStat)
        anno.df <- ldply(anno)
        categoryColumn <- ".id"
    } else {
        stop("peakAnno should be a data.frame, a GRanges object or a named list of data.frame.")
    }
    
    p <- ggplot(anno.df, aes_string(x = categoryColumn,
                                    fill = "Feature",
                                    y = "Frequency"))
    
    p <- p + geom_bar(stat="identity") + coord_flip() + theme_bw()
    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

    if (categoryColumn == 1) {
        p <- p + scale_x_continuous(breaks=NULL)
    }

    ## p <- p+scale_fill_manual(values=getCols(nrow(anno.df)))
    
    return(p)
}

##' pieplot from peak genomic annotation
##'
##' 
##' @title plotAnnoPie
##' @param peakAnno peakAnno
##' @param ndigit number of digit to round
##' @param cex label cex
##' @param col color
##' @param pie3D plot in 3D or not
##' @param ... extra parameter
##' @return pie plot of peak genomic feature annotation
##' @examples
##' ## example not run
##' ## require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' ## txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' ## peakfile <- system.file("extdata", "sample_peaks.txt", package="chipseeker")
##' ## peakAnno <- annotatePeak(peakfile, TranscriptDb=txdb)
##' ## plotAnnoPie(peakAnno)
##' @seealso \code{\link{annotatePeak}} \code{\link{plotAnnoBar}}
##' @export
##' @author G Yu
plotAnnoPie <- function(peakAnno,
                        ndigit=2,
                        cex=0.9,
                        col=NA,
                        pie3D=FALSE,
                        ...){
    
    anno.df <- getGenomicAnnoStat(peakAnno)
    if (is.na(col)) {
        col <- getCols(nrow(anno.df))
    }
    
    if (pie3D)
        annoPie3D(anno.df, ndigit=ndigit, cex=cex, col=col, ...)
    
    annoPie(anno.df, ndigit=ndigit, cex=cex, col=col, ...)
 }

##' @importFrom RColorBrewer brewer.pal
##' @importFrom grDevices colorRampPalette
annoPie <- function(anno.df, ndigit=2, cex=0.9, col=NA, ...) {
    if ( ! all(c("Feature", "Frequency") %in% colnames(anno.df))) {
        stop("check your input...")
    }
    
    pie(anno.df$Frequency,
        ## labels=paste(round(anno.df$Frequency/sum(anno.df$Frequency)*100, 2), "%", sep=""),
        labels=paste(anno.df$Feature, " (",
            round(anno.df$Frequency/sum(anno.df$Frequency)*100, ndigit),
            "%)", sep=""),
        cex=cex,
        col=col,
        ...
        )
    ## legend(legend=anno.df$Feature, fill=col, "topright")
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




