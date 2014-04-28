
##' plot the Peak Regions over Chromosomes
##'
##' 
##' @title plotChrCov
##' @param peak peak file or GRanges object
##' @param weightCol weight column of peak
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @return ggplot2 object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 facet_grid
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##' @importFrom plyr ldply
##' @importFrom GenomicRanges seqlengths
##' @export
##' @author G Yu
plotChrCov <- function(peak, weightCol=NULL,
                       xlab = "Chromosome Size (bp)",
                       ylab = "",
                       title = "ChIP Peaks over Chromosomes") {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        peak.gr <- readPeakFile(peak, as="GRanges")
    } else {
        stop("peak should be a GRanges object or a peak file...")
    }
    
    tm <- getChrCov(peak.gr=peak.gr, weightCol=weightCol)
    
    pos <- cnt <- chr <- NULL
    p <- ggplot(tm, aes(pos,cnt))
    p <- p + geom_segment(aes(x=pos, y=0, xend=pos, yend=cnt))
    
    p <- p + facet_grid(chr ~., scales="free")
    p <- p + theme_classic()
    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
    
    p <- p + theme(strip.text.y=element_text(angle=360))

    return(p)
}

##' @importFrom GenomicRanges GRanges
##' @importFrom GenomicRanges elementMetadata
##' @importFrom GenomicRanges seqnames
##' @importFrom IRanges elementLengths
##' @importFrom IRanges IRanges
##' @importFrom IRanges Views
##' @importFrom IRanges viewApply
##' @importFrom IRanges as.vector
##' @importFrom Matrix Matrix
##' @importFrom Matrix summary
getChrCov <- function(peak.gr, weightCol) {

    if ( is.null(weightCol)) {
        peak.cov <- coverage(peak.gr)
    } else {
        weight <- elementMetadata(peak.gr)[[weightCol]]
        peak.cov <- coverage(peak.gr, weight=weight)
    }
    seqLen <- lapply(peak.cov, length)
    chrs <- GRanges(seqnames=names(seqLen),
                    ranges=IRanges(rep(1, length(seqLen)),unlist(seqLen)),
                    strand="*")
    
    peakView <- Views(peak.cov, as(chrs, "RangesList"))
    tagMatrixList <- lapply(peakView, function(x) t(viewApply(x, as.vector)))

    tm <- list()
    nn <- names(tagMatrixList)
    for (i in 1:length(tagMatrixList)) {
        cat(">> processing chromosome ", nn[i])
        if (nchar(nn[i]) > 4) {
            cat("\t\t")
        } else {
            cat("\t")
        }
        cat(format(Sys.time(), "%Y-%m-%d %X"), "\n")
        
        M <- summary(Matrix(tagMatrixList[[i]]))
        tm[[i]] <- data.frame(chr=nn[i],
                              pos=M$j,
                              cnt=M$x)
    }

    tmd <- do.call("rbind", tm)

    ## sort chromosome name
    chr.sorted <- sortChrName(names(seqLen))
    tmd$chr <- factor(tmd$chr, levels=chr.sorted)
  
    return(tmd)
}

sortChrName <- function(chr.name) {
    ## X, Y and M will cause warnings when change to number.
    noChr <- suppressWarnings(as.numeric(sub("chr", "", chr.name)))
    ## index of chromosome name are character, such as X, Y
    ch.idx <- which(is.na(noChr))
    
    n.idx <- which(!is.na(noChr))
    chr.n <- noChr[n.idx]

    chr.sorted <- chr.name[n.idx][order(chr.n)]
    if (length(ch.idx) != 0) {
        chr.sorted <- c(chr.sorted, chr.name[ch.idx])
    }
         
    return(chr.sorted)
}
