
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
##' @export
##' @author G Yu
plotChrCov <- function(peak, weightCol=NULL,
                    xlab = "Chromosome Size (bp)",
                    ylab = "",
                    title = "ChIP Peaks over Chromosomes"){
    .Deprecated("covplot")
}


##' plot peak coverage
##'
##' 
##' @title covplot
##' @param peak peak file or GRanges object
##' @param weightCol weight column of peak
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @param chrs selected chromosomes to plot, all chromosomes by default
##' @param xlim ranges to plot, default is whole chromosome
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
##' @importFrom GenomeInfoDb seqlengths
##' @importFrom data.table data.table
##' @export
##' @author G Yu
covplot <- function(peak, weightCol=NULL,
                    xlab = "Chromosome Size (bp)",
                    ylab = "",
                    title = "ChIP Peaks over Chromosomes",
                    chrs = NULL,
                    xlim = NULL) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        peak.gr <- readPeakFile(peak, as="GRanges")
    } else {
        stop("peak should be a GRanges object or a peak file...")
    }
    
    tm <- getChrCov(peak.gr=peak.gr, weightCol=weightCol, chrs, xlim)
    bin <- floor((max(tm$end) - min(tm$start))/1000)
    if (bin < 1) {
        bin <- 1
    }
    
    tml <- lapply(1:nrow(tm), function(i) {
        x <- tm[i,]
        data.table(chr=x$chr, pos=seq(x$start, x$end, by=bin), cnt=x[,4])
    })
    
    tm2 <- do.call("rbind", tml)
    tm2 <- ddply(tm2, .(chr, pos), transform, value=sum(cnt))
    tm2 <- tm2[,-3]
    tm2 <- unique(tm2)
    tm2 <- as.data.frame(tm2)
    ##colnames(tm2) <- c("chr", "pos", "value")
    
    pos <- chr <- value <- NULL
    
    p <- ggplot(tm2, aes(pos, value))
    p <- p + geom_segment(aes(x=pos, y=0, xend=pos, yend= value))
    if(length(unique(tm$chr)) > 1) {
        p <- p + facet_grid(chr ~., scales="free")
    }
    p <- p + theme_classic()
    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
    p <- p + scale_y_continuous(expand=c(0,0))
    p <- p + theme(strip.text.y=element_text(angle=360))
    
    return(p)
}

##' @importFrom GenomicRanges elementMetadata
##' @importFrom IRanges slice
##' @importFrom S4Vectors runValue
getChrCov <- function(peak.gr, weightCol, chrs, xlim) {

    if ( is.null(weightCol)) {
        peak.cov <- coverage(peak.gr)
    } else {
        weight <- elementMetadata(peak.gr)[[weightCol]]
        peak.cov <- coverage(peak.gr, weight=weight)
    }

    cov <- lapply(peak.cov, slice, lower=1)

    ldf <- lapply(1:length(cov), function(i) {
        x <- cov[[i]]
        data.frame(chr=names(cov[i]),
                   start=start(x),
                   end = end(x),
                   cnt = sapply(x, runValue)
                                        # value <- x@subject@values
                                        # value <- value[value != 0]
                   )
    })
    
    df <- do.call("rbind", ldf)
    chr.sorted <- sortChrName(as.character(unique(df$chr)))
    df$chr <- factor(df$chr, levels=chr.sorted)
    if (!is.null(chrs) && !is.na(chrs) && chrs %in% chr.sorted) {
        df <- df[df$chr %in% chrs, ]
    }
    if (!is.null(xlim) && !is.na(xlim) && is.numeric(xlim) && length(xlim) == 2) {
        df <- df[df$start >= xlim[1] & df$end <= xlim[2],]
    }

    ##colnames(df) <- c("chr", "start", "end", "cnt")
    
    df2 <- ddply(df, .(chr, start, end), transform, value=sum(cnt))
    df2 <- df2[,-4]
    df2 <- unique(df2)
    return(df)
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
        chr.sorted <- c(chr.sorted, sort(chr.name[ch.idx]))
    }
         
    return(chr.sorted)
}
