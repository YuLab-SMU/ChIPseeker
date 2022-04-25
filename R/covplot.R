
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
##' @param lower lower cutoff of coverage signal
##' @return ggplot2 object
##' @import GenomeInfoDb
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 geom_blank
##' @importFrom ggplot2 geom_rect
##' @importFrom ggplot2 facet_grid
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 ggtitle
##' @export
##' @author G Yu
covplot <- function(peak, weightCol=NULL,
                    xlab  = "Chromosome Size (bp)",
                    ylab  = "",
                    title = "ChIP Peaks over Chromosomes",
                    chrs  = NULL,
                    xlim  = NULL,
                    lower = 1) {

    isList <- FALSE
    if(is(peak, "GRanges") || length(peak) == 1) {
        tm <- getChrCov(peak=peak, weightCol=weightCol, chrs, xlim, lower=lower)
    } else {
        isList <- TRUE
        ltm <- lapply(peak, getChrCov, weightCol=weightCol, chrs=chrs, xlim=xlim, lower=lower)
        if (is.null(names(ltm))) {
            nn <- paste0("peak", seq_along(ltm))
            warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
            names(ltm) <- nn
        }
        tm <- list_to_dataframe(ltm)
        chr.sorted <- sortChrName(as.character(unique(tm$chr)))
        tm$chr <- factor(tm$chr, levels = chr.sorted)
    }
    
    chr <- start <- end <- value <- .id <- NULL
   
    if(length(tm$chr) == 0){
        p <- ggplot(data.frame(x = 1)) + geom_blank()
    } else {
        p <- ggplot(tm, aes(start, value))
        
        ## p <- p + geom_segment(aes(x=start, y=0, xend=end, yend= value))
        if (isList) {
            p <- p + geom_rect(aes(xmin=start, ymin=0, xmax=end, ymax=value, fill=.id, color=.id)) 
        } else {
            p <- p + geom_rect(aes(xmin=start, ymin=0, xmax=end, ymax=value), fill='black', color='black')
        }
        if(length(unique(tm$chr)) > 1) {
            p <- p + facet_grid(chr ~., scales="free")
        }
    }
    
    p <- p + theme_classic()
    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
    p <- p + scale_y_continuous(expand=c(0,0))
    p <- p + theme(strip.text.y=element_text(angle=360))

    if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) && length(xlim) == 2) {
        p <- p + xlim(xlim)
    }

    return(p)
}

##' @import S4Vectors IRanges
##' @importFrom dplyr group_by
##' @importFrom dplyr summarise
##' @importFrom magrittr %>%
getChrCov <- function(peak, weightCol, chrs, xlim, lower=1) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        peak.gr <- readPeakFile(peak, as="GRanges")
    } else {
        stop("peak should be a GRanges object or a peak file...")
    }

    if ( is.null(weightCol)) {
        peak.cov <- coverage(peak.gr)
    } else {
        weight <- mcols(peak.gr)[[weightCol]]
        peak.cov <- coverage(peak.gr, weight=weight)
    }

    cov <- lapply(peak.cov, slice, lower=lower)

    get.runValue <- function(x) {
        y <- runValue(x)
        sapply(y@listData, mean)
        ## value <- x@subject@values
        ## value[value != 0]
    }

    chr <- start <- end <- cnt <- NULL
    
    ldf <- lapply(1:length(cov), function(i) {
        x <- cov[[i]]
        if (length(x@ranges) == 0) {
            msg <- paste0(names(cov[i]),
                          " dosen't contain signal higher than ",
                          lower)
            message(msg)
            return(NA)
        }
        data.frame(chr   = names(cov[i]),
                   start = start(x),
                   end   = end(x),
                   cnt   = get.runValue(x)
                                        # the following versions are more slower
                                        # unlist(runValue(x)) 
                                        # sapply(x, runValue)
                   )
    })

    ldf <- ldf[!is.na(ldf)]
    df <- do.call("rbind", ldf)
    
    chr.sorted <- sortChrName(as.character(unique(df$chr)))
    df$chr <- factor(df$chr, levels=chr.sorted)
    if (!is.null(chrs) && !all(is.na(chrs)) && all(chrs %in% chr.sorted)) {
        df <- df[df$chr %in% chrs, ]
    }
    if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) && length(xlim) == 2) {
        df <- df[df$start >= xlim[1] & df$end <= xlim[2],]
    }

    df2 <- group_by(df, chr, start, end) %>% summarise(value=sum(cnt))
    return(df2)
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
