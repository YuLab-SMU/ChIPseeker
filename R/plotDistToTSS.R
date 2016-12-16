
##' plot feature distribution based on the distances to the TSS
##'
##'
##' @title plotDistToTSS.data.frame
##' @param peakDist peak annotation
##' @param distanceColumn column name of the distance from peak to nearest gene
##' @param xlab x label
##' @param ylab y lable
##' @param title figure title
##' @param categoryColumn category column
##' @return bar plot that summarize distance from peak to
##' TSS of the nearest gene.
##' @importFrom magrittr %<>%
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_fill_brewer
##' @importFrom ggplot2 scale_fill_hue
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 geom_text
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, TxDb=txdb)
##' plotDistToTSS(peakAnno)
##' }
##' @seealso \code{\link{annotatePeak}}
##' @author Guangchuang Yu \url{http://ygc.name}
plotDistToTSS.data.frame <- function(peakDist,
                                     distanceColumn="distanceToTSS",
                                     xlab="",
                                     ylab="Binding sites (%) (5'->3')",
                                     title="Distribution of transcription factor-binding loci relative to TSS",
                                     categoryColumn) {

    ## to satisfy codetools
    Feature <- freq <- .id <- NULL

    ## assign Feature according to the distancetoFeature
    peakDist$Feature <- NA
    limit <- c(0, 1000, 3000, 5000, 10000, 100000)
    lbs <- c("0-1kb", "1-3kb", "3-5kb", "5-10kb", "10-100kb", ">100kb")
    for (i in 1:length(limit)) {
        if (i < length(limit)) {
            peakDist$Feature[ abs(peakDist[, distanceColumn]) >= limit[i] & abs(peakDist[,distanceColumn]) < limit[i+1] ] <- lbs[i]
	} else {
            peakDist$Feature[abs(peakDist[,distanceColumn]) > limit[i]] <- lbs[i]
	}
    }
    peakDist$Feature <- factor(peakDist$Feature, levels=lbs)

    ## sign containing -1 and 1 for upstream and downstream
    peakDist$sign <- sign(peakDist[,distanceColumn])

    ## count frequencies
    if (categoryColumn == 1) {
        peakDist %<>% group_by(Feature, sign) %>%
            summarise(freq = length(Feature))

        peakDist$freq = peakDist$freq/sum(peakDist$freq)
        peakDist$freq = peakDist$freq * 100
    } else {
        peakDist %<>% group_by(.id, Feature, sign) %>%
            summarise(freq = length(Feature)) %>%
                group_by(.id) %>%
                    mutate(freq = freq/sum(freq) * 100)
    }

    if (any(peakDist$sign == 0)) {
        zeroDist <- peakDist[peakDist$sign == 0,]
        zeroDist$freq <- zeroDist$freq/2
        zeroDist$sign <- -1
        peakDist[peakDist$sign == 0,] <- zeroDist
        zeroDist$sign <- 1
        peakDist <- rbind(peakDist, zeroDist)
    }

    if (categoryColumn == 1) {
        peakDist %<>% group_by(Feature, sign) %>%
            summarise(freq = sum(freq))

        totalFreq <- peakDist %>% group_by(sign) %>%
            summarise(total = sum(freq))
    } else {
        peakDist %<>% group_by(.id, Feature, sign) %>%
            summarise(freq = sum(freq))
        totalFreq <- peakDist %>% group_by(.id, sign) %>%
            summarise(total = sum(freq))
    }


    ## preparing ylim and y tick labels
    ds = max(totalFreq$total[totalFreq$sign == 1])
    dslim = ceiling(ds/10) * 10
    us = max(totalFreq$total[totalFreq$sign == -1])
    uslim = ceiling(us/10) * 10
    ybreaks <- seq(-uslim, dslim, by=10)
    ylbs <- abs(ybreaks)
    ylbs[ylbs == 0] <- "TSS"

    peakDist$Feature <- factor(peakDist$Feature, levels=rev(levels(peakDist$Feature)))
    if (categoryColumn == 1) {
        p <- ggplot(peakDist, aes(x=1, fill=Feature))
    } else {
        p <- ggplot(peakDist, aes_string(x=categoryColumn, fill="Feature"))
    }

    p <- p + geom_bar(data=subset(peakDist, sign==1), aes(y=freq), stat="identity") +
        geom_bar(data=subset(peakDist, sign==-1), aes(y=-freq), stat="identity")

    p <- p + geom_hline(yintercept = 0, colour = "black") +
        coord_flip() + theme_bw() +
            scale_y_continuous(breaks=ybreaks,labels=ylbs)

    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

    if (categoryColumn == 1) {
        p <- p + scale_x_continuous(breaks=NULL)
    }

    cols <- c("#9ecae1", "#3182bd", "#C7A76C", "#86B875", "#39BEB1", "#CD99D8")
    ## p <- p + scale_fill_hue("Feature", breaks=lbs, labels=lbs)
    ## p <- p + scale_fill_manual(values=getCols(length(lbs)), breaks=lbs, labels=lbs)
    p <- p + scale_fill_manual(values=rev(cols), breaks=rev(lbs), labels=rev(lbs), guide=guide_legend(reverse=TRUE))

    return(p)
}


