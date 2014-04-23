##' plot the profile of peaks that align to flank sequences of TSS 
##'
##' 
##' @title plotPeakProf
##' @param peak peak file or GRanges object
##' @param TranscriptDb TranscriptDb object
##' @param upstream upstream position
##' @param downstream downstream position
##' @param xlab xlab
##' @param ylab ylab
##' @return ggplot object
##' @export
##' @author G Yu
plotPeakProf <- function(peak, TranscriptDb=NULL,
                         upstream=1000, downstream=1000,
                         xlab="Genomic Region (5'->3')",
                         ylab="Read count Per Million mapped reads") {
   
    
    promoter <- getPromoters(TranscriptDb=TranscriptDb,
                             upstream=upstream, downstream=downstream)

    tagMatrix <- getTagMatrix(peak, promoter)

    p <- plotPeakProf.internal(tagMatrix, xcoord=c(-upstream, downstream),
                               xlab=xlab, ylab=ylab)
    return(p)
}

##' plot the heatmap of peaks align to flank sequences of TSS
##'
##' 
##' @title plotPeakHeatmap
##' @param peak peak file or GRanges object
##' @param TranscriptDb TranscriptDb object
##' @param upstream upstream position
##' @param downstream downstream position
##' @param color color
##' @return figure
##' @export
##' @author G Yu
plotPeakHeatmap <- function(peak, TranscriptDb=NULL,
                            upstream=1000, downstream=1000,
                            color="red") {
    promoter <- getPromoters(TranscriptDb=TranscriptDb,
                             upstream=upstream, downstream=downstream)

    tagMatrix <- getTagMatrix(peak, promoter)

    plotPeakHeatmap.internal(tagMatrix, color)
}

##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
plotPeakProf.internal <- function(tagMatrix,
                         xcoord=c(-1000,1000),
                         xlab="Genomic Region (5'->3')",
                         ylab="Read count Per Million mapped reads") {
    if ( (xcoord[2]-xcoord[1]+1) != ncol(tagMatrix) )
        stop("please specify appropreate xcoord...")
    
    ss <- colSums(tagMatrix)/1e6
    ## plot(1:length(ss), ss, type="l", xlab=xlab, ylab=ylab)
    pos <- value <- NULL
    dd <- data.frame(pos=c(xcoord[1]:xcoord[2]), value=ss)
    p <- ggplot(dd, aes(pos, value)) + geom_line()
    if ( 0 > xcoord[1] && 0 < xcoord[2] ) {
        p <- p + geom_vline(xintercept=0,
                            linetype="longdash")
        p <- p + scale_x_continuous(breaks=c(xcoord[1], floor(xcoord[1]/2),
                                       0,
                                       floor(xcoord[2]/2), xcoord[2]),
                                   labels=c(xcoord[1], floor(xcoord[1]/2),
                                       "TSS",
                                       floor(xcoord[2]/2), xcoord[2]))
    }
    p <- p+xlab(xlab)+ylab(ylab)
    p <- p + theme_bw()
    return(p)
}

## @importFrom stats kmeans
##' @importFrom grDevices colorRampPalette
##' @importFrom pheatmap pheatmap
plotPeakHeatmap.internal <- function(tagMatrix, color="red") {
    k <- kmeans(tagMatrix, 3)
    ii <- order(rowSums(tagMatrix), decreasing=TRUE)
    kc <- k$cluster[ii]    
    
    ## tdf <- as.data.frame(tagMatrix)
    ## tdf <- cbind(tdf[order(kc),], id=seq(nrow(tdf)))    
    ## tdf$idsort <- tdf$id[order(tdf$cluster)]
    ## tdfm <- melt(tdf, id.vars="id")
    ## p <- ggplot(tdfm, aes(x=variable, y=id))
    ## p <- p + geom_tile(aes(fill=value), color="white")
    ## p <- p + scale_fill_gradient(low="white", high="red")

    hmcols <- colorRampPalette(c("white",color))(50)    
    
    pheatmap(tagMatrix[order(kc),], color=hmcols,
             cluster_rows=F,cluster_cols=FALSE, legend=F,
             show_rownames=FALSE,show_colnames=FALSE)
}


##' extract promoter positions from TranscriptDb.
##'
##' 
##' @title getPromoters
##' @param TranscriptDb TranscriptDb 
##' @param upstream upstream position
##' @param downstream downstream position
##' @return GRanges object
##' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
##' @importFrom BiocGenerics unique
##' @importFrom GenomicRanges GRanges
##' @importFrom GenomicRanges unlist
##' @importFrom GenomicRanges strand
##' @importFrom GenomicFeatures transcriptsBy
##' @importFrom IRanges start
##' @importFrom IRanges end
##' @importFrom IRanges IRanges
##' @author G Yu
getPromoters <- function(TranscriptDb=NULL, upstream=1000, downstream=1000) {
    
    if ( is.null(TranscriptDb) ) {
        TranscriptDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    
    .ChIPseekerEnv(TranscriptDb)
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

    if ( exists("upstream", envir=ChIPseekerEnv, inherits=FALSE) &&
        exists("downstream", envir=ChIPseekerEnv, inherits=FALSE) ) {
        us <- get("upstream", envir=ChIPseekerEnv)
        ds <- get("downstream", envir=ChIPseekerEnv)
        if (us == upstream && ds == downstream &&
            exists("promoters", envir=ChIPseekerEnv, inherits=FALSE) ){
            promoters <- get("promoters", envir=ChIPseekerEnv)
            return(promoters)
        }
    }

    assign("upstream", upstream, envir=ChIPseekerEnv)
    assign("downstream", downstream, envir=ChIPseekerEnv)

    if ( exists("Transcripts", envir=ChIPseekerEnv, inherits=FALSE) ) {
        Transcripts <- get("Transcripts", envir=ChIPseekerEnv)
    } else {
        Transcripts <- transcriptsBy(TranscriptDb)
        Transcripts <- unlist(Transcripts)
        assign("Transcripts", Transcripts, envir=ChIPseekerEnv)
    }
    ## get start position based on strand
    tss <- ifelse(strand(Transcripts) == "+", start(Transcripts), end(Transcripts))
    promoters <- GRanges(seqnames=seqnames(Transcripts),
                         ranges=IRanges(tss-upstream, tss+downstream),
                         strand=strand(Transcripts))
    promoters <- unique(promoters)
    assign("promoters", promoters, envir=ChIPseekerEnv)
    
    return(promoters)
}

##' @importFrom IRanges subsetByOverlaps
##' @importFrom IRanges elementLengths
##' @importFrom IRanges width
##' @importFrom IRanges coverage
##' @importFrom IRanges IRanges
##' @importFrom IRanges Views
##' @importFrom IRanges viewApply
##' @importFrom IRanges as.vector
##' @importFrom IRanges as.factor
##' @importFrom GenomicRanges GRanges
##' @importFrom GenomicRanges seqnames
##' @importFrom BiocGenerics intersect
##' @importFrom BiocGenerics unique
getTagMatrix <- function(peak, windows) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        peak.gr <- readPeakFile(peak, as="GRanges")
    } else {
        stop("peak should be a GRanges object or a peak file...")
    }
    
    if (! is(windows, "GRanges")) {
        stop("windows should be a GRanges object...")
    }
    if (length(unique(width(windows))) != 1) {
        stop("width of windows should be equal...")
    }

    if (!exists("ChIPseekerEnv", envir = .GlobalEnv)) {
        assign("ChIPseekerEnv", new.env(), .GlobalEnv)
    }
    ChIPseekerEnv <- get("ChIPseekerEnv", envir = .GlobalEnv)
    
    if (exists("peak", envir=ChIPseekerEnv, inherits=FALSE) &&
        exists("promoters", envir=ChIPseekerEnv, inherits=FALSE) &&
        exists("tagMatrix", envir=ChIPseekerEnv, inherits=FALSE) ) {

        pp <- get("peak", envir=ChIPseekerEnv)
        promoters <- get("promoters", envir=ChIPseekerEnv)

        if (all(pp == peak)) {
            if (all(windows == promoters)) {
                tagMatrix <- get("tagMatrix", envir=ChIPseekerEnv)
                return(tagMatrix)
            } else {
                assign("promoters", windows)
                ## make sure it is not conflict with getPromoters
                if ( exists("upstream", envir=ChIPseekerEnv, inherits=FALSE))
                    rm("upstream", envir=ChIPseekerEnv)
            }
        } else {
            assign("peak", peak, envir=ChIPseekerEnv)
        }
        
    }

    if ( !exists("peak", envir=ChIPseekerEnv, inherits=FALSE)) {
        assign("peak", peak, envir=ChIPseekerEnv)
    }

    if ( !exists("promoters", envir=ChIPseekerEnv, inherits=FALSE)) {
        assign("promoters", windows, envir=ChIPseekerEnv)
    }
        
    peak.cov <- coverage(peak.gr)
    cov.len <- elementLengths(peak.cov)
    cov.width <- GRanges(seqnames=names(cov.len),
                         IRanges(start=rep(1, length(cov.len)),
                                 end=cov.len))
    windows <- subsetByOverlaps(windows, cov.width,
                                type="within", ignore.strand=TRUE)

    chr.idx <- intersect(names(peak.cov),
                         unique(as.character(seqnames(windows))))
    
    peakView <- Views(peak.cov[chr.idx], as(windows, "RangesList")[chr.idx])
    tagMatrixList <- lapply(peakView, function(x) t(viewApply(x, as.vector)))
    tagMatrix <- do.call("rbind", tagMatrixList)

    ## get the index of windows, that are reorganized by as(windows, "RangesList")
    idx.list <- split(1:length(windows),  as.factor(seqnames(windows)))
    idx <- do.call("c", idx.list)
    
    rownames(tagMatrix) <- idx
    tagMatrix <- tagMatrix[order(idx),]
    
    ## minus strand
    minus.idx <- which(as.character(strand(windows)) == "-")
    tagMatrix[minus.idx,] <- tagMatrix[minus.idx, ncol(tagMatrix):1]
    tagMatrix <- tagMatrix[rowSums(tagMatrix)!=0,]
    assign("tagMatrix", tagMatrix, envir=ChIPseekerEnv)
    return(tagMatrix)
}
