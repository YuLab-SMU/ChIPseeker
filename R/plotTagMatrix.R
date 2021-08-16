##' plot the profile of peaks
##'
##'
##' @title plotAvgProf
##' @param tagMatrix tagMatrix or a list of tagMatrix
##' @param xlim xlim
##' @param xlab x label
##' @param ylab y label
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param ... additional parameter
##' @return ggplot object
##' @export
##' @author G Yu; Y Yan
plotAvgProf <- function(tagMatrix, xlim,
                        xlab="Genomic Region (5'->3')",
                        ylab = "Peak Count Frequency",
                        conf,
                        facet="none", free_y = TRUE, ...) {
    ## S4Vectors change the behavior of ifelse
    ## see https://support.bioconductor.org/p/70871/
    ##
    ## conf <- ifelse(missingArg(conf), NA, conf)
    ##
    conf <- if(missingArg(conf)) NA else conf

    if (!(missingArg(conf) || is.na(conf))){
        p <- plotAvgProf.internal(tagMatrix, conf = conf, xlim = xlim,
                                  xlab = xlab, ylab = ylab,
                                  facet = facet, free_y = free_y, ...
                                  )
    } else {
        p <- plotAvgProf.internal(tagMatrix, xlim = xlim,
                                  xlab = xlab, ylab = ylab,
                                  facet = facet, free_y = free_y, ...)
    }
    return(p)
}

##' plot the profile of peaks that align to flank sequences of TSS
##'
##'
##' @title plotAvgProf
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight
##' @param TxDb TxDb object
##' @param upstream upstream position
##' @param downstream downstream position
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param verbose print message or not
##' @param ... additional parameter
##' @return ggplot object
##' @export
##' @author G Yu
plotAvgProf2 <- function(peak, weightCol = NULL, TxDb = NULL,
                         upstream = 1000, downstream = 1000,
                         xlab = "Genomic Region (5'->3')",
                         ylab = "Peak Count Frequency",
                         conf,
                         facet = "none",
                         free_y = TRUE,
                         verbose = TRUE, ...) {

    if (verbose) {
        cat(">> preparing promoter regions...\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    promoter <- getPromoters(TxDb=TxDb,
                             upstream=upstream,
                             downstream=downstream)

    if (verbose) {
        cat(">> preparing tag matrix...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    if ( is(peak, "list") ) {
        tagMatrix <- lapply(peak, getTagMatrix,
                            weightCol=weightCol, windows=promoter)
    } else {
        tagMatrix <- getTagMatrix(peak, weightCol, promoter)
    }

    if (verbose) {
        cat(">> plotting figure...\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }

    if (!(missingArg(conf) || is.na(conf))){
        p <- plotAvgProf.internal(tagMatrix,
                                  xlim = c(-upstream, downstream),
                                  xlab = xlab, ylab = ylab, conf = conf,
                                  facet = facet, free_y = free_y, ...)
    } else {
        p <- plotAvgProf.internal(tagMatrix,
                                  xlim=c(-upstream, downstream),
                                  xlab=xlab, ylab=ylab,
                                  facet = facet, free_y = free_y, ...)
    }
    return(p)
}

##' plot the heatmap of tagMatrix
##'
##'
##' @title tagHeatmap
##' @param tagMatrix tagMatrix or a list of tagMatrix
##' @param xlim xlim
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @param color color
##' @return figure
##' @export
##' @author G Yu
tagHeatmap <- function(tagMatrix, xlim, xlab="", ylab="", title=NULL, color="red") {
    listFlag <- FALSE
    if (is(tagMatrix, "list")) {
        listFlag <- TRUE
    }
    peakHeatmap.internal2(tagMatrix, xlim, listFlag, color, xlab, ylab, title)
}

##' plot the heatmap of peaks align to flank sequences of TSS
##'
##'
##' @title peakHeatmap
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight
##' @param TxDb TxDb object
##' @param upstream upstream position
##' @param downstream downstream position
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @param color color
##' @param verbose print message or not
##' @return figure
##' @export
##' @author G Yu
peakHeatmap <- function(peak, weightCol=NULL, TxDb=NULL,
                            upstream=1000, downstream=1000,
                            xlab="", ylab="", title=NULL,
                            color=NULL, verbose=TRUE) {
    listFlag <- FALSE
    if ( is(peak, "list") ) {
        listFlag <- TRUE
        if (is.null(names(peak)))
            stop("peak should be a peak file or a name list of peak files...")
    }

    if (verbose) {
        cat(">> preparing promoter regions...\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    promoter <- getPromoters(TxDb=TxDb,
                             upstream=upstream, downstream=downstream)

    if (verbose) {
        cat(">> preparing tag matrix...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    if (listFlag) {
        tagMatrix <- lapply(peak, getTagMatrix, weightCol=weightCol, windows=promoter)
    } else {
        tagMatrix <- getTagMatrix(peak, weightCol, promoter)
    }

    if (verbose) {
        cat(">> generating figure...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }

    xlim=c(-upstream, downstream)

    peakHeatmap.internal2(tagMatrix, xlim, listFlag, color, xlab, ylab, title)

    if (verbose) {
        cat(">> done...\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    invisible(tagMatrix)
}

peakHeatmap.internal2 <- function(tagMatrix, xlim, listFlag, color, xlab, ylab, title) {
    if ( is.null(xlab) || is.na(xlab))
        xlab <- ""
    if ( is.null(ylab) || is.na(ylab))
        ylab <- ""

    if (listFlag) {
        nc <- length(tagMatrix)
        if ( is.null(color) || is.na(color) ) {
            cols <- getCols(nc)
        } else if (length(color) != nc) {
            cols <- rep(color[1], nc)
        } else {
            cols <- color
        }

        if (is.null(title) || is.na(title))
            title <- names(tagMatrix)
        if (length(xlab) != nc) {
            xlab <- rep(xlab[1], nc)
        }
        if (length(ylab) != nc) {
            ylab <- rep(ylab[1], nc)
        }
        if (length(title) != nc) {
            title <- rep(title[1], nc)
        }
        par(mfrow=c(1, nc))
        for (i in 1:nc) {
            peakHeatmap.internal(tagMatrix[[i]], xlim, cols[i], xlab[i], ylab[i], title[i])
        }
    } else {
        if (is.null(color) || is.na(color))
            color <- "red"
        if (is.null(title) || is.na(title))
            title <- ""
        peakHeatmap.internal(tagMatrix, xlim, color, xlab, ylab, title)
    }
}


##' @import BiocGenerics
##' @importFrom grDevices colorRampPalette
peakHeatmap.internal <- function(tagMatrix, xlim=NULL, color="red", xlab="", ylab="", title="") {
    tagMatrix <- t(apply(tagMatrix, 1, function(x) x/max(x)))
    ii <- order(rowSums(tagMatrix))
    tagMatrix <- tagMatrix[ii,]
    cols <- colorRampPalette(c("white",color))(200)
    if (is.null(xlim)) {
        xlim <- 1:ncol(tagMatrix)
    } else if (length(xlim) == 2) {
        xlim <- seq(xlim[1], xlim[2])
    }
    image(x=xlim, y=1:nrow(tagMatrix),z=t(tagMatrix),useRaster=TRUE, col=cols, yaxt="n", ylab="", xlab=xlab, main=title)
}



##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_ribbon
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 facet_grid
plotAvgProf.internal <- function(tagMatrix, conf,
                                 xlim = c(-3000,3000),
                                 xlab = "Genomic Region (5'->3')",
                                 ylab = "Peak Count Frequency",
                                 facet="none", free_y = TRUE,
                                 origin_label = "TSS", ...) {

    listFlag <- FALSE
    if (is(tagMatrix, "list")) {
        if ( is.null(names(tagMatrix)) ) {
            nn <- paste0("peak", seq_along(tagMatrix))
            warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
            names(tagMatrix) <- nn
            ## stop("tagMatrix should be a named list...")
        }
        listFlag <- TRUE
    }

    if ( listFlag ) {
        facet <- match.arg(facet, c("none", "row", "column"))
        if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix[[1]]) ) {
            stop("please specify appropreate xcoordinations...")
        }
    } else {
        if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix) ) {
            stop("please specify appropreate xcoordinations...")
        }
    }

    ## S4Vectors change the behavior of ifelse
    ## see https://support.bioconductor.org/p/70871/
    ##
    ## conf <- ifelse(missingArg(conf), NA, conf)
    ##
    conf <- if(missingArg(conf)) NA else conf

    pos <- value <- .id <- Lower <- Upper <- NULL

    if ( listFlag ) {
        tagCount <- lapply(tagMatrix, function(x) getTagCount(x, xlim = xlim, conf = conf, ...))
        tagCount <- list_to_dataframe(tagCount)
        tagCount$.id <- factor(tagCount$.id, levels=names(tagMatrix))
        p <- ggplot(tagCount, aes(pos, group=.id, color=.id))
        if (!(is.na(conf))) {
            p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = .id),
                                linetype = 0, alpha = 0.2)
        }
    } else {
        tagCount <- getTagCount(tagMatrix, xlim = xlim, conf = conf, ...)
        p <- ggplot(tagCount, aes(pos))
        if (!(is.na(conf))) {
            p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper),
                                 linetype = 0, alpha = 0.2)
        }
    }

    p <- p + geom_line(aes(y = value))

    if ( 0 > xlim[1] && 0 < xlim[2] ) {
        p <- p + geom_vline(xintercept=0,
                            linetype="longdash")
        p <- p + scale_x_continuous(breaks=c(xlim[1], floor(xlim[1]/2),
                                       0,
                                       floor(xlim[2]/2), xlim[2]),
                                   labels=c(xlim[1], floor(xlim[1]/2),
                                       origin_label, 
                                       floor(xlim[2]/2), xlim[2]))
    }

    if (listFlag) {
        cols <- getCols(length(tagMatrix))
        p <- p + scale_color_manual(values=cols)
        if (facet == "row") {
            if (free_y) {
                p <- p + facet_grid(.id ~ ., scales = "free_y")
            } else {
                p <- p + facet_grid(.id ~ .)
            }
        } else if (facet == "column") {
            if (free_y) {
                p <-  p + facet_grid(. ~ .id, scales = "free_y")
            } else {
                p <-  p + facet_grid(. ~ .id)
            }
        }
    }
    p <- p+xlab(xlab)+ylab(ylab)
    p <- p + theme_bw() + theme(legend.title=element_blank())
    if(facet != "none") {
        p <- p + theme(legend.position="none")
    }
    return(p)
}


##' plot the profile of peaks of the gene regions
##' 
##' 
##' Title plotGeneBody
##'
##' @param bodymatrix tagMatrix tagMatrix or a list of tagMatrix
##' @param xlab x label
##' @param ylab y label 
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled 
##' @param ... 
##' @return ggplot object
##' @export
plotGeneBody <- function(bodymatrix, 
                         xlab= "Scaled Genomic Region (5'->3')",
                         ylab = "Peak Count Frequency",
                         conf,
                         facet="none", free_y = TRUE, 
                         ...) {
    
    ## S4Vectors change the behavior of ifelse
    ## see https://support.bioconductor.org/p/70871/
    ##
    ## conf <- ifelse(missingArg(conf), NA, conf)
    ##
    
    conf <- if(missingArg(conf)) NA else conf
    
    box <- dim(bodymatrix)[2]
    xlim <- c(1,box)
    
    if (!(missingArg(conf) || is.na(conf))){
        p <- plotGeneBody.internal(bodymatrix , conf = conf, xlim = xlim,
                                   xlab = xlab, ylab = ylab,
                                   facet = facet, free_y = free_y,
                                   box=box, ...)
    } else {
        p <- plotGeneBody.internal(bodymatrix , xlim = xlim,
                                   xlab = xlab, ylab = ylab,
                                   facet = facet, free_y = free_y, 
                                   box=box, ...)
    }
    return(p)
}


##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_ribbon
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 facet_grid
plotGeneBody.internal <- function(bodymatrix, conf,
                                  xlim,
                                  xlab = "Scaled Genomic Region (5'->3')",
                                  ylab = "Peak Count Frequency",
                                  facet="none", free_y = TRUE,
                                  box,...) {
    
    listFlag <- FALSE
    if (is(bodymatrix, "list")) {
        if ( is.null(names(bodymatrix )) ) {
            nn <- paste0("peak", seq_along(bodymatrix ))
            warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
            names(bodymatrix) <- nn
            ## stop("tagMatrix should be a named list...")
        }
        listFlag <- TRUE
    }
    
    if ( listFlag ) {
        facet <- match.arg(facet, c("none", "row", "column"))
    }
    
    ## S4Vectors change the behavior of ifelse
    ## see https://support.bioconductor.org/p/70871/
    ##
    ## conf <- ifelse(missingArg(conf), NA, conf)
    ##
    conf <- if(missingArg(conf)) NA else conf
    
    pos <- value <- .id <- Lower <- Upper <- NULL
    
    if ( listFlag ) {
        tagCount <- lapply(bodymatrix , function(x) getTagCount(x, xlim = xlim, conf = conf, ...))
        tagCount <- list_to_dataframe(tagCount)
        tagCount$.id <- factor(tagCount$.id, levels=names(bodymatrix ))
        p <- ggplot(tagCount, aes(pos, group=.id, color=.id))
        if (!(is.na(conf))) {
            p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = .id),
                                 linetype = 0, alpha = 0.2)
        }
    } else {
        tagCount <- getTagCount(bodymatrix , xlim = xlim, conf = conf, ...)
        p <- ggplot(tagCount, aes(pos))
        if (!(is.na(conf))) {
            p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper),
                                 linetype = 0, alpha = 0.2)
        }
    }
    
    p <- p + geom_line(aes(y = value))
    
    p <- p + scale_x_continuous(breaks=c(1, 
                                         floor(box*0.25),
                                         floor(box*0.5),
                                         floor(box*0.75),
                                         box),
                                labels=c("1 box", 
                                         paste(floor(box*0.25)*binsize,"box",seq=" "),
                                         paste(floor(box*0.5)*binsize,"box",seq=" "),
                                         paste(floor(box*0.75)*binsize,"box",seq=" "),
                                         paste(box,"box",seq=" ")))
    
    
    if (listFlag) {
        cols <- getCols(length(bodymatrix))
        p <- p + scale_color_manual(values=cols)
        if (facet == "row") {
            if (free_y) {
                p <- p + facet_grid(.id ~ ., scales = "free_y")
            } else {
                p <- p + facet_grid(.id ~ .)
            }
        } else if (facet == "column") {
            if (free_y) {
                p <-  p + facet_grid(. ~ .id, scales = "free_y")
            } else {
                p <-  p + facet_grid(. ~ .id)
            }
        }
    }
    p <- p+xlab(xlab)+ylab(ylab)
    p <- p + theme_bw() + theme(legend.title=element_blank())
    if(facet != "none") {
        p <- p + theme(legend.position="none")
    }
    return(p)
}



##' plot the profile of peaks that align to gene regions
##' 
##' 
##' Title plotGeneBody2
##'
##' @param peak peak file or GRanges object
##' @param weightCol weightCol column name of weight
##' @param TxDb TxDb object
##' @param type one of "genes", "exon", "intron", "promoters"
##' @param xlab x label
##' @param ylab y label
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param verbose print message or not
##' @param box the amount of boxes needed to be splited which should be less than min_body_length
##' @param min_body_length the minimum length that each gene region should be 
##' @param downstream downstream, when use this function to draw TSS region
##' @param upstream upstream, when use this function to draw TSS region
##' @return ggplot object
##' @export
plotGeneBody2 <- function(peak, weightCol = NULL, TxDb = NULL,
                          type = "genes",
                          xlab = "Scaled Genomic Region (5'->3')",
                          ylab = "Peak Count Frequency",
                          conf,
                          facet = "none",
                          free_y = TRUE,
                          verbose = TRUE, 
                          box=500,
                          min_body_length=1000, 
                          upstream=3000,
                          downstream=3000,...) {
    
    if (verbose) {
        cat(">> preparing genebody regions...\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    
    ## TSS region is also supported
    type <- match.arg(type, c("genes", "exon", "intron", "promoters"))
    if(type=="promoters"){
        genebody <- getPromoters(TxDb=txdb, 
                                 upstream=upstream, downstream=downstream)
    }else{
        genebody <- getGeneBody(TxDb = txdb, type = type)
    }
    
    
    
    if (verbose) {
        cat(">> preparing genebody matrix...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    
    if ( is(peak, "list") ) {
        bodymatrix <- lapply(peak, getGenebodyMatrix,
                             weightCol=weightCol, 
                             windows=genebody,
                             box=box,
                             min_body_length=min_body_length)
        
    } else {
        bodymatrix <- getGenebodyMatrix(peak, weightCol, genebody,
                                        box,
                                        min_body_length)
    }
    
    
    xlim <- c(1,box)
    
    if (verbose) {
        cat(">> plotting figure...\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    
    
    if (!(missingArg(conf) || is.na(conf))){
        p <- plotGeneBody.internal(bodymatrix,
                                   xlim = xlim,
                                   xlab = xlab, ylab = ylab, conf = conf,
                                   facet = facet, free_y = free_y, ...)
    } else {
        p <- plotGeneBody.internal(bodymatrix,
                                   xlim=xlim,
                                   xlab=xlab, ylab=ylab,
                                   facet = facet, free_y = free_y, ...)
    }
    return(p)
}

