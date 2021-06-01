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
##' @param color color like c("RdYlBu","RdYlGn","Spectral","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")
##' @return figure
##' @export
##' @author G Yu
tagHeatmap <- function(tagMatrix, xlim, xlab="", ylab="", title=NULL, color="RdYlBu") {
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
##' @param color color like c("RdYlBu","RdYlGn","Spectral","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")
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
            cols <- getHeatmapColor(nc)
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
      
        ## create a vector x for the layout
        ## the vector is like 1,1,2,3,3,4.....
        x <- vector(mode="numeric",length=0)
        v <- 1

        for(i in 1:(nc*2)){
          if(i%%2!=0){
            x[v] <- i
            v <- v+1
            x[v] <- i
            v <- v+1
          }else{
            x[v] <- i
            v <- v+1
          }
        }
      
        ## set the layout
        layout(matrix(x,ncol = nc*3))
        layout.show(nc*2)
      
        for (i in 1:nc) {
            ## this par is for the heatmap
            par(mar=c(3,2,2,2))
            peakHeatmap.internal(tagMatrix[[i]], xlim, cols[i], xlab[i], ylab[i], title[i], listFlag)
        }
    } else {
        if (is.null(color) || is.na(color))
            color <- "RdYlBu"
        if (is.null(title) || is.na(title))
            title <- ""
        
        layout(matrix(c(1,1,2),ncol = 3))
        layout.show(2)
      
        peakHeatmap.internal(tagMatrix, xlim, color, xlab, ylab, title, listFlag)
    }
}


##' @import BiocGenerics
##' @importFrom RColorBrewer brewer.pal
##' @importFrom grDevices colorRampPalette
##' @importFrom sinkr imageScale
peakHeatmap.internal <- function(tagMatrix, xlim=NULL, color="RdYlBu", xlab="", ylab="", title="", listFlag=FALSE) {
  
    tagMatrix <- t(apply(tagMatrix, 1, function(x) (x/max(x))*10))
    ii <- order(rowSums(tagMatrix))
    tagMatrix <- tagMatrix[ii,]
                         
    breaks <- seq(0,10,length.out = 1000)
    cols <- colorRampPalette(brewer.pal(8,color))
                                             
    if (is.null(xlim)) {
        xlim <- 1:ncol(tagMatrix)
    } else if (length(xlim) == 2) {
        xlim <- seq(xlim[1], xlim[2])
    }
    image(x=xlim, y=1:nrow(tagMatrix),z=t(tagMatrix),useRaster=TRUE, col=cols(length(breaks)-1), yaxt="n", ylab="", xlab=xlab, main=title)
    
    ## this par is set for the scale                 
    if(listFlag){
      par(pin=c(0.07,1.7))
    }else{
      par(pin=c(0.1,1.7)) 
    }
    
    ## a scale is added to graph 
    imageScale(tagMatrix, col=cols(length(breaks)-1), breaks=breaks, axis.pos=4, add.axis=FALSE)
    axis(4,at=c(0,2,4,6,8,10), las=2)    
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

