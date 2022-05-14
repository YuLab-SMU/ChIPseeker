##' plot the profile of peaks
##' 
##' this function combined previous function plotAvgProf()
##'
##' @title plotPeakProf
##' @param tagMatrix tagMatrix or a list of tagMatrix
##' @param conf confidence interval
##' @param xlab x label
##' @param ylab y label
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param ... additional parameter
##' @return ggplot object
##' @importFrom ggplot2 rel
##' @export
plotPeakProf <- function(tagMatrix,
                         conf,
                         xlab="Genomic Region (5'->3')",
                         ylab = "Peak Count Frequency",
                         facet="none", 
                         free_y = TRUE,
                         ...){
  
  if(is(tagMatrix, "list")){
    upstream <- attr(tagMatrix[[1]], 'upstream')
    downstream <- attr(tagMatrix[[1]], 'downstream')
    label <- attr(tagMatrix[[1]], 'label')
    attr(tagMatrix, 'type') <- attr(tagMatrix[[1]], 'type')
    attr(tagMatrix, 'is.binning') <- attr(tagMatrix[[1]], 'is.binning')
    
  }else{
    upstream <- attr(tagMatrix, 'upstream')
    downstream <- attr(tagMatrix, 'downstream')
    label <- attr(tagMatrix, 'label')
  }
  
  
  if(attr(tagMatrix, 'is.binning')){
    
    plotAvgProf.binning(tagMatrix = tagMatrix, 
                        xlab = xlab,
                        ylab = ylab,
                        conf = conf,
                        facet = facet, 
                        free_y = free_y,
                        upstream = upstream,
                        downstream = downstream,
                        label = label,
                        ...)
    
    
  }else{
    
    xlim <- c(-upstream, downstream)
    plotAvgProf (tagMatrix = tagMatrix, 
                 xlim = xlim,
                 xlab = xlab,
                 ylab = ylab,
                 conf = conf,
                 facet = facet, 
                 free_y = free_y,
                 origin_label = label,
                 ...)
    
  }
  
  
}


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
##' @param origin_label label of the center
##' @param verbose print message or not
##' @param ... additional parameter
##' @return ggplot object
##' @author G Yu; Y Yan
##' @export
plotAvgProf <- function(tagMatrix, xlim,
                        xlab="Genomic Region (5'->3')",
                        ylab = "Peak Count Frequency",
                        conf,
                        facet="none", 
                        free_y = TRUE, 
                        origin_label = "TSS",
                        verbose = TRUE,
                        ...) {
  
  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)

  if (verbose) {
      cat(">> plotting figure...\t\t\t",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  
  conf <- if(missingArg(conf)) NA else conf
  
  if (!(missingArg(conf) || is.na(conf))){
    p <- plotAvgProf.internal(tagMatrix = tagMatrix, 
                              conf = conf, 
                              xlim = xlim,
                              xlab = xlab, 
                              ylab = ylab,
                              facet = facet, 
                              free_y = free_y, 
                              origin_label = origin_label,
                              ...)
  } else {
    p <- plotAvgProf.internal(tagMatrix, 
                              xlim = xlim,
                              xlab = xlab, 
                              ylab = ylab,
                              facet = facet, 
                              free_y = free_y, 
                              origin_label = origin_label,
                              ...)
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
plotAvgProf.internal <- function(tagMatrix, conf,
                                 xlim = c(-3000,3000),
                                 xlab = "Genomic Region (5'->3')",
                                 ylab = "Peak Count Frequency",
                                 facet="none", 
                                 free_y = TRUE,
                                 origin_label, 
                                 ...) {
  
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
                                labels=c(paste0(xlim[1],"bp"), paste0(floor(xlim[1]/2),"bp"),
                                         origin_label, 
                                         paste0(floor(xlim[2]/2),"bp"), paste0(xlim[2], "bp")))
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
##' @param ignore_strand ignore the strand information or not
##' @param ... additional parameter
##' @return ggplot object
##' @export
##' @author G Yu, Ming L
plotAvgProf2 <- function(peak, weightCol = NULL, TxDb = NULL,
                         upstream = 1000, downstream = 1000,
                         xlab = "Genomic Region (5'->3')",
                         ylab = "Peak Count Frequency",
                         conf,
                         facet = "none",
                         free_y = TRUE,
                         verbose = TRUE, 
                         ignore_strand = FALSE,
                         ...) {
  
  plotPeakProf2(peak = peak, 
                upstream = upstream, 
                downstream = downstream,
                conf,
                by = "gene",
                type = "start_site",
                weightCol = weightCol, 
                TxDb = TxDb,
                xlab = xlab,
                ylab = ylab,
                facet = facet,
                free_y = free_y,
                verbose = verbose, 
                nbin = 800,
                ignore_strand = ignore_strand,
                ...)
  
}

##' plot the profile of peaks  by binning
##' 
##' 
##' @title plotAvgProf.binning
##' @param tagMatrix tagMatrix or a list of tagMatrix
##' @param xlab x label
##' @param ylab y label 
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled 
##' @param upstream rel object reflects the percentage of flank extension, e.g rel(0.2)
##'                 integer reflects the actual length of flank extension or TSS region
##'                 NULL reflects the gene body with no extension
##' @param downstream rel object reflects the percentage of flank extension, e.g rel(0.2)
##'                   integer reflects the actual length of flank extension or TSS region
##'                   NULL reflects the gene body with no extension
##' @param label label
##' @param ... additional parameter
##' @return ggplot object
##' @importFrom ggplot2 rel
plotAvgProf.binning <- function(tagMatrix, 
                                xlab = "Genomic Region (5'->3')",
                                ylab = "Peak Count Frequency",
                                conf,
                                facet ="none", 
                                free_y = TRUE,
                                upstream = NULL,
                                downstream = NULL,
                                label,
                                ...) {
  
  ## S4Vectors change the behavior of ifelse
  ## see https://support.bioconductor.org/p/70871/
  ##
  ## conf <- ifelse(missingArg(conf), NA, conf)
  conf <- if(missingArg(conf)) NA else conf
  
  if (!(missingArg(conf) || is.na(conf))){
    p <- plotAvgProf.binning.internal(tagMatrix , 
                                      conf = conf, 
                                      xlab = xlab, 
                                      ylab = ylab,
                                      facet = facet, 
                                      free_y = free_y,
                                      upstream = upstream,
                                      downstream = downstream,
                                      label = label,
                                      ...)
  } else {
    p <- plotAvgProf.binning.internal(tagMatrix , 
                                      xlab = xlab, 
                                      ylab = ylab,
                                      facet = facet, 
                                      free_y = free_y, 
                                      upstream = upstream,
                                      downstream = downstream,
                                      label = label,
                                      ...)
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
##' @importFrom ggplot2 rel
plotAvgProf.binning.internal <- function(tagMatrix, 
                                         conf,
                                         xlab = "Genomic Region (5'->3')",
                                         ylab = "Peak Count Frequency",
                                         facet="none", 
                                         free_y = TRUE,
                                         upstream = NULL,
                                         downstream = NULL,
                                         label,
                                         ...) {
  
  listFlag <- FALSE
  if (is(tagMatrix, "list")) {
    if ( is.null(names(tagMatrix )) ) {
      nn <- paste0("peak", seq_along(tagMatrix ))
      warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
      names(tagMatrix) <- nn
      ## stop("tagMatrix should be a named list...")
    }
    listFlag <- TRUE
  }
  
  if(listFlag){
    nbin <- dim(tagMatrix[[1]])[2]
  }else{
    nbin <- dim(tagMatrix)[2]
  }
  xlim <- c(1,nbin)
  
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
    tagCount <- lapply(tagMatrix , function(x) getTagCount(x, xlim = xlim, conf = conf, ...))
    tagCount <- list_to_dataframe(tagCount)
    tagCount$.id <- factor(tagCount$.id, levels=names(tagMatrix ))
    p <- ggplot(tagCount, aes(pos, group=.id, color=.id))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = .id),
                           linetype = 0, alpha = 0.2)
    }
  } else {
    tagCount <- getTagCount(tagMatrix , xlim = xlim, conf = conf, ...)
    p <- ggplot(tagCount, aes(pos))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper),
                           linetype = 0, alpha = 0.2)
    }
  }
  
  p <- p + geom_line(aes(y = value))
  
  ## x_scale for genebody
  if(attr(tagMatrix, 'type') == 'body'){
    ## x_scale for gene body with no flank extension
    if(is.null(upstream)){
      p <- p + scale_x_continuous(breaks=c(1, 
                                           floor(nbin*0.25),
                                           floor(nbin*0.5),
                                           floor(nbin*0.75),
                                           nbin),
                                  labels=c(label[1], 
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2]))
    }
    
    
    ## x_scale for flank extension by relative value
    if(inherits(upstream, 'rel')){
      
      p <- p + scale_x_continuous(breaks=c(1, 
                                           floor(nbin*(as.numeric(upstream)*100/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+25)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+50)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+75)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           floor(nbin*((as.numeric(upstream)*100+100)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                                           nbin),
                                  labels=c(paste0("-",as.numeric(upstream)*100,"%"), 
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0("+",as.numeric(downstream)*100,"%")))
      p <- p + geom_vline(xintercept=floor(nbin*(as.numeric(upstream)*100/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                          linetype="longdash")
      
      p <- p + geom_vline(xintercept=floor(nbin*((as.numeric(upstream)*100+100)/(100+(as.numeric(upstream)+as.numeric(downstream))*100))),
                          linetype="longdash")
    }
    
    ## x_scale for flank extension by absolute value
    if(!is.null(upstream) & !inherits(upstream, 'rel')){
      
      upstreamPer <- floor(upstream/1000)*0.1
      downstreamPer <- floor(downstream/1000)*0.1
      
      p <- p + scale_x_continuous(breaks=c(1, 
                                           floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.25)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.5)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+0.75)/(1+upstreamPer+downstreamPer))),
                                           floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))),
                                           nbin),
                                  labels=c(paste0("-",upstream,"bp"), 
                                           label[1],
                                           "25%",
                                           "50%",
                                           "75%",
                                           label[2],
                                           paste0(downstream,"bp")))
      p <- p + geom_vline(xintercept=floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
                          linetype="longdash")
      
      p <- p + geom_vline(xintercept=floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))),
                          linetype="longdash")
    }
  }
  
  
  ## x_scale for start region
  if(attr(tagMatrix, 'type') != 'body'){
    
    p <- p + scale_x_continuous(breaks=c(1, 
                                         floor(nbin*0.25),
                                         floor(nbin*0.5),
                                         floor(nbin*0.75),
                                         nbin),
                                labels=c(paste0("-",upstream,"bp"), 
                                         paste0("-",floor(upstream*0.5),"bp"),
                                         label,
                                         paste0(floor(downstream*0.5),"bp"),
                                         paste0(downstream,"bp")))
    
    p <- p + geom_vline(xintercept=floor(nbin*0.5),
                        linetype="longdash")
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


##' plot the profile of peaks automatically
##'
##'`\code{plotPeakProf2()} will call \code{getTagMatrix()} function.
##' \code{getTagMatrix()} function can produce the matrix for visualization.
##' \code{peak} stands for the peak file. \code{window} stands for a collection of regions
##' that users want to look into. Users can use \code{window} to capture the peak of interest.
##' There are two ways to input \code{window}. 
##' 
##' The first way is that users can use
##' \code{getPromoters()/getBioRegion()/makeBioRegionFromGranges()} to get \code{window} and
##' put it into \code{getTagMatrix()}. 
##' 
##' The second way is that users can use \code{getTagMatrix()} to
##' call \code{getPromoters()/getBioRegion()/makeBioRegionFromGranges()}. In this way
##' users do not need to input \code{window} parameter but they need to input
##' \code{txdb} or \code{gr (self-made granges object)}. 
##' 
##' \code{txdb} is a set of packages contained annotation 
##' of regions of different genomes. Users can
##' get the regions of interest through specific functions. These specific functions
##' are built in \code{getPromoters()/getBioRegion()}. Many regions can not be gain
##' through \code{txdb}, like insulator and enhancer regions. Users can provide these
##' regions in the form of granges object. These self-made granges object will be passed
##' to \code{makeBioRegionFromGranges()} to produce the \code{window}.
##' 
##' Details see \code{\link{getPromoters}},\code{\link{getBioRegion}} and \code{\link{makeBioRegionFromGranges}}
##' 
##' \code{upstream} and \code{downstream} parameter have different usages:
##' 
##' (1) \code{window} parameter is provided, 
##' 
##' if \code{type == 'body'}, \code{upstream} and \code{downstream} can use to extend 
##' the flank of body region.
##' 
##' if \code{type == 'start_site'/'end_site'}, \code{upstream} and \code{downstream} do not
##' play a role in \code{getTagMatrix()} function.
##' 
##' (2) \code{window} parameter is missing,
##' 
##' if \code{type == 'body'}, \code{upstream} and \code{downstream} can use to extend 
##' the flank of body region.
##' 
##' if \code{type == 'start_site'/'end_site'}, \code{upstream} and \code{downstream} refer to
##' the upstream and downstream of the start_site or the end_site.
##' 
##' \code{weightCol} refers to column in peak file. This column acts as a weight vaule. Details
##' see \url{https://github.com/YuLab-SMU/ChIPseeker/issues/15}
##' 
##' \code{nbin} refers to the number of bins. \code{getTagMatrix()} provide a binning method
##' to get the tag matrix.
##'
##' @title plotPeakProf2
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight
##' @param TxDb TxDb object
##' @param gr self-made granges object
##' @param upstream upstream position
##' @param downstream downstream position
##' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR'
##' @param type one of "start_site", "end_site", "body"
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param verbose print message or not
##' @param nbin the amount of nbines 
##' @param ignore_strand ignore the strand information or not
##' @param ... additional parameter
##' @return ggplot object
##' @export
##' @author G Yu, Ming Li
plotPeakProf2 <- function(peak, 
                          upstream, 
                          downstream,
                          conf,
                          by,
                          type,
                          weightCol = NULL, 
                          TxDb = NULL,
                          gr,
                          xlab = "Genomic Region (5'->3')",
                          ylab = "Peak Count Frequency",
                          facet = "none",
                          free_y = TRUE,
                          verbose = TRUE, 
                          nbin = NULL,
                          ignore_strand = FALSE,
                          ...){
  
  if ( is(peak, "list") ) {
    tagMatrix <- lapply(peak, getTagMatrix, 
                        upstream = upstream,
                        downstream = downstream, 
                        type = type,
                        TxDb = TxDb,
                        gr = gr,
                        by = by,
                        weightCol = weightCol, 
                        nbin = nbin,
                        verbose = verbose,
                        ignore_strand = ignore_strand)
  } else {
    tagMatrix <- getTagMatrix(peak = peak, 
                              upstream = upstream,
                              downstream = downstream, 
                              type = type,
                              by = by,
                              gr = gr,
                              TxDb = TxDb,
                              weightCol = weightCol, 
                              nbin = nbin,
                              verbose = verbose,
                              ignore_strand = ignore_strand)
  }
  
  
  if (!(missingArg(conf) || is.na(conf))){
    p <- plotPeakProf(tagMatrix = tagMatrix,
                      conf = conf,
                      xlab = xlab,
                      ylab = ylab,
                      facet = facet, 
                      free_y = free_y,
                      ...)
    
  } else {
    p <- plotPeakProf(tagMatrix = tagMatrix,
                      xlab = xlab,
                      ylab = ylab,
                      facet= facet, 
                      free_y = free_y,
                      ...)
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
  promoter <- getBioRegion(TxDb=TxDb,
                           upstream=upstream,
                           downstream=downstream,
                           by="gene",
                           type="start_site")
  
  if (verbose) {
    cat(">> preparing tag matrix...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  if (listFlag) {
    tagMatrix <- lapply(peak, getTagMatrix, weightCol=weightCol, windows=promoter)
  } else {
    tagMatrix <- getTagMatrix(peak, weightCol = weightCol, windows = promoter)
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

