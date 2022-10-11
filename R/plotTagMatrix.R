##' plot the profile of peaks
##'`
##' \code{plotPeakProf_MultiWindows()} is almost the same as \code{plotPeakProf2()}, having
##' the main difference of accepting two or more granges objects. Accepting more
##' granges objects can help compare the same peaks in different windows.
##' 
##' \code{TxDb} parameter can accept txdb object.
##' But many regions can not be obtained by txdb object. In this case,
##' Users can provide self-made granges served the same role 
##' as txdb object and pass to \code{TxDb} object.
##' 
##' \code{by} the features of interest. 
##' 
##' (1) if users use \code{txdb}, \code{by} can be one of 'gene', 'transcript', 'exon', 
##' 'intron' , '3UTR' , '5UTR', 'UTR'. These features can be obtained by functions from txdb object.
##' 
##' (2) if users use self-made granges object, \code{by} can be everything. Because this \code{by}
##' will not pass to functions to get features, which is different from the case of using 
##' txdb object. This \code{by} is only used to made labels showed in picture.
##' 
##' \code{type} means the property of the region. one of the "start site",
##' "end site" and "body".
##' 
##' \code{upstream} and \code{downstream} parameter have different usages:
##' 
##' (1) if \code{type == 'body'}, \code{upstream} and \code{downstream} can use to extend 
##' the flank of body region.
##' 
##' (2) if \code{type == 'start_site'/'end_site'}, \code{upstream} and \code{downstream} refer to
##' the upstream and downstream of the start_site or the end_site.
##' 
##' \code{weightCol} refers to column in peak file. This column acts as a weight value. Details
##' see \url{https://github.com/YuLab-SMU/ChIPseeker/issues/15}
##' 
##' \code{nbin} refers to the number of bins. \code{getTagMatrix()} provide a binning method
##' to get the tag matrix.
##' 
##' There are two ways input a list of window.
##' 
##' (1) Users can input a list of self-made granges objects
##' 
##' (2) Users can input a list of \code{by} and only one \code{type}. In this way, 
##' \code{plotPeakProf_MultiWindows()} can made a list of window from txdb object based on \code{by} and \code{type}.
##' 
##' Warning: 
##' 
##' (1) All of these window should be the same type. It means users can only
##' compare a list of "start site"/"end site"/"body region" with the same upstream
##' and downstream.
##' 
##' (2) So it will be only one \code{type} and several \code{by}.
##' 
##' (3) Users can make window by txdb object or self-made granges object. Users can only
##' choose one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR' or 'UTR' in the
##' way of using txdb object. User can input any \code{by} in the way of using 
##' self-made granges object.
##' 
##' (4) Users can mingle the \code{by} designed for the two ways. \code{plotPeakProf_MultiWindows} can
##' accpet the hybrid \code{by}. But the above rules should be followed.
##' 
##' \url{https://github.com/YuLab-SMU/ChIPseeker/issues/189}
##'
##' @title plotPeakProf_MultiWindows
##' 
##' @param tagMatrix tagMatrix or a list of tagMatrix
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight
##' @param TxDb TxDb object or self-made granges objects
##' @param upstream upstream position
##' @param downstream downstream position
##' @param by feature of interest
##' @param type one of "start_site", "end_site", "body"
##' @param windows_name the name for each window, which will also be showed in the picture as labels
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param verbose print message or not
##' @param nbin the amount of bines 
##' @param ignore_strand ignore the strand information or not
##' @param ... additional parameter
##' @return ggplot object
##' @export
plotPeakProf <- function(tagMatrix = NULL,
                         peak,
                         upstream,
                         downstream,
                         conf,
                         by,
                         type,
                         windows_name = NULL,
                         weightCol = NULL,
                         TxDb = NULL,
                         xlab = "Genomic Region (5'->3')",
                         ylab = "Peak Count Frequency",
                         facet = "row",
                         free_y = TRUE,
                         verbose = TRUE,
                         nbin = NULL,
                         ignore_strand = FALSE,
                         ...){
  
  if(is.null(tagMatrix)){
    
    conf <- if(missingArg(conf)) NA else conf
    upstream <- if(missingArg(upstream)) NULL else upstream
    downstream <- if(missingArg(downstream)) NULL else downstream
    
    if(length(by) == 1){
      
      plotPeakProf2(peak = peak, 
                    upstream = upstream, 
                    downstream = downstream,
                    conf = conf,
                    by = by,
                    type = type,
                    weightCol = weightCol, 
                    TxDb = TxDb,
                    xlab = xlab,
                    ylab = ylab,
                    facet = facet,
                    free_y = free_y,
                    verbose = verbose, 
                    nbin = nbin,
                    ignore_strand = ignore_strand,
                    ...)
      
    }else{
      
      if(is.null(windows_name) && !is.null(names(TxDb)))
        windows_name <- names(TxDb)
      
      plotPeakProf_MultiWindows(peak = peak,
                                upstream = upstream,
                                downstream = downstream,
                                conf = conf,
                                by = by,
                                type = type,
                                windows_name = windows_name,
                                weightCol = weightCol,
                                TxDb = TxDb,
                                xlab = xlab,
                                ylab = ylab,
                                facet = facet,
                                free_y = free_y,
                                verbose = verbose,
                                nbin = nbin,
                                ignore_strand = ignore_strand,
                                ...)
      
    }
    
  }else{
    
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
      
      if (!(missingArg(conf) || is.na(conf))){
        
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
        
        plotAvgProf.binning(tagMatrix = tagMatrix, 
                            xlab = xlab,
                            ylab = ylab,
                            facet = facet, 
                            free_y = free_y,
                            upstream = upstream,
                            downstream = downstream,
                            label = label,
                            ...)
        
      }
      
      
    }else{
      
      xlim <- c(-upstream, downstream)
      
      if (!(missingArg(conf) || is.na(conf))){
        
        plotAvgProf (tagMatrix = tagMatrix, 
                     xlim = xlim,
                     xlab = xlab,
                     ylab = ylab,
                     conf = conf,
                     facet = facet, 
                     free_y = free_y,
                     origin_label = label,
                     ...)
        
      }else{
        
        plotAvgProf (tagMatrix = tagMatrix, 
                     xlim = xlim,
                     xlab = xlab,
                     ylab = ylab,
                     facet = facet, 
                     free_y = free_y,
                     origin_label = label,
                     ...)
        
      }
      
      
    }
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
##' This function is the old function of \code{plotPeakProf2}. It can
##' only plot the start site region of gene.
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
##' \code{peak} stands for the peak file. 
##' 
##' \code{by} the features of interest. 
##' 
##' (1) if users use \code{txdb}, \code{by} can be one of 'gene', 'transcript', 'exon', 
##' 'intron' , '3UTR' , '5UTR', 'UTR'. These features can be obtained by functions from txdb object.
##' 
##' (2) if users use self-made granges object, \code{by} can be everything. Because this \code{by}
##' will not pass to functions to get features, which is different from the case of using 
##' txdb object. This \code{by} is only used to made labels showed in picture.
##' 
##' \code{type} means the property of the region. one of the "start site",
##' "end site" and "body".
##' 
##' \code{upstream} and \code{downstream} parameter have different usages:
##' 
##' (1) if \code{type == 'body'}, \code{upstream} and \code{downstream} can use to extend 
##' the flank of body region.
##' 
##' (2) if \code{type == 'start_site'/'end_site'}, \code{upstream} and \code{downstream} refer to
##' the upstream and downstream of the start_site or the end_site.
##' 
##' \code{weightCol} refers to column in peak file. This column acts as a weight vaule. Details
##' see \url{https://github.com/YuLab-SMU/ChIPseeker/issues/15}
##' 
##' \code{nbin} refers to the number of bins, providing a binning method
##' to get the tag matrix.
##' 
##' \code{TxDb} parameter can accept txdb object.
##' But many regions can not be obtained by txdb object. In this case,
##' Users can provide self-made granges served the same role 
##' as txdb object and pass to \code{TxDb} object.
##' 
##' \code{plotPeakProf2()} is different from the \code{plotPeakProf()}. \code{plotPeakProf2()} do not
##' need to provide \code{window} parameter, which means \code{plotPeakProf2()} will call relevent
##' functions to make \code{window} automatically.
##'
##' @title plotPeakProf2
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight
##' @param TxDb TxDb object, or self-made granges object
##' @param upstream upstream position
##' @param downstream downstream position
##' @param by e.g. 'gene', 'transcript', 'exon' or features of interest(e.g. "enhancer")
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
                          xlab = "Genomic Region (5'->3')",
                          ylab = "Peak Count Frequency",
                          facet = "none",
                          free_y = TRUE,
                          verbose = TRUE, 
                          nbin = NULL,
                          ignore_strand = FALSE,
                          ...){
  
  conf <- if(missingArg(conf)) NA else conf
  upstream <- if(missingArg(upstream)) NULL else upstream
  downstream <- if(missingArg(downstream)) NULL else downstream
  
  if ( is(peak, "list") ) {
    tagMatrix <- lapply(peak, getTagMatrix, 
                        upstream = upstream,
                        downstream = downstream, 
                        type = type,
                        TxDb = TxDb,
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


##' plot the profile of peaks in two or more windows
##'
##'
##' This function comes from \url{https://github.com/YuLab-SMU/ChIPseeker/issues/189}
##'`
##' \code{plotPeakProf_MultiWindows()} is almost the same as \code{plotPeakProf2()}, having
##' the main difference of accepting two or more granges objects. Accepting more
##' granges objects can help compare the same peaks in different windows.
##' 
##' \code{TxDb} parameter can accept txdb object.
##' But many regions can not be obtained by txdb object. In this case,
##' Users can provide self-made granges served the same role 
##' as txdb object and pass to \code{TxDb} object.
##' 
##' \code{by} the features of interest. 
##' 
##' (1) if users use \code{txdb}, \code{by} can be one of 'gene', 'transcript', 'exon', 
##' 'intron' , '3UTR' , '5UTR', 'UTR'. These features can be obtained by functions from txdb object.
##' 
##' (2) if users use self-made granges object, \code{by} can be everything. Because this \code{by}
##' will not pass to functions to get features, which is different from the case of using 
##' txdb object. This \code{by} is only used to made labels showed in picture.
##' 
##' \code{type} means the property of the region. one of the "start site",
##' "end site" and "body".
##' 
##' \code{upstream} and \code{downstream} parameter have different usages:
##' 
##' (1) if \code{type == 'body'}, \code{upstream} and \code{downstream} can use to extend 
##' the flank of body region.
##' 
##' (2) if \code{type == 'start_site'/'end_site'}, \code{upstream} and \code{downstream} refer to
##' the upstream and downstream of the start_site or the end_site.
##' 
##' \code{weightCol} refers to column in peak file. This column acts as a weight value. Details
##' see \url{https://github.com/YuLab-SMU/ChIPseeker/issues/15}
##' 
##' \code{nbin} refers to the number of bins. \code{getTagMatrix()} provide a binning method
##' to get the tag matrix.
##' 
##' There are two ways input a list of window.
##' 
##' (1) Users can input a list of self-made granges objects
##' 
##' (2) Users can input a list of \code{by} and only one \code{type}. In this way, 
##' \code{plotPeakProf_MultiWindows()} can made a list of window from txdb object based on \code{by} and \code{type}.
##' 
##' Warning: 
##' 
##' (1) All of these window should be the same type. It means users can only
##' compare a list of "start site"/"end site"/"body region" with the same upstream
##' and downstream.
##' 
##' (2) So it will be only one \code{type} and several \code{by}.
##' 
##' (3) Users can make window by txdb object or self-made granges object. Users can only
##' choose one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR' or 'UTR' in the
##' way of using txdb object. User can input any \code{by} in the way of using 
##' self-made granges object.
##' 
##' (4) Users can mingle the \code{by} designed for the two ways. \code{plotPeakProf_MultiWindows} can
##' accpet the hybrid \code{by}. But the above rules should be followed.
##' 
##'
##' @title plotPeakProf_MultiWindows
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight
##' @param TxDb TxDb object or self-made granges objects
##' @param upstream upstream position
##' @param downstream downstream position
##' @param by feature of interest
##' @param type one of "start_site", "end_site", "body"
##' @param windows_name the name for each window, which will also be showed in the picture as labels
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param verbose print message or not
##' @param nbin the amount of bines 
##' @param ignore_strand ignore the strand information or not
##' @param ... additional parameter
##' @return ggplot object
plotPeakProf_MultiWindows <- function(peak,
                                      upstream,
                                      downstream,
                                      conf,
                                      by,
                                      type,
                                      windows_name = NULL,
                                      weightCol = NULL,
                                      TxDb = NULL,
                                      xlab = "Genomic Region (5'->3')",
                                      ylab = "Peak Count Frequency",
                                      facet = "row",
                                      free_y = TRUE,
                                      verbose = TRUE,
                                      nbin = NULL,
                                      ignore_strand = FALSE,
                                      ...){
  
  conf <- if(missingArg(conf)) NA else conf
  upstream <- if(missingArg(upstream)) NULL else upstream
  downstream <- if(missingArg(downstream)) NULL else downstream
  
  ## check type
  if(length(type) != 1){
    stop("It should be only one type...")
  }
  
  ## make the window name
  if (is.null(windows_name)) {
    nn <- by
    warning("set the name automatically to ", paste(nn, collapse=' '))
    windows_name <- nn
  }else{
    if (length(windows_name) != length(by)) {
      stop("the length of the window name and the by should be equal...")
    }
  }
  
  
  if ( is(peak, "list") ) {
    tagMatrix <- lapply(peak, getTagMatrix2,
                        upstream=upstream,
                        downstream=downstream,
                        windows_name=windows_name,
                        type=type,
                        by=by,
                        TxDb=TxDb,
                        weightCol = weightCol, 
                        nbin = nbin,
                        verbose = verbose,
                        ignore_strand= ignore_strand)
  } else {
    tagMatrix <- getTagMatrix2(peak=peak, 
                               upstream=upstream,
                               downstream=downstream,
                               windows_name=windows_name,
                               type=type,
                               by=by,
                               TxDb=TxDb,
                               weightCol = weightCol, 
                               nbin = nbin,
                               verbose = verbose,
                               ignore_strand= ignore_strand)
  }
  
  if (!(missingArg(conf) || is.na(conf))){
    p <- plotMultiProf(tagMatrix = tagMatrix,
                       conf = conf,
                       xlab = xlab,
                       ylab = ylab,
                       facet = facet, 
                       free_y = free_y,
                       ...)
    
  } else {
    p <- plotMultiProf(tagMatrix = tagMatrix,
                       xlab = xlab,
                       ylab = ylab,
                       facet= facet, 
                       free_y = free_y,
                       ...)
  }
  return(p)
  
}


##' internal function for plotPeakProf_MultiWindows
##' 
##' @param tagMatrix tagMatrix
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param ... additional parameter
plotMultiProf <- function(tagMatrix,
                          conf,
                          xlab="Genomic Region (5'->3')",
                          ylab = "Peak Count Frequency",
                          facet="none", 
                          free_y = TRUE,
                          ...){
  
  
  if(is(tagMatrix[[1]][[1]],"matrix")){
    upstream <- attr(tagMatrix[[1]][[1]], 'upstream')
    downstream <- attr(tagMatrix[[1]][[1]], 'downstream')
    # attr(tagMatrix, 'type') <- attr(tagMatrix[[1]][[1]], 'type')
    # attr(tagMatrix, 'is.binning') <- attr(tagMatrix[[1]][[1]], 'is.binning')
    binFlag <- attr(tagMatrix[[1]][[1]], 'is.binning')
    type <- attr(tagMatrix[[1]][[1]], 'type')
    
  }else{
    upstream <- attr(tagMatrix[[1]], 'upstream')
    downstream <- attr(tagMatrix[[1]], 'downstream')
    binFlag <- attr(tagMatrix[[1]], 'is.binning')
    type <- attr(tagMatrix[[1]], 'type')
  }
  
  if(type == "body"){
    
    label <- c("SS","TS")
    
  }else if(type == "start_site"){
    
    label <- "SS"
    
  }else{
    
    label <- "TS"
    
  }
  
  
  if(binFlag){
    
    if (!(missingArg(conf) || is.na(conf))){
      
      plotMultiProf.binning(tagMatrix = tagMatrix, 
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
      
      plotMultiProf.binning(tagMatrix = tagMatrix, 
                            xlab = xlab,
                            ylab = ylab,
                            facet = facet, 
                            free_y = free_y,
                            upstream = upstream,
                            downstream = downstream,
                            label = label,
                            ...)
    }
    
    
  }else{
    
    xlim <- c(-upstream, downstream)
    
    if (!(missingArg(conf) || is.na(conf))){
      
      plotMultiProf.normal(tagMatrix = tagMatrix, 
                           xlim = xlim,
                           xlab = xlab,
                           ylab = ylab,
                           conf = conf,
                           facet = facet, 
                           free_y = free_y,
                           origin_label = label,
                           ...)
      
    }else{
      
      plotMultiProf.normal(tagMatrix = tagMatrix, 
                           xlim = xlim,
                           xlab = xlab,
                           ylab = ylab,
                           facet = facet, 
                           free_y = free_y,
                           origin_label = label,
                           ...)
    }
  }
  
}

##' internal function
##' 
##' @param tagMatrix tagMatrix
##' @param xlim xlim
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param origin_label the label of the center
##' @param verbose print message or not
##' @param ... additional parameter
plotMultiProf.normal <- function(tagMatrix, xlim,
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
    
    p <- plotMultiProf.normal.internal(tagMatrix = tagMatrix, 
                                       conf = conf, 
                                       xlim = xlim,
                                       xlab = xlab, 
                                       ylab = ylab,
                                       facet = facet, 
                                       free_y = free_y, 
                                       origin_label = origin_label,
                                       ...)
    
    
  } else {
    
    p <- plotMultiProf.normal.internal(tagMatrix, 
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

##' internal function
##' 
##' 
##' @param tagMatrix tagMatrix
##' @param xlim xlim
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param origin_label the label of the center
##' @param ... additional parameter
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
plotMultiProf.normal.internal <- function(tagMatrix, conf,
                                          xlim = c(-3000,3000),
                                          xlab = "Genomic Region (5'->3')",
                                          ylab = "Peak Count Frequency",
                                          facet="row", 
                                          free_y = TRUE,
                                          origin_label, 
                                          ...) {
  
  listFlag <- FALSE
  if (is.null(attr(tagMatrix[[1]],'upstream'))) {
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
    if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix[[1]][[1]]) ) {
      stop("please specify appropreate xcoordinations...")
    }
  } else {
    if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix[[1]]) ) {
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
    
    tagCount <- lapply(as.list(names(tagMatrix)), function(x){
      
      tmp <- tagMatrix[[x]]
      tagCount_tmp <- lapply(as.list(names(tmp)),function(x){
        result <- getTagCount(tmp[[x]], xlim = xlim, conf = conf, ...)
        result$type <- x
        
        return(result)
      })
      tagCount_tmp <- list_to_dataframe(tagCount_tmp)
      return(tagCount_tmp)
      
    })
    
    names(tagCount) <- names(tagMatrix)
    tagCount <- list_to_dataframe(tagCount)
    tagCount$.id <- factor(tagCount$.id, levels=names(tagMatrix))
    p <- ggplot(tagCount, aes(pos, group=type, color=type))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = type),
                           linetype = 0, alpha = 0.2)
    }
    
  } else {
    
    tagCount <- lapply(as.list(names(tagMatrix)), function(x){
      
      result <- getTagCount(tagMatrix[[x]], xlim = xlim, conf = conf, ...)
      result$type <- x
      
      return(result)
    })
    
    tagCount <- do.call("rbind",tagCount)
    
    p <- ggplot(tagCount, aes(x = pos))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper,fill = type),
                           linetype = 0, alpha = 0.2)
    }
  }
  
  p <- p + geom_line(aes(y = value,color = type))
  
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
    # cols <- getCols(length(tagMatrix[[1]]))
    # p <- p + scale_color_manual(values=cols)
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

  # if(facet != "none") {
  #   p <- p + theme(legend.position="none")
  # }
  
  return(p)
}

##' internal function
##' 
##' @param tagMatrix tagMatrix
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param upstream the upstream extension
##' @param downstream the downstream extension
##' @param label the label of the center
##' @param ... additional parameter
plotMultiProf.binning <- function(tagMatrix, 
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
    p <- plotMultiProf.binning.internal(tagMatrix , 
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
    p <- plotMultiProf.binning.internal(tagMatrix , 
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

##' internal function
##' 
##' @param tagMatrix tagMatrix
##' @param xlab xlab
##' @param ylab ylab
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param upstream the upstream extension
##' @param downstream the downstream extension
##' @param label the label of the center
##' @param ... additional parameter
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
plotMultiProf.binning.internal <- function(tagMatrix, 
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
  if (is(tagMatrix[[1]][[1]],"matrix")) {
    if ( is.null(names(tagMatrix)) ) {
      nn <- paste0("peak", seq_along(tagMatrix))
      warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
      names(tagMatrix) <- nn
      ## stop("tagMatrix should be a named list...")
    }
    listFlag <- TRUE
  }
  
  if(listFlag){
    nbin <- dim(tagMatrix[[1]][[1]])[2]
    type <- attr(tagMatrix[[1]][[1]], 'type')
  }else{
    nbin <- dim(tagMatrix[[1]])[2]
    type <- attr(tagMatrix[[1]], 'type')
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
    
    tagCount <- lapply(as.list(names(tagMatrix)), function(x){
      
      tmp <- tagMatrix[[x]]
      tagCount_tmp <- lapply(as.list(names(tmp)),function(x){
        result <- getTagCount(tmp[[x]], xlim = xlim, conf = conf, ...)
        result$type <- x
        
        return(result)
      })
      tagCount_tmp <- list_to_dataframe(tagCount_tmp)
      return(tagCount_tmp)
      
    })
    
    names(tagCount) <- names(tagMatrix)
    tagCount <- list_to_dataframe(tagCount)
    tagCount$.id <- factor(tagCount$.id, levels=names(tagMatrix))
    p <- ggplot(tagCount, aes(pos, group=type, color=type))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = type),
                           linetype = 0, alpha = 0.2)
    }
    
  } else {
    
    tagCount <- lapply(as.list(names(tagMatrix)), function(x){
      
      result <- getTagCount(tagMatrix[[x]], xlim = xlim, conf = conf, ...)
      result$type <- x
      
      return(result)
    })
    
    tagCount <- do.call("rbind",tagCount)
    
    p <- ggplot(tagCount, aes(pos,group=type,color=type))
    if (!(is.na(conf))) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper,fill = type),
                           linetype = 0, alpha = 0.2)
    }
  }
  
  p <- p + geom_line(aes(y = value,color = type))
  
  ## x_scale for genebody
  if(type == 'body'){
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
  if(type != 'body'){
    
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
  # if(facet != "none") {
  #   p <- p + theme(legend.position="none")
  # }
  return(p)
}


##' plot the heatmap of tagMatrix
##'
##'
##' @title tagHeatmap
##' @param tagMatrix tagMatrix or a list of tagMatrix
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @param palette palette to be filled in,details see \link[ggplot2]{scale_colour_brewer}
##' @param nrow the nrow of plotting a list of peak
##' @param ncol the ncol of plotting a list of peak
##' @return figure
##' @export
##' @author G Yu
tagHeatmap <- function(tagMatrix, 
                       xlab="", 
                       ylab="", 
                       title=NULL, 
                       palette="RdBu",
                       nrow = NULL,
                       ncol = NULL) {
  listFlag <- FALSE
  if (is(tagMatrix, "list")) {
    listFlag <- TRUE
  }
  peakHeatmap.internal2(tagMatrix = tagMatrix, 
                        listFlag = listFlag, 
                        palette = palette, 
                        xlab = xlab, 
                        ylab = ylab, 
                        title = title,
                        ncol = ncol,
                        nrow = nrow)
}

##' plot the heatmap of peaks 
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
##' @param palette palette to be filled in,details see \link[ggplot2]{scale_colour_brewer}
##' @param verbose print message or not
##' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR', 'UTR'
##' @param type one of "start_site", "end_site", "body"
##' @param nbin the amount of nbines 
##' @param ignore_strand ignore the strand information or not
##' @param windows a collection of region
##' @param nrow the nrow of plotting a list of peak
##' @param ncol the ncol of plotting a list of peak
##' @return figure
##' @export
##' @author G Yu
peakHeatmap <- function(peak, weightCol=NULL, TxDb=NULL,
                        upstream=1000, downstream=1000,
                        xlab="", ylab="", title=NULL,
                        palette=NULL, verbose=TRUE,
                        by="gene", type="start_site",
                        nbin = NULL,ignore_strand = FALSE,
                        windows,ncol = NULL, nrow = NULL) {
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
  
  if (verbose) {
    cat(">> preparing tag matrix...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  
  if(missing(windows)){
    windows <- getBioRegion(TxDb=TxDb,
                            upstream=upstream,
                            downstream=downstream,
                            by=by,
                            type=type)
  }
  
  
  if (listFlag) {
    tagMatrix <- lapply(peak, getTagMatrix, 
                        weightCol=weightCol, 
                        windows = windows,
                        upstream=upstream,
                        downstream=downstream,
                        TxDb = TxDb,
                        nbin = nbin,
                        verbose = verbose,
                        ignore_strand= ignore_strand)
    
    names(tagMatrix) <- names(peak)
    
  } else {
    tagMatrix <- getTagMatrix(peak, 
                              weightCol=weightCol, 
                              windows = windows,
                              TxDb = TxDb,
                              upstream=upstream,
                              downstream=downstream,
                              nbin = nbin,
                              verbose = verbose,
                              ignore_strand= ignore_strand)
  }
  
  if (verbose) {
    cat(">> generating figure...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  
  xlim <- NULL
  
  p <- peakHeatmap.internal2(tagMatrix = tagMatrix,
                             listFlag = listFlag, 
                             palette = palette, 
                             xlab = xlab,
                             ylab = ylab, 
                             title = title,
                             nrow = nrow,
                             ncol = ncol)
  
  if (verbose) {
    cat(">> done...\t\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  invisible(tagMatrix)
  p
}

##' @importFrom aplot plot_list
peakHeatmap.internal2 <- function(tagMatrix, 
                                  listFlag, 
                                  palette, 
                                  xlab, 
                                  ylab, 
                                  title,
                                  nrow,
                                  ncol) {
  if ( is.null(xlab) || is.na(xlab))
    xlab <- ""
  if ( is.null(ylab) || is.na(ylab))
    ylab <- ""
  
  if (listFlag) {
    nc <- length(tagMatrix)
    if ( is.null(palette) || is.na(palette) ) {
      palette <- getPalette(nc)
    } else if (length(palette) != nc) {
      palette <- rep(palette[1], nc)
    } else {
      palette <- palette
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
    
    tmp <- list()
    
    for (i in 1:nc) {
      
      p <- peakHeatmap.internal(tagMatrix = tagMatrix[[i]], 
                                palette = palette[i], 
                                xlab = xlab[i], 
                                ylab = ylab[i], 
                                title= title[i])
      
      p <- p + theme(plot.title = element_text(hjust = 0.5))
      
      tmp[[i]] <- p
    }
    
    if(is.null(nrow) && is.null(ncol))
      nrow <- 1
    
    p <- plot_list(gglist = tmp,
                   ncol = ncol,
                   nrow = nrow)
    return(p)
    
  } else {
    if (is.null(palette) || is.na(palette))
      palette <- "RdBu"
    if (is.null(title) || is.na(title))
      title <- ""
    peakHeatmap.internal(tagMatrix = tagMatrix, 
                         palette = palette, 
                         xlab = xlab, 
                         ylab = ylab, 
                         title = title)
  }
}


##' @import BiocGenerics
##' @importFrom tibble rownames_to_column
##' @importFrom tidyr pivot_longer
##' @importFrom magrittr %>%
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 scale_fill_distiller
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 scale_x_continuous
peakHeatmap.internal <- function(tagMatrix, 
                                 palette="RdBu", 
                                 xlab="", 
                                 ylab="",
                                 title="") {
  
  upstream <- attr(tagMatrix, "upstream")
  downstream <- attr(tagMatrix, "downstream")
  binning_Flag <- attr(tagMatrix,"is.binning")
  type <- attr(tagMatrix,"type")
  
  body_Flag <- FALSE
  if(type == "body"){
    body_Flag <- TRUE
    label <- attr(tagMatrix,"label")
  }
  
  if(binning_Flag){
    nbin <- dim(tagMatrix)[2]
  }
  
  tagMatrix <- t(apply(tagMatrix, 1, function(x) x/max(x)))
  ii <- order(rowSums(tagMatrix))
  tagMatrix <- tagMatrix[ii,]
  
  colnames(tagMatrix) <- seq_len(dim(tagMatrix)[2])
  rownames(tagMatrix) <- seq_len(dim(tagMatrix)[1])
  tagMatrix <- tagMatrix %>% as.data.frame() %>% 
    rownames_to_column("sample_ID") %>%
    pivot_longer(-c(sample_ID),names_to = "coordinate", 
                 values_to = "values")
  tagMatrix$coordinate <- as.numeric(tagMatrix$coordinate)

  sample_ID <- coordinate <- NULL
  
  p <- ggplot(tagMatrix, aes(x = coordinate,y = sample_ID)) + 
    geom_tile(aes(fill = values)) +
    scale_fill_distiller(palette = palette)  +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          panel.grid=element_blank(),
          panel.background = element_blank()) +
    labs(x = xlab, y = ylab, title = title)

  if(body_Flag){
    
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
    }
    
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
    
    if(!is.null(upstream) && !inherits(upstream, 'rel')){
      
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
    }
    
    return(p)
    
  }
  
  if(binning_Flag){
    
    p <- p + scale_x_continuous(breaks = c(1,
                                           floor(nbin*(downstream*0.5/(downstream+upstream))),
                                           floor(nbin*(downstream/(downstream+upstream))),
                                           floor(nbin*((downstream + upstream*0.5)/(downstream+upstream))),
                                           nbin),
                                labels = c((-1*downstream),
                                           floor(-1*downstream*0.5),
                                           0,
                                           floor(upstream*0.5),
                                           upstream))
  }else{
    
    p <- p + scale_x_continuous(breaks = c(1,
                                           floor(downstream*0.5),
                                           (downstream + 1),
                                           (downstream + 1 + floor(upstream * 0.5)), 
                                           upstream+downstream+1),
                                labels = c((-1*downstream),
                                           floor(-1*downstream*0.5),
                                           0,
                                           floor(upstream*0.5),
                                           upstream))    
    
  }
  
  p
}

##' plot the heatmap of peaks align to a sets of regions
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
##' @param palette palette to be filled in,details see \link[ggplot2]{scale_colour_brewer}
##' @param verbose print message or not
##' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR', 'UTR'
##' @param type one of "start_site", "end_site", "body"
##' @param nbin the amount of nbines 
##' @param ignore_strand ignore the strand information or not
##' @param windows_name the name for each window, which will also be showed in the picture as labels
##' @param nrow the nrow of plotting a list of peak
##' @param ncol the ncol of plotting a list of peak
##' @param facet_label_text_size the size of facet label text
##' @importFrom tibble rownames_to_column
##' @importFrom tidyr pivot_longer
##' @importFrom magrittr %>%
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 scale_fill_distiller
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 scale_x_continuous
##' @return figure
##' @export
peakHeatmap_multiple_Sets <- function(peak, 
                                      weightCol=NULL,
                                      TxDb=NULL,
                                      upstream=1000, 
                                      downstream=1000,
                                      xlab="", 
                                      ylab="", 
                                      title=NULL,
                                      palette=NULL, 
                                      verbose=TRUE,
                                      by="gene", 
                                      type="start_site",
                                      nbin = NULL,
                                      ignore_strand = FALSE,
                                      windows_name = NULL,
                                      ncol = NULL,
                                      nrow = NULL,
                                      facet_label_text_size = 12){
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
  
  
  ## check type
  if(length(type) != 1){
    stop("It should be only one type...")
  }
  
  if(is.null(windows_name) && !is.null(names(TxDb)))
    windows_name <- names(TxDb)
  
  ## make the window name
  if (is.null(windows_name)) {
    nn <- by
    warning("set the name automatically to ", paste(nn, collapse=' '))
    windows_name <- nn
  }else{
    if (length(windows_name) != length(by)) {
      stop("the length of the window name and the by should be equal...")
    }
  }
  
  if ( is(peak, "list") ) {
    tagMatrix <- lapply(peak, getTagMatrix2,
                        upstream=upstream,
                        downstream=downstream,
                        windows_name=windows_name,
                        type=type,
                        by=by,
                        TxDb=TxDb,
                        weightCol = weightCol, 
                        nbin = nbin,
                        verbose = verbose,
                        ignore_strand= ignore_strand)
  } else {
    tagMatrix <- getTagMatrix2(peak=peak, 
                               upstream=upstream,
                               downstream=downstream,
                               windows_name=windows_name,
                               type=type,
                               by=by,
                               TxDb=TxDb,
                               weightCol = weightCol, 
                               nbin = nbin,
                               verbose = verbose,
                               ignore_strand= ignore_strand)
  }
  
  if(listFlag){
    
    nc <- length(tagMatrix)
    if ( is.null(palette) || is.na(palette) ) {
      palette <- getPalette(nc)
    } else if (length(palette) != nc) {
      palette <- rep(palette[1], nc)
    } else {
      palette <- palette
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
    
    tmp <- list()
    
    for (i in 1:nc) {
      
      p <- peakHeatmap_multiple_Sets.internal(tagMatrix = tagMatrix[[i]],
                                              upstream=upstream, 
                                              downstream=downstream,
                                              xlab=xlab[[i]], 
                                              ylab=ylab[[i]], 
                                              title=title[[i]],
                                              palette=palette[[i]], 
                                              ncol = ncol,
                                              nrow = nrow,
                                              facet_label_text_size = facet_label_text_size)
      
      p <- p + theme(plot.title = element_text(hjust = 0.5))
      
      tmp[[i]] <- p
    }
    
    if(is.null(nrow) && is.null(ncol))
      nrow <- 1
    
    p <- plot_list(gglist = tmp,
                   ncol = ncol,
                   nrow = nrow)
    
  }else{
    
    if (is.null(palette) || is.na(palette))
      palette <- "RdBu"
    if (is.null(title) || is.na(title))
      title <- ""
    
    p <- peakHeatmap_multiple_Sets.internal(tagMatrix = tagMatrix,
                                            upstream=upstream, 
                                            downstream=downstream,
                                            xlab=xlab, 
                                            ylab=ylab, 
                                            title=title,
                                            palette=palette, 
                                            ncol = ncol,
                                            nrow = nrow,
                                            facet_label_text_size = facet_label_text_size)
    
  }
  
  return(p)
  
}


##' @importFrom tibble rownames_to_column
##' @importFrom tidyr pivot_longer
##' @importFrom magrittr %>%
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 scale_fill_distiller
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 facet_grid
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 element_blank
peakHeatmap_multiple_Sets.internal <- function(tagMatrix,
                                               upstream=1000, 
                                               downstream=1000,
                                               xlab="", 
                                               ylab="", 
                                               title=NULL,
                                               palette=NULL, 
                                               ncol = NULL,
                                               nrow = NULL,
                                               facet_label_text_size = 12){

  binning_Flag <- attr(tagMatrix[[1]],"is.binning")
  if(binning_Flag) nbin <- dim(tagMatrix[[1]])[2]
  
  type <- attr(tagMatrix,"type")
  body_Flag <- FALSE
  if(attr(tagMatrix[[1]],"type") == "body"){
    body_Flag <- TRUE
    label <- attr(tagMatrix,"label")
  }
  
  name_of_list <- as.list(names(tagMatrix))
  
  peak_list <- lapply(name_of_list,function(x){
    
    tagMatrix[[x]] <- t(apply(tagMatrix[[x]], 1, function(x) x/max(x)))
    ii <- order(rowSums(tagMatrix[[x]]))
    tagMatrix[[x]] <- tagMatrix[[x]][ii,]
    
    colnames(tagMatrix[[x]]) <- seq_len(dim(tagMatrix[[x]])[2])
    rownames(tagMatrix[[x]]) <- seq_len(dim(tagMatrix[[x]])[1])
    tagMatrix[[x]] <- tagMatrix[[x]] %>% as.data.frame() %>% 
      rownames_to_column("sample_ID") %>%
      pivot_longer(-c(sample_ID),names_to = "coordinate", 
                   values_to = "values")
    tagMatrix[[x]]$coordinate <- as.numeric(tagMatrix[[x]]$coordinate)
    tagMatrix[[x]]$sample <- x
    return(tagMatrix[[x]])
  })
  
  peak_df <- list_to_dataframe(peak_list)
  
  sample_ID <- coordinate <- NULL
  
  p <- ggplot(peak_df, aes(x = coordinate,y = sample_ID)) + 
    geom_tile(aes(fill = values)) +
    scale_fill_distiller(palette = palette)  +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          panel.grid=element_blank(),
          panel.background = element_blank()) +
    labs(x = xlab, y = ylab, title = title)
  
  if(body_Flag){
    
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
    }
    
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
    
    if(!is.null(upstream) && !inherits(upstream, 'rel')){
      
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
    }
    
    p <-  p + facet_grid(sample ~ .,switch = "y",space = "free_y",scales = "free_y") +
      theme(strip.text.y.left = element_text(color = "black",face = "bold",
                                             size = facet_label_text_size),
            strip.background = element_blank())
    
    return(p)
    
  }
  
  if(binning_Flag){
    
    p <- p + scale_x_continuous(breaks = c(1,
                                           floor(nbin*(downstream*0.5/(downstream+upstream))),
                                           floor(nbin*(downstream/(downstream+upstream))),
                                           floor(nbin*((downstream + upstream*0.5)/(downstream+upstream))),
                                           nbin),
                                labels = c((-1*downstream),
                                           floor(-1*downstream*0.5),
                                           0,
                                           floor(upstream*0.5),
                                           upstream))
  }else{
    
    p <- p + scale_x_continuous(breaks = c(1,
                                           floor(downstream*0.5),
                                           (downstream + 1),
                                           (downstream + 1 + floor(upstream * 0.5)), 
                                           upstream+downstream+1),
                                labels = c((-1*downstream),
                                           floor(-1*downstream*0.5),
                                           0,
                                           floor(upstream*0.5),
                                           upstream))    
    
  }
  
  p <-  p + facet_grid(sample ~ .,switch = "y",scales = "free_y",space = "free") +
    theme(strip.text.y.left = element_text(color = "black",face = "bold",
                                           size = facet_label_text_size),
          strip.background = element_blank())
  
  return(p)
  
}




##' plot peak heatmap and profile in a picture
##' 
##' 
##' @title peak_Profile_Heatmap
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight
##' @param TxDb TxDb object
##' @param upstream upstream position
##' @param downstream downstream position
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @param palette palette to be filled in,details see \link[ggplot2]{scale_colour_brewer}
##' @param verbose print message or not
##' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR', 'UTR'
##' @param type one of "start_site", "end_site", "body"
##' @param nbin the amount of nbines 
##' @param ignore_strand ignore the strand information or not
##' @param windows_name the name for each window, which will also be showed in the picture as labels
##' @param nrow the nrow of plotting a list of peak
##' @param ncol the ncol of plotting a list of peak
##' @param facet_label_text_size the size of facet label text
##' @param conf confidence interval
##' @param facet one of 'none', 'row' and 'column'
##' @param free_y if TRUE, y will be scaled by AvgProf
##' @param height_proportion the proportion of profiling picture and heatmap
##' @importFrom aplot insert_bottom
##' @importFrom aplot plot_list
##' @export
peak_Profile_Heatmap <- function(peak, 
                                 weightCol=NULL,
                                 TxDb=NULL,
                                 upstream=1000, 
                                 downstream=1000,
                                 xlab="", 
                                 ylab="", 
                                 title=NULL,
                                 palette=NULL, 
                                 verbose=TRUE,
                                 by="gene", 
                                 type="start_site",
                                 nbin = NULL,
                                 ignore_strand = FALSE,
                                 windows_name = NULL,
                                 ncol = NULL,
                                 nrow = NULL,
                                 facet_label_text_size = 12,
                                 conf,
                                 facet = "row",
                                 free_y = TRUE,
                                 height_proportion = 4){
  
  conf <- if(missingArg(conf)) NA else conf
  
  if(is(peak, "list")){
    
    nc <- length(peak)
    
    tmp <- list()
    
    if ( is.null(names(peak)) ) {
      nn <- paste0("peak", seq_along(peak))
      warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
      names(peak) <- nn
      ## stop("tagMatrix should be a named list...")
    }
    
    if(is.null(palette)) palette <- getPalette(nc)
    
    if(is.null(title)) title_of_plot <- names(peak)
    
    for (i in 1:nc) {
      peak_profile <- plotPeakProf(peak = peak[[i]],
                                   upstream = upstream,
                                   downstream = downstream,
                                   conf = conf,
                                   by = by,
                                   type = type,
                                   windows_name = windows_name,
                                   weightCol = weightCol,
                                   TxDb = TxDb,
                                   xlab = xlab,
                                   ylab = ylab,
                                   facet = facet,
                                   free_y = free_y,
                                   verbose = verbose,
                                   nbin = nbin,
                                   ignore_strand = ignore_strand)
      
      peak_profile <- peak_profile + labs(title = title_of_plot[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      
      if(length(by) != 1){
        peak_heatmap <- peakHeatmap_multiple_Sets(peak = peak[[i]], 
                                                  weightCol=weightCol,
                                                  TxDb=TxDb,
                                                  upstream=upstream, 
                                                  downstream=downstream,
                                                  xlab=xlab, 
                                                  ylab=ylab, 
                                                  title=title,
                                                  palette=palette[[i]], 
                                                  verbose=verbose,
                                                  by=by, 
                                                  type=type,
                                                  nbin = nbin,
                                                  ignore_strand = ignore_strand,
                                                  windows_name = windows_name,
                                                  ncol = ncol,
                                                  nrow = nrow,
                                                  facet_label_text_size = facet_label_text_size)
      }else{
        
        peak_heatmap <- peakHeatmap(peak[[i]], 
                                    weightCol=weightCol, 
                                    TxDb=TxDb,
                                    upstream=upstream, 
                                    downstream=downstream,
                                    xlab=xlab, 
                                    ylab=ylab, 
                                    title=title,
                                    palette=palette[[i]], 
                                    verbose=verbose,
                                    by=by, 
                                    type=type,
                                    nbin = nbin,
                                    ignore_strand = ignore_strand,
                                    ncol = ncol,
                                    nrow = nrow)
        
      }
      
      p <- peak_profile %>% 
        insert_bottom(peak_heatmap,height = height_proportion)
      
      tmp[[i]] <- p
    }
    
    if (is.null(ncol) && is.null(nrow))
      nrow <- 1
    
    p <- plot_list(gglist = tmp,
                   ncol = ncol,
                   nrow = nrow)
    
    return(p)
    
  }
  
  peak_profile <- plotPeakProf(peak = peak,
                               upstream = upstream,
                               downstream = downstream,
                               conf = conf,
                               by = by,
                               type = type,
                               windows_name = windows_name,
                               weightCol = weightCol,
                               TxDb = TxDb,
                               xlab = xlab,
                               ylab = ylab,
                               facet = facet,
                               free_y = free_y,
                               verbose = verbose,
                               nbin = nbin,
                               ignore_strand = ignore_strand)
  
  
  if(length(by) != 1){
    peak_heatmap <- peakHeatmap_multiple_Sets(peak = peak, 
                                              weightCol=weightCol,
                                              TxDb=TxDb,
                                              upstream=upstream, 
                                              downstream=downstream,
                                              xlab=xlab, 
                                              ylab=ylab, 
                                              title=title,
                                              palette=palette, 
                                              verbose=verbose,
                                              by=by, 
                                              type=type,
                                              nbin = nbin,
                                              ignore_strand = ignore_strand,
                                              windows_name = windows_name,
                                              ncol = ncol,
                                              nrow = nrow,
                                              facet_label_text_size = facet_label_text_size)
  }else{
    
    peak_heatmap <- peakHeatmap(peak = peak, 
                                weightCol=weightCol, 
                                TxDb=TxDb,
                                upstream=upstream, 
                                downstream=downstream,
                                xlab=xlab, 
                                ylab=ylab, 
                                title=title,
                                palette=palette, 
                                verbose=verbose,
                                by=by, 
                                type=type,
                                nbin = nbin,
                                ignore_strand = ignore_strand,
                                ncol = ncol,
                                nrow = nrow)
    
  }
  
  p <- peak_profile %>% 
    insert_bottom(peak_heatmap,height = height_proportion)
  

  return(p)
}