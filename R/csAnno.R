##' Class "csAnno"
##' This class represents the output of ChIPseeker Annotation
##'
##'
##' @name csAnno-class
##' @aliases csAnno-class
##' show,csAnno-method vennpie,csAnno-method
##' plotDistToTSS,csAnno-method plotAnnoBar,csAnno-method
##' plotAnnoPie,csAnno-method upsetplot,csAnno-method
##' subset,csAnno-method
##'
##' @docType class
##' @slot anno annotation
##' @slot tssRegion TSS region
##' @slot level transcript or gene
##' @slot hasGenomicAnnotation logical
##' @slot detailGenomicAnnotation Genomic Annotation in detail
##' @slot annoStat annotation statistics
##' @slot peakNum number of peaks
##' @exportClass csAnno
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @seealso \code{\link{annotatePeak}}
##' @keywords classes
setClass("csAnno",
         representation=representation(
             anno = "GRanges",
             tssRegion = "numeric",
             level = "character",
             hasGenomicAnnotation = "logical",
             detailGenomicAnnotation="data.frame",
             annoStat="data.frame",
             peakNum="numeric"
             ))


##' convert csAnno object to GRanges
##'
##'
##' @title as.GRanges
##' @param x csAnno object
##' @return GRanges object
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @export
as.GRanges <- function(x) {
    if (!is(x, "csAnno"))
        stop("not supported...")
    return(x@anno)
}

##' getting status of annotation
##' 
##' 
##' @title getAnnoStat
##' @param x csAnno object
##' @export
getAnnoStat <- function(x) {
    if (!is(x, "csAnno"))
        stop("not supported...")
    return(x@annoStat)
}



##' Combine csAnno Object
##'
##'
##' https://github.com/YuLab-SMU/ChIPseeker/issues/157
##' @title combine_csAnno
##' @param x csAnno object
##' @param ... csAnno objects
##' @return csAnno object
##' @export
combine_csAnno <- function(x, ...){
    z <- list(x, ...)
    
    if(sum(vapply(z, function(x) !is(x, "csAnno"), FUN.VALUE = logical(1))) != 0){
        stop("not supported...")
    }
    
    if(length(z)<2){
        stop("need two or more csAnno object...")
    }
    
    
    if(sum(!duplicated(lapply(z, function(x) x@tssRegion[1]))) != 1 
       && sum(!duplicated(lapply(z, function(x) x@tssRegion[2]))) != 1){
        stop("the tss regions of different csAnno objects should be the same...")
    }
    
    if(sum(!duplicated(lapply(z, function(x) x@level))) != 1){
        stop("the level of different csAnno object should be the same...")
    }
    
    if(sum(!duplicated(lapply(z, function(x) x@hasGenomicAnnotation))) != 1){
        stop("the status of GenomicAnnotation should be the same...")
    }
    
    combine_tssRegion <- x@tssRegion
    combine_level <- x@level
    combine_hasGenomicAnnotation <- x@hasGenomicAnnotation
    
    combine_anno <- x@anno
    for(i in 2:length(z)){
        combine_anno <- c(combine_anno,z[[i]]@anno)
    }
    
    combine_detailGenomicAnnotation <- lapply(z, function(x) x@detailGenomicAnnotation)
    combine_detailGenomicAnnotation <- do.call("rbind",combine_detailGenomicAnnotation)
    
    combine_peakNum <- x@peakNum
    for(i in 2:length(z)){
        combine_peakNum <- combine_peakNum+z[[i]]@peakNum
    }
    
    feature <- x@annoStat$Feature
    for(i in 2:length(z)){
        if(length(feature)<length(z[[i]]@annoStat$Feature)){
            feature_levels <- levels(z[[i]]@annoStat$Feature)
            feature <- c(as.vector(feature),as.vector(z[[i]]@annoStat$Feature))
            feature <- feature[!duplicated(feature)]
            feature <- factor(feature, 
                              levels = feature_levels)
            feature <- sort(feature)
        }else{
            feature_levels <- levels(feature)
            feature <- c(as.vector(feature),as.vector(z[[i]]@annoStat$Feature))
            feature <- feature[!duplicated(feature)]
            feature <- factor(feature, 
                              levels = feature_levels)
            feature <- sort(feature)
        }
    }
    
    combine_annoStat <- data.frame(Feature=feature)
    
    for(i in 1:length(z)){
        combine_annoStat <- merge(combine_annoStat, z[[i]]@annoStat, 
                                  by = "Feature", all = T, sort = F)
        combine_annoStat[is.na(combine_annoStat)] <- 0
        combine_annoStat <- combine_annoStat[order(combine_annoStat$Feature),]
    }
    
    total <- (ncol(combine_annoStat)-1)*100
    combine_annoStat$sum <- rowSums(combine_annoStat[, 2:ncol(combine_annoStat)])
    
    
    for (i in 1:length(combine_annoStat$sum)) {
        combine_annoStat$result[i] <- (combine_annoStat$sum[i]/total)*100
    }
    
    annoStat_result <- data.frame(Feature=combine_annoStat[,1],Frequency=combine_annoStat[,ncol(combine_annoStat)])
    
    res <- new("csAnno",
               anno = combine_anno,
               tssRegion = combine_tssRegion,
               level = combine_level,
               hasGenomicAnnotation = combine_hasGenomicAnnotation,
               detailGenomicAnnotation = combine_detailGenomicAnnotation,
               annoStat = annoStat_result,
               peakNum = combine_peakNum
    )
    
    return(res)
}

##' vennpie method generics
##'
##' @name vennpie
##' @docType methods
##' @rdname vennpie-methods
##'
##' @title vennpie method
##' @param x A \code{csAnno} instance
##' @param r initial radius
##' @param ... additional parameter
##' @return plot
##' @usage vennpie(x, r=0.2, ...)
##' @exportMethod vennpie
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("vennpie", signature(x="csAnno"),
          function(x, r=0.2, ...) {
              vennpie.csAnno(x, r, ...)
          }
          )


##' upsetplot method generics
##'
##' @name upsetplot
##' @docType methods
##' @rdname upsetplot-methods
##'
##' @title upsetplot method
##' @param x A \code{csAnno} instance
##' @param ... additional parameter
##' @return plot
##' @usage upsetplot(x, ...)
##' @importFrom enrichplot upsetplot
##' @exportMethod upsetplot
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("upsetplot", signature(x="csAnno"),
          function(x, ...) {
              upsetplot.csAnno(x, ...)
          }
          )

##' convert csAnno object to data.frame
##'
##'
##' @title as.data.frame.csAnno
##' @param x csAnno object
##' @param row.names row names
##' @param optional should be omitted.
##' @param ... additional parameters
##' @return data.frame
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @method as.data.frame csAnno
##' @export
as.data.frame.csAnno <- function(x, row.names=NULL, optional=FALSE, ...) {
    y <- as.GRanges(x)
    if (!(is.null(row.names) || is.character(row.names)))
        stop("'row.names' must be NULL or a character vector")
    df <- as.data.frame(y)
    rownames(df) <- row.names
    return(df)
}

##' show method for \code{csAnno} instance
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##' @aliases show,csAnno,ANY-method
##' @title show method
##' @param object A \code{csAnno} instance
##' @return message
##' @importFrom methods show
##' @exportMethod show
##' @usage show(object)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("show", signature(object="csAnno"),
          function(object) {
              cat("Annotated peaks generated by ChIPseeker\n")
              cat(paste(length(object@anno), object@peakNum, sep="/"),
                  " peaks were annotated\n")
              if (object@hasGenomicAnnotation) {
                  cat("Genomic Annotation Summary:\n")
                  print(object@annoStat)
              }
          }
          )

##' plotAnnoBar method for list of \code{csAnno} instances
##'
##' @name plotAnnoBar
##' @docType methods
##' @rdname plotAnnoBar-methods
##' @aliases plotAnnoBar,list-method
##' @exportMethod plotAnnoBar
setMethod("plotAnnoBar", signature(x="list"),
          function(x,
                   xlab="",
                   ylab='Percentage(%)',
                   title="Feature Distribution",
                   ...) {
              if (is.null(names(x))) {
                  nn <- paste0("Peak", seq_along(x))
                  warning("input is not a named list, set the name automatically to ", paste(nn, collapse = " "))
                  names(x) <- nn
                  ## stop("input object should be a named list...")
              }
              anno <- lapply(x, getAnnoStat)
              ## anno.df <- ldply(anno)
              anno.df <- list_to_dataframe(anno)
              categoryColumn <- ".id"
              plotAnnoBar.data.frame(anno.df, xlab, ylab, title, categoryColumn)
          })

##' plotAnnoBar method for \code{csAnno} instance
##'
##' @name plotAnnoBar
##' @docType methods
##' @rdname plotAnnoBar-methods
##' @aliases plotAnnoBar,csAnno,ANY-method
##' @title plotAnnoBar method
##' @param x \code{csAnno} instance
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @param ... additional paramter
##' @return plot
##' @exportMethod plotAnnoBar
##' @usage plotAnnoBar(x, xlab="", ylab='Percentage(\%)',title="Feature Distribution", ...)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("plotAnnoBar", signature(x="csAnno"),
          function(x,
                   xlab="",
                   ylab="Percentage(%)",
                   title="Feature Distribution",
                   ...) {
              anno.df <- getAnnoStat(x)
              categoryColumn <- 1
              plotAnnoBar.data.frame(anno.df, xlab, ylab, title, categoryColumn)
          })



##' plotAnnoPie method for \code{csAnno} instance
##'
##' @name plotAnnoPie
##' @docType methods
##' @rdname plotAnnoPie-methods
##' @aliases plotAnnoPie,csAnno,ANY-method
##' @title plotAnnoPie method
##' @param x \code{csAnno} instance
##' @param ndigit number of digit to round
##' @param cex label cex
##' @param col color
##' @param legend.position topright or other.
##' @param pie3D plot in 3D or not
##' @param radius radius of the pie
##' @param ... extra parameter
##' @return plot
##' @exportMethod plotAnnoPie
##' @usage plotAnnoPie(x,ndigit=2,cex=0.9,col=NA,legend.position="rightside",pie3D=FALSE,radius=0.8,...)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("plotAnnoPie", signature(x="csAnno"),
          function(x,
                   ndigit=2,
                   cex=0.9,
                   col=NA,
                   legend.position="rightside",
                   pie3D=FALSE,
                   radius=0.8,
                   ...){
              plotAnnoPie.csAnno(x, ndigit, cex, col, legend.position, pie3D, radius, ...)
          })



##' plotDistToTSS method for list of \code{csAnno} instances
##'
##' @name plotDistToTSS
##' @docType methods
##' @rdname plotDistToTSS-methods
##' @aliases plotDistToTSS,list-method
##' @exportMethod plotDistToTSS
setMethod("plotDistToTSS", signature(x="list"),
          function(x, distanceColumn="distanceToTSS",
                                     xlab="", ylab="Binding sites (%) (5'->3')",
                                     title="Distribution of transcription factor-binding loci relative to TSS", ...) {
              if (is.null(names(x))) {
                  nn <- paste0("Peak", seq_along(x))
                  warning("input is not a named list, set the name automatically to ", paste(nn, collapse = " "))
                  names(x) <- nn
                  ## stop("input object should be a named list...")
              }

              peakAnno <- lapply(x, as.data.frame)
              ## peakDist <- ldply(peakAnno)
              peakDist <- list_to_dataframe(peakAnno)
              categoryColumn <- ".id"
              plotDistToTSS.data.frame(peakDist, distanceColumn,
                                       xlab, ylab, title, categoryColumn)
          })


##' plotDistToTSS method for \code{csAnno} instance
##'
##' @name plotDistToTSS
##' @docType methods
##' @rdname plotDistToTSS-methods
##' @aliases plotDistToTSS,csAnno,ANY-method
##' @title plotDistToTSS method
##' @param distanceColumn distance column name
##' @param x \code{csAnno} instance
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @param ... additional parameter
##' @return plot
##' @exportMethod plotDistToTSS
##' @usage plotDistToTSS(x,distanceColumn="distanceToTSS", xlab="",
##' ylab="Binding sites (\%) (5'->3')",
##' title="Distribution of transcription factor-binding loci relative to TSS",...)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("plotDistToTSS", signature(x="csAnno"),
          function(x, distanceColumn="distanceToTSS",
                                     xlab="", ylab="Binding sites (%) (5'->3')",
                                     title="Distribution of transcription factor-binding loci relative to TSS", ...) {
              peakDist <- as.data.frame(x)
              categoryColumn <- 1
              plotDistToTSS.data.frame(peakDist, distanceColumn,
                                       xlab, ylab, title, categoryColumn)
          })

