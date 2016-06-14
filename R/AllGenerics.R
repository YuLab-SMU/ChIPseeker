##' vennpie method generics
##'
##'
##' @docType methods
##' @name vennpie
##' @rdname vennpie-methods
##' @export
setGeneric("vennpie", function(x, r=0.2, ...) standardGeneric("vennpie"))


##' plotDistToTSS method generics
##'
##'
##' @docType methods
##' @name plotDistToTSS
##' @rdname plotDistToTSS-methods
##' @export
setGeneric("plotDistToTSS", function(x, distanceColumn="distanceToTSS",
                                     xlab="", ylab="Binding sites (%) (5'->3')",
                                     title="Distribution of transcription factor-binding loci relative to TSS", ...)
           standardGeneric("plotDistToTSS"))

##' plotAnnoBar method generics
##'
##'
##' @docType methods
##' @name plotAnnoBar
##' @rdname plotAnnoBar-methods
##' @export
setGeneric("plotAnnoBar", function(x,
                                   xlab="",
                                   ylab="Percentage(%)",
                                   title="Feature Distribution",
                                   ...)
           standardGeneric("plotAnnoBar"))


##' plotAnnoPie method generics
##'
##'
##' @docType methods
##' @name plotAnnoPie
##' @rdname plotAnnoPie-methods
##' @export
setGeneric("plotAnnoPie", 
           function(x, 
                    ndigit=2,
                    cex=0.9,
                    col=NA,
                    legend.position="rightside",
                    pie3D=FALSE,
                    ...)
           standardGeneric("plotAnnoPie"))
