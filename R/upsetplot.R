##' @importFrom UpSetR upset
##' @importFrom grid viewport
##' @importFrom grid pushViewport
##' @importFrom grid popViewport
##' @importFrom gridBase gridPLT
##' @importFrom graphics plot.new
##' @author Guangchuang Yu
upsetplot.csAnno <- function(x, sets=NULL, order.matrix = "freq", sets.bar.color=NULL, vennpie=FALSE, ...) {
    y <- x@detailGenomicAnnotation
    y <- as.matrix(y)
    y[y] <- 1
    y <- as.data.frame(y)
    ## cn <- colnames(y)
    ## cn[cn == "fiveUTR"] <- "5 UTR"
    ## cn[cn == "threeUTR"] <- "3 UTR"
    ## colnames(y) <- cn

    if (is.null(sets)) {
        sets <- c("distal_intergenic", "downstream",
                  "threeUTR", "fiveUTR", "Intron",
                  "Exon", "Promoter")
        if (vennpie && is.null(sets.bar.color)) {
            sets.bar.color <- c("#d95f0e", "#fee0d2", "#98D277",
                                "#6F9E4C", "#fc9272", "#9ecae1", "#ffeda0")
        } 
    }
    
    if (is.null(sets.bar.color)) {
        sets.bar.color <- "black"
    }

    if (vennpie) {
        plot.new()
        # grid.rect(gp = gpar(fill="white"))
        upset(y, sets=sets, sets.bar.color=sets.bar.color,
              order.matrix = order.matrix, ...)
        pushViewport(viewport(x=.6, y=.7, width=.6, height=.6))
        par(plt=gridPLT(), new=TRUE)
        vennpie(x)
        popViewport()
    } else {
        upset(y, sets=sets,sets.bar.color=sets.bar.color,
              order.matrix = order.matrix, ...)
    }
}
