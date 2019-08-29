## @importFrom UpSetR upset
## @importFrom grid viewport
## @importFrom grid pushViewport
## @importFrom grid popViewport
## @importFrom gridBase gridPLT
## @importFrom graphics plot.new
##' @importFrom ggplot2 coord_fixed
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_minimal
##' @author Guangchuang Yu
upsetplot.csAnno <- function(x, order_by = "freq", vennpie=FALSE, vp = list(x=.6, y=.7, width=.8, height=.8)) {
    y <- x@detailGenomicAnnotation
    nn <- names(y)
    y <- as.matrix(y)

    res <- tibble::tibble(anno = lapply(1:nrow(y), function(i) nn[y[i,]]))
    g <- ggplot(res, aes_(x = ~anno)) + geom_bar() +
        xlab(NULL) + ylab(NULL) + theme_minimal() +
        ggupset::scale_x_upset(n_intersections = 20, order_by = order_by) 

    if (!vennpie) return(g)

    f <- function() vennpie(x, cex = .9)

    p <- ggplotify::as.ggplot(f) + coord_fixed() 

    ggplotify::as.ggplot(g) +
        ggimage::geom_subview(subview = p, x = vp$x, y = vp$y, width = vp$width, height = vp$height)


    ## y[y] <- 1
    ## y <- as.data.frame(y)
    ## ## cn <- colnames(y)
    ## ## cn[cn == "fiveUTR"] <- "5 UTR"
    ## ## cn[cn == "threeUTR"] <- "3 UTR"
    ## ## colnames(y) <- cn

    ## if (is.null(sets)) {
    ##     sets <- c("distal_intergenic", "downstream",
    ##               "threeUTR", "fiveUTR", "Intron",
    ##               "Exon", "Promoter")
    ##     if (vennpie && is.null(sets.bar.color)) {
    ##         sets.bar.color <- c("#d95f0e", "#fee0d2", "#98D277",
    ##                             "#6F9E4C", "#fc9272", "#9ecae1", "#ffeda0")
    ##     }
    ## }

    ## if (is.null(sets.bar.color)) {
    ##     sets.bar.color <- "black"
    ## }

    ## if (vennpie) {
    ##     plot.new()
    ##     # grid.rect(gp = gpar(fill="white"))
    ##     upset(y, sets=sets, sets.bar.color=sets.bar.color,
    ##           order.by = order.by, ...)
    ##     pushViewport(vp)
    ##     ##par(plt=gridPLT(), new=TRUE)
    ##     vennpie(x)
    ##     popViewport()
    ## } else {
    ##     upset(y, sets=sets,sets.bar.color=sets.bar.color,
    ##           order.by = order.by, ...)
    ## }
}
