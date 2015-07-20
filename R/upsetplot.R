##' @importFrom UpSetR upset
##' @author Guangchuang Yu
upsetplot.csAnno <- function(x, sets=NULL, order.matrix = "freq", ...) {
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
    }
    upset(y, sets=sets, order.matrix = order.matrix, ...)
}
