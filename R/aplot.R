##' @method print aplot
##' @importFrom patchwork plot_layout
##' @importFrom patchwork plot_spacer
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 theme_void
##' @export
print.aplot <- function(x, ...) {
  grid.draw(x)
}


##' @importFrom grid grid.draw
##' @importFrom aplot as.patchwork
##' @method grid.draw aplot
##' @export
grid.draw.aplot <- function(x, recoding = TRUE) {
  attr(x,"class") <- "aplot"
  grid::grid.draw(as.patchwork(x,align="y"))
}

##' @importFrom ggplot2 ggplotGrob
##' @importFrom patchwork patchworkGrob
##' @importFrom aplot as.patchwork
aplotGrob <- function(x) {
  res <- as.patchwork(x,align="y")
  patchworkGrob(res)
}