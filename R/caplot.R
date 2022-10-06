as.caplot <- function(object) {
  if (inherits(object, "caplot"))
    return(object)
  
  if (!inherits(object, 'aplot')) {
    stop("input should be a 'aplot' object.")
  }
  
  attr(object,"class") <- "caplot"
  
  return(object)
}


##' @method print caplot
##' @importFrom patchwork plot_layout
##' @importFrom patchwork plot_spacer
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 theme_void
##' @export
print.caplot <- function(x, ...) {
  grid.draw(x)
}


##' @importFrom grid grid.draw
##' @importFrom aplot as.patchwork
##' @method grid.draw caplot
##' @export
grid.draw.caplot <- function(x, recoding = TRUE) {
  attr(x,"class") <- "aplot"
  grid::grid.draw(as.patchwork(x,align="y"))
}
