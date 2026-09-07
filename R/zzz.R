##' Package initialization function
##'
##' This function is automatically called when the ChIPseeker package is attached.
##' It displays a startup message and sets default package options.
##'
##' @description
##' The \code{.onAttach} function is a standard R package hook that runs when
##' \code{library(ChIPseeker)} or \code{require(ChIPseeker)} is called. It
##' initializes package-specific options that control various behaviors throughout
##' the package.
##'
##' @details
##' The function performs the following initialization:
##' \enumerate{
##'   \item Displays a package startup message using \code{yulab_msg()}
##'   \item Sets default package options:
##'     \itemize{
##'       \item \code{ChIPseeker.downstreamDistance}: Default distance (in base pairs)
##'             for downstream annotation. Default is 300bp
##'       \item \code{ChIPseeker.ignore_1st_exon}: Whether to ignore first exon
##'             in annotation. Default is FALSE
##'       \item \code{ChIPseeker.ignore_1st_intron}: Whether to ignore first intron
##'             in annotation. Default is FALSE
##'       \item \code{ChIPseeker.ignore_downstream}: Whether to ignore downstream
##'             regions in annotation. Default is FALSE
##'       \item \code{ChIPseeker.ignore_promoter_subcategory}: Whether to ignore
##'             promoter subcategories in annotation. Default is FALSE
##'     }
##'   \item Sets the \code{aplot_align} option to 'y' for alignment in plotting
##'         functions
##' }
##'
##' These options can be modified by users using \code{options()} to customize
##' package behavior. For example:
##' \code{options(ChIPseeker.downstreamDistance = 5000)}
##'
##' @param libname character, the library directory where the package was found
##' @param pkgname character, the name of the package
##' @return No return value. Called for its side effects (setting options and
##'   displaying message)
##' @seealso \code{\link[base]{options}} for modifying package options,
##'   \code{\link[base]{getOption}} for retrieving option values
##' @importFrom yulab.utils yulab_msg
##' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(yulab_msg(pkgname))

  options(ChIPseeker.downstreamDistance = 300)
  options(ChIPseeker.ignore_1st_exon = FALSE)
  options(ChIPseeker.ignore_1st_intron = FALSE)
  options(ChIPseeker.ignore_downstream = FALSE)
  options(ChIPseeker.ignore_promoter_subcategory= FALSE)

  options(aplot_align = 'y')

}

