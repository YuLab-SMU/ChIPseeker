##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "  ",
                "For help: https://guangchuangyu.github.io/", pkgname, "\n\n")

  citation <- paste0("If you use ", pkgname, " in published research, please cite:\n",
                     "Guangchuang Yu, Li-Gen Wang, Qing-Yu He. ",
                     "ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. ",
                     "Bioinformatics 2015, 31(14):2382-2383")

  packageStartupMessage(paste0(msg, citation))

}
