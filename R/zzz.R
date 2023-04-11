##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "\n\n") ## "  ",
              ##  "For help: https://guangchuangyu.github.io/software/", pkgname, "\n\n")

  citation <- paste0("If you use ", pkgname, " in published research, please cite:\n",
                     "Qianwen Wang, Ming Li, Tianzhi Wu, Li Zhan, Lin Li, Meijun Chen, Wenqin Xie, Zijing Xie, Erqiang Hu, Shuangbin Xu, Guangchuang Yu. ",
                     "Exploring epigenomic datasets by ChIPseeker. ",
                     "Current Protocols 2022, 2(10): e585")


  packageStartupMessage(paste0(msg, citation))

  options(ChIPseeker.downstreamDistance = 300)
  options(ChIPseeker.ignore_1st_exon = FALSE)
  options(ChIPseeker.ignore_1st_intron = FALSE)
  options(ChIPseeker.ignore_downstream = FALSE)
  options(ChIPseeker.ignore_promoter_subcategory= FALSE)
  
  options(aplot_align = 'y')

}

