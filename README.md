#  ChIPseeker for ChIP peak Annotation, Comparison, and Visualization #

This package implements functions to retrieve the nearest genes around the peak, annotate genomic region of the peak. Statstical methods for estimate the significance of overlap among ChIP peak data sets, and incorporate GEO database for user to compare the own dataset with those deposited in database. The comparison can be used to infer cooperative regulation and thus can be used to generate hypotheses. Several visualization functions are implemented to summarize the coverage of the peak experiment, average profile and heatmap of peaks binding to TSS regions, genomic annotation, distance to TSS, and overlap of peaks or genes.

## Authors ##

Guangchuang YU, School of Public Health, The University of Hong Kong [http://ygc.name](http://ygc.name)

## License ## 

All source code is copyright, under the Artistic-2.0 License.
For more information on Artistic-2.0 License see [http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)

## Installation ##

To install:
 * the latest released version:
   `require(BiocInstaller)
   biocLite("ChIPseeker")
   `
 * the latest development version:
   `require(devtools)
   install_github("GuangchuangYu/ChIPseeker") 
   `

Find out more at [http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html](http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html) and check out the vignettes.

