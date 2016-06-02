# ChIP peak Annotation, Comparison, and Visualization #

<!--[![Build Status](https://travis-ci.org/GuangchuangYu/ChIPseeker.svg?branch=master)](https://travis-ci.org/GuangchuangYu/ChIPseeker)-->
[![platform](http://www.bioconductor.org/shields/availability/devel/ChIPseeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/ChIPseeker.html#archives)
[![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/ChIPseeker.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/ChIPseeker/)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/ChIPseeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/ChIPseeker.html#since)
[![post](http://www.bioconductor.org/shields/posts/ChIPseeker.svg)](https://support.bioconductor.org/t/ChIPseeker/)
[![commit](http://www.bioconductor.org/shields/commits/bioc/ChIPseeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/ChIPseeker.html#svn_source)
[![download](http://www.bioconductor.org/shields/downloads/ChIPseeker.svg)](https://bioconductor.org/packages/stats/bioc/ChIPseeker)


This package implements functions to retrieve the nearest genes around the peak, annotate genomic region of the peak, statstical methods for estimate the significance of overlap among ChIP peak data sets, and incorporate GEO database for user to compare their own dataset with those deposited in database. The comparison can be used to infer cooperative regulation and thus can be used to generate hypotheses. Several visualization functions are implemented to summarize the coverage of the peak experiment, average profile and heatmap of peaks binding to TSS regions, genomic annotation, distance to TSS, and overlap of peaks or genes.

## Authors ##

Guangchuang YU, School of Public Health, The University of Hong Kong <http://guangchuangyu.github.io>

## Citation ##

Please cite the following article when using `ChIPseeker`:

```
Yu G, Wang LG and He QY.
ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization.
Bioinformatics, 2015, 31(14):2382-2383.
doi: 10.1093/bioinformatics/btv145
```

URL: [http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btv145](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btv145)

## License ##

All source code is copyright, under the Artistic-2.0 License.
For more information on Artistic-2.0 License see [http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)

## Installation ##

To install:
 * the latest released version:
   `biocLite("ChIPseeker")`
 * the latest development version:
   `devtools::install_github("GuangchuangYu/ChIPseeker")`

## Documentation ##

Find out more at [http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html](http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html) and check out the [vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html).

To view the vignette of `ChIPseeker` installed in your system, start R and enter:
```r
vignette("ChIPseeker", "ChIPseeker")
```

More documents can be found in <http://guangchuangyu.github.io/tags/chipseeker>.

## Bugs/Feature requests ##

 - If you have any, [let me know](https://github.com/GuangchuangYu/ChIPseeker/issues). Thx!


