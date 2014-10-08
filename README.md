```
  ____ _     ___ ____                _             
 / ___| |__ |_ _|  _ \ ___  ___  ___| | _____ _ __ 
| |   | '_ \ | || |_) / __|/ _ \/ _ \ |/ / _ \ '__|
| |___| | | || ||  __/\__ \  __/  __/   <  __/ |   
 \____|_| |_|___|_|   |___/\___|\___|_|\_\___|_|   
                                                                                                      
```
# ChIP peak Annotation, Comparison, and Visualization #

This package implements functions to retrieve the nearest genes around the peak, annotate genomic region of the peak, statstical methods for estimate the significance of overlap among ChIP peak data sets, and incorporate GEO database for user to compare the own dataset with those deposited in database. The comparison can be used to infer cooperative regulation and thus can be used to generate hypotheses. Several visualization functions are implemented to summarize the coverage of the peak experiment, average profile and heatmap of peaks binding to TSS regions, genomic annotation, distance to TSS, and overlap of peaks or genes.

## Authors ##

Guangchuang YU, School of Public Health, The University of Hong Kong [http://ygc.name](http://ygc.name)

## License ##

All source code is copyright, under the Artistic-2.0 License.
For more information on Artistic-2.0 License see [http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)

## Installation ##

To install:
 * the latest released version:
   `biocLite("ChIPseeker")`
 * the latest development version:
   `install_github("GuangchuangYu/ChIPseeker")`

## Documentation ##

+ [Peak Annotation](http://ygc.name/2014/04/13/chipseeker-for-chip-peak-annotation/)
+ [Multiple Annotation](http://ygc.name/2014/10/01/multiple-annotation-in-chipseeker/)
+ [Visualization Methods](http://ygc.name/2014/04/30/visualization-methods-in-chipseeker/)

Find out more at [http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html](http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html) and check out the [vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.pdf).

To view the vignette of `ChIPseeker` installed in your system, start R and enter:
```r
browseVignettes("ChIPseeker")
```

## Bugs/Feature requests ##

 - If you have any, [let me know](https://github.com/GuangchuangYu/ChIPseeker/issues). Thx!


