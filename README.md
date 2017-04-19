ChIPseeker: ChIP peak Annotation, Comparison, and Visualization
===============================================================

[![](https://img.shields.io/badge/release%20version-1.10.3-green.svg?style=flat)](https://bioconductor.org/packages/ChIPseeker) [![](https://img.shields.io/badge/devel%20version-1.11.3-green.svg?style=flat)](https://github.com/guangchuangyu/ChIPseeker) [![Bioc](http://www.bioconductor.org/shields/years-in-bioc/ChIPseeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/ChIPseeker.html#since) [![](https://img.shields.io/badge/download-13514/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ChIPseeker) [![](https://img.shields.io/badge/download-423/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ChIPseeker)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![codecov](https://codecov.io/gh/GuangchuangYu/ChIPseeker/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/ChIPseeker/) [![Last-changedate](https://img.shields.io/badge/last%20change-2017--04--19-green.svg)](https://github.com/GuangchuangYu/ChIPseeker/commits/master) [![commit](http://www.bioconductor.org/shields/commits/bioc/ChIPseeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/ChIPseeker.html#svn_source) [![GitHub forks](https://img.shields.io/github/forks/GuangchuangYu/ChIPseeker.svg)](https://github.com/GuangchuangYu/ChIPseeker/network) [![GitHub stars](https://img.shields.io/github/stars/GuangchuangYu/ChIPseeker.svg)](https://github.com/GuangchuangYu/ChIPseeker/stargazers)

[![platform](http://www.bioconductor.org/shields/availability/devel/ChIPseeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/ChIPseeker.html#archives) [![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/ChIPseeker.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/ChIPseeker/) [![Linux/Mac Travis Build Status](https://img.shields.io/travis/GuangchuangYu/ChIPseeker/master.svg?label=Mac%20OSX%20%26%20Linux)](https://travis-ci.org/GuangchuangYu/ChIPseeker) [![AppVeyor Build Status](https://img.shields.io/appveyor/ci/Guangchuangyu/ChIPseeker/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/ChIPseeker) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-green.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-chipseeker/README.html)

This package implements functions to retrieve the nearest genes around the peak, annotate genomic region of the peak, statstical methods for estimate the significance of overlap among ChIP peak data sets, and incorporate GEO database for user to compare their own dataset with those deposited in database. The comparison can be used to infer cooperative regulation and thus can be used to generate hypotheses. Several visualization functions are implemented to summarize the coverage of the peak experiment, average profile and heatmap of peaks binding to TSS regions, genomic annotation, distance to TSS, and overlap of peaks or genes.

[![Twitter](https://img.shields.io/twitter/url/https/github.com/GuangchuangYu/ChIPseeker.svg?style=social)](https://twitter.com/intent/tweet?hashtags=ChIPseeker&url=http://bioinformatics.oxfordjournals.org/content/31/14/2382&screen_name=guangchuangyu)

------------------------------------------------------------------------

Please cite the following article when using `ChIPseeker`:

***Yu G***, Wang LG and He QY<sup>\*</sup>. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. ***Bioinformatics*** 2015, 31(14):2382-2383.

[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/btv145-green.svg?style=flat)](http://dx.doi.org/10.1093/bioinformatics/btv145) [![](https://img.shields.io/badge/Altmetric-31-green.svg?style=flat)](https://www.altmetric.com/details/3781087) [![citation](https://img.shields.io/badge/cited%20by-51-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=12053363057899219488) [![](https://img.shields.io/badge/ESI-Highly%20Cited%20Paper-green.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000358173500022)

------------------------------------------------------------------------

For details, please visit our project website, <https://guangchuangyu.github.io/ChIPseeker>.

-   [Documentation](https://guangchuangyu.github.io/ChIPseeker/documentation/)
-   [Featured Articles](https://guangchuangyu.github.io/ChIPseeker/featuredArticles/)
-   [Feedback](https://guangchuangyu.github.io/ChIPseeker/#feedback)

### Citation

[![citation](https://img.shields.io/badge/cited%20by-51-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=12053363057899219488) [![](https://img.shields.io/badge/ESI-Highly%20Cited%20Paper-green.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000358173500022)

       +-+------------+-----------+------------+-----------+---+
    30 +                          *                            +
       |                                                       |
       |                                                       |
    25 +                                                       +
       |                                                       |
    20 +                                                       +
       |                                                       |
    15 +                                                   *   +
       |                                                       |
    10 +                                                       +
       |                                                       |
     5 + *                                                     +
       +-+------------+-----------+------------+-----------+---+
       2015        2015.5       2016        2016.5       2017   

### Download stats

[![download](http://www.bioconductor.org/shields/downloads/ChIPseeker.svg)](https://bioconductor.org/packages/stats/bioc/ChIPseeker) [![](https://img.shields.io/badge/download-13514/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ChIPseeker) [![](https://img.shields.io/badge/download-423/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ChIPseeker)

         +--------+--------------+-------------+--------------+-------------+--------------+-----------+
         |                                                                              *              |
    1000 +                                                                  *                          +
         |                                                                                             |
         |                                                                       *                     |
         |                                                                                             |
         |                                                           *    *   *                        |
     800 +                                                                                             +
         |                                                                                   *         |
         |                                               *             *                   *      *    |
         |                                                 *    * *                                    |
         |                                       *  * *       *                       *        *       |
     600 +                                                                                             +
         |                           *  * *  *                                     *                   |
         |        *                                                                                    |
         |               *         *                                                                   |
         |                                                                                             |
     400 +      *      *                       *                                                       +
         |                  * *                                                                        |
         |           *                                                                                 |
         |                       *                                                                     |
         |                                                                                             |
     200 +   *                                                                                         +
         +--------+--------------+-------------+--------------+-------------+--------------+-----------+
               2014.5          2015         2015.5          2016         2016.5          2017
