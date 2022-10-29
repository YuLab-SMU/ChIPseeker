# ChIPseeker: ChIP peak Annotation, Comparison, and Visualization

<img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/ChIPseeker/ChIPseeker.png" height="200" align="right" />

[![](https://img.shields.io/badge/release%20version-1.32.1-green.svg)](https://www.bioconductor.org/packages/ChIPseeker)
[![](https://img.shields.io/badge/devel%20version-1.33.4-green.svg)](https://github.com/guangchuangyu/ChIPseeker)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/ChIPseeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/ChIPseeker.html#since)
[![Say
Thanks!](https://img.shields.io/badge/Say%20Thanks-!-1EAEDB.svg)](https://saythanks.io/to/GuangchuangYu)

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/GuangchuangYu/ChIPseeker/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/ChIPseeker/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2022--10--29-green.svg)](https://github.com/GuangchuangYu/ChIPseeker/commits/master)
[![GitHub
forks](https://img.shields.io/github/forks/GuangchuangYu/ChIPseeker.svg)](https://github.com/GuangchuangYu/ChIPseeker/network)
[![GitHub
stars](https://img.shields.io/github/stars/GuangchuangYu/ChIPseeker.svg)](https://github.com/GuangchuangYu/ChIPseeker/stargazers)

[![platform](http://www.bioconductor.org/shields/availability/devel/ChIPseeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/ChIPseeker.html#archives)
[![Build
Status](http://www.bioconductor.org/shields/build/devel/bioc/ChIPseeker.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/ChIPseeker/)
[![Linux/Mac Travis Build
Status](https://img.shields.io/travis/GuangchuangYu/ChIPseeker/master.svg?label=Mac%20OSX%20%26%20Linux)](https://travis-ci.org/GuangchuangYu/ChIPseeker)
[![AppVeyor Build
Status](https://img.shields.io/appveyor/ci/Guangchuangyu/ChIPseeker/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/ChIPseeker)

This package implements functions to retrieve the nearest genes around
the peak, annotate genomic region of the peak, statstical methods for
estimate the significance of overlap among ChIP peak data sets, and
incorporate GEO database for user to compare the own dataset with those
deposited in database. The comparison can be used to infer cooperative
regulation and thus can be used to generate hypotheses. Several
visualization functions are implemented to summarize the coverage of the
peak experiment, average profile and heatmap of peaks binding to TSS
regions, genomic annotation, distance to TSS, and overlap of peaks or
genes.

## :writing_hand: Authors

Guangchuang YU

School of Basic Medical Sciences, Southern Medical University

<https://yulab-smu.top>

If you use [ChIPseeker](http://bioconductor.org/packages/ChIPseeker) in
published research, please cite:

-   Q Wang<sup>\#</sup>, M Li<sup>\#</sup>, T Wu, L Zhan, L Li, M Chen,
    W Xie, Z Xie, E Hu, S Xu, **G Yu**<sup>\*</sup>. [Exploring
    epigenomic datasets by
    ChIPseeker](https://onlinelibrary.wiley.com/share/author/GYJGUBYCTRMYJFN2JFZZ?target=10.1002/cpz1.585).
    ***Current Protocols***, 2022, 2(10): e585.
-   **G Yu**<sup>\*</sup>, LG Wang, QY He<sup>\*</sup>. [ChIPseeker: an
    R/Bioconductor package for ChIP peak annotation, comparision and
    visualization](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btv145).
    ***Bioinformatics***. 2015, 31(14):2382-2383.

## :arrow_double_down: Installation

Get the released version from Bioconductor:

``` r
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("ChIPseeker")
```

Or the development version from github:

``` r
## install.packages("devtools")
devtools::install_github("YuLab-SMU/ChIPseeker")
```

## Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
