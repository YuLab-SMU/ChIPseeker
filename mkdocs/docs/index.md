<!-- AddToAny BEGIN -->
<div class="a2a_kit a2a_kit_size_32 a2a_default_style">
<a class="a2a_dd" href="//www.addtoany.com/share"></a>
<a class="a2a_button_facebook"></a>
<a class="a2a_button_twitter"></a>
<a class="a2a_button_google_plus"></a>
<a class="a2a_button_pinterest"></a>
<a class="a2a_button_reddit"></a>
<a class="a2a_button_sina_weibo"></a>
<a class="a2a_button_wechat"></a>
<a class="a2a_button_douban"></a>
</div>
<script async src="//static.addtoany.com/menu/page.js"></script>
<!-- AddToAny END -->

<link rel="stylesheet" href="https://guangchuangyu.github.io/css/font-awesome.min.css">

[![releaseVersion](https://img.shields.io/badge/release%20version-1.8.9-blue.svg?style=flat)](https://bioconductor.org/packages/ChIPseeker)
[![develVersion](https://img.shields.io/badge/devel%20version-1.11.0-blue.svg?style=flat)](https://github.com/GuangchuangYu/ChIPseeker)
[![total](https://img.shields.io/badge/downloads-17954/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ChIPseeker)
[![month](https://img.shields.io/badge/downloads-555/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ChIPseeker)

This package implements functions to retrieve the nearest genes around the peak, annotate genomic region of the peak, statstical methods for estimate the significance of overlap among ChIP peak data sets, and incorporate GEO database for user to compare their own dataset with those deposited in database. The comparison can be used to infer cooperative regulation and thus can be used to generate hypotheses. Several visualization functions are implemented to summarize the coverage of the peak experiment, average profile and heatmap of peaks binding to TSS regions, genomic annotation, distance to TSS, and overlap of peaks or genes.

`ChIPseeker` is released within the [Bioconductor](https://www.bioconductor.org/packages/ChIPseeker) project and the source code is hosted on <a href="https://github.com/GuangchuangYu/ChIPseeker"><i class="fa fa-github fa-lg"></i> GitHub</a>.

## <i class="fa fa-user"></i> Author

Guangchuang Yu, School of Public Health, The University of Hong Kong.

## <i class="fa fa-book"></i> Citation

Please cite the following article when using `ChIPseeker`:

[![doi](https://img.shields.io/badge/doi-10.1093/bioinformatics/btv145-blue.svg?style=flat)](http://dx.doi.org/10.1093/bioinformatics/btv145)
[![citation](https://img.shields.io/badge/cited%20by-20-blue.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=12053363057899219488)
[![Altmetric](https://img.shields.io/badge/Altmetric-23-blue.svg?style=flat)](https://www.altmetric.com/details/3781087)

__Yu G__, Wang LG and He QY<sup>*</sup>. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. **_Bioinformatics_**, 2015, 31(14):2382-2383.


## <i class="fa fa-pencil"></i> Featured Articles

<img src="https://guangchuangyu.github.io/featured_img/ChIPseeker/heatmap2016.gif" width="650">


<i class="fa fa-hand-o-right"></i> Find out more on <i class="fa fa-pencil"></i> [Featured Articles](https://guangchuangyu.github.io/ChIPseeker/featuredArticles/).

## <i class="fa fa-download"></i> Installation

Install `ChIPseeker` is easy, follow the guide on the [Bioconductor page](https://bioconductor.org/packages/ChIPseeker):


```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("ChIPseeker")
```

## <i class="fa fa-cogs"></i> Overview

#### <i class="fa fa-angle-double-right"></i> Annotation

+ retrieve the nearest genes around the peak
+ annotate genomic region of the peak

#### <i class="fa fa-angle-double-right"></i> Comparison

+ estimate the significance of overlap among ChIP peak data sets
+ incorporate GEO database for users to compare their own dataset with those deposited in the database

#### <i class="fa fa-angle-double-right"></i> Visualization

+ summarize the coverage of the peak experiment
+ average profile and heatmap of peaks binding to TSS regions
+ genomic annotation
+ distance to TSS
+ overlap of peaks or genes

<!--

## <i class="fa fa-code-fork"></i> Projects that depend on _ChIPseeker_

-->


## <i class="fa fa-comment"></i> Feedback

<ul class="fa-ul">
	<li><i class="fa-li fa fa-hand-o-right"></i> Please make sure you <a href="https://guangchuangyu.github.io/2016/07/how-to-bug-author">follow the guide</a> before posting any issue/question</li>
	<li><i class="fa-li fa fa-bug"></i> For bugs or feature requests, please post to <i class="fa fa-github-alt"></i> <a href="https://github.com/GuangchuangYu/ChIPseeker/issues">github issue</a></li>
	<li><i class="fa-li fa fa-question"></i> For user questions, please post to <i class="fa fa-support"></i> <a href="https://support.bioconductor.org">Bioconductor support site</a> or <a href="https://www.biostars.org">Biostars</a></li>
	<li><i class="fa-li fa fa-commenting"></i> Join the group chat in <a href="https://twitter.com/hashtag/ChIPseeker"><i class="fa fa-twitter fa-lg"></i></a> and <a href="http://huati.weibo.com/k/ChIPseeker"><i class="fa fa-weibo fa-lg"></i></a></li>
</ul>


