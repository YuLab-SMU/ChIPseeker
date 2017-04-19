---
html_preview: False
output:
  md_document:
    variant: markdown
---

ChIPseeker: ChIP peak Annotation, Comparison, and Visualization
===============================================================

<!-- AddToAny BEGIN -->
<div class="a2a_kit a2a_kit_size_32 a2a_default_style">

<a class="a2a_dd" href="//www.addtoany.com/share"></a>
<a class="a2a_button_facebook"></a> <a class="a2a_button_twitter"></a>
<a class="a2a_button_google_plus"></a>
<a class="a2a_button_pinterest"></a> <a class="a2a_button_reddit"></a>
<a class="a2a_button_sina_weibo"></a> <a class="a2a_button_wechat"></a>
<a class="a2a_button_douban"></a>

</div>

<script async src="//static.addtoany.com/menu/page.js"></script>
<!-- AddToAny END -->
<link rel="stylesheet" href="https://guangchuangyu.github.io/css/font-awesome.min.css">
<link rel="stylesheet" href="https://guangchuangyu.github.io/css/academicons.min.css">

[![](https://img.shields.io/badge/release%20version-1.10.3-blue.svg?style=flat)](https://bioconductor.org/packages/ChIPseeker)
[![](https://img.shields.io/badge/devel%20version-1.11.3-blue.svg?style=flat)](https://github.com/guangchuangyu/ChIPseeker)
[![](https://img.shields.io/badge/download-13514/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ChIPseeker)
[![](https://img.shields.io/badge/download-423/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ChIPseeker)

This package implements functions to retrieve the nearest genes around
the peak, annotate genomic region of the peak, statstical methods for
estimate the significance of overlap among ChIP peak data sets, and
incorporate GEO database for user to compare their own dataset with
those deposited in database. The comparison can be used to infer
cooperative regulation and thus can be used to generate hypotheses.
Several visualization functions are implemented to summarize the
coverage of the peak experiment, average profile and heatmap of peaks
binding to TSS regions, genomic annotation, distance to TSS, and overlap
of peaks or genes.

`ChIPseeker` is released within the
[Bioconductor](https://www.bioconductor.org/packages/ChIPseeker) project
and the source code is hosted on
<a href="https://github.com/GuangchuangYu/ChIPseeker"><i class="fa fa-github fa-lg"></i>
GitHub</a>.

<i class="fa fa-user"></i> Author
---------------------------------

Guangchuang Yu, School of Public Health, The University of Hong Kong.

<a href="https://twitter.com/guangchuangyu"><i class="fa fa-twitter fa-3x"></i></a>
<a href="https://guangchuangyu.github.io/blog_images/biobabble.jpg"><i class="fa fa-wechat fa-3x"></i></a>
<a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=Guangchuang+Yu[Author+-+Full]"><i class="ai ai-pubmed ai-3x"></i></a>
<a href="https://scholar.google.com.hk/citations?user=DO5oG40AAAAJ&hl=en"><i class="ai ai-google-scholar ai-3x"></i></a>
<a href="https://orcid.org/0000-0002-6485-8781"><i class="ai ai-orcid ai-3x"></i></a>
<a href="https://impactstory.org/u/0000-0002-6485-8781"><i class="ai ai-impactstory ai-3x"></i></a>

<i class="fa fa-book"></i> Citation
-----------------------------------

Please cite the following article when using `ChIPseeker`:

[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/btv145-blue.svg?style=flat)](http://dx.doi.org/10.1093/bioinformatics/btv145)
[![](https://img.shields.io/badge/Altmetric-31-blue.svg?style=flat)](https://www.altmetric.com/details/3781087)
[![citation](https://img.shields.io/badge/cited%20by-51-blue.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=12053363057899219488)
[![](https://img.shields.io/badge/ESI-Highly%20Cited%20Paper-blue.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000358173500022)
[![](https://img.shields.io/badge/cited%20in%20Web%20of%20Science%20Core%20Collection-29-green.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000358173500022)

**Yu G**, Wang LG and He QY<sup>\*</sup>. ChIPseeker: an R/Bioconductor
package for ChIP peak annotation, comparison and visualization.
***Bioinformatics***, 2015, 31(14):2382-2383.

<i class="fa fa-pencil"></i> Featured Articles
----------------------------------------------

<img src="https://guangchuangyu.github.io/featured_img/ChIPseeker/heatmap2016.gif" width="650">

<i class="fa fa-hand-o-right"></i> Find out more on
<i class="fa fa-pencil"></i> [Featured
Articles](https://guangchuangyu.github.io/ChIPseeker/featuredArticles/).

<i class="fa fa-download"></i> Installation
-------------------------------------------

Install `ChIPseeker` is easy, follow the guide on the [Bioconductor
page](https://bioconductor.org/packages/ChIPseeker):

``` {.r}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("ChIPseeker")
```

<i class="fa fa-cogs"></i> Overview
-----------------------------------

#### <i class="fa fa-angle-double-right"></i> Annotation

-   retrieve the nearest genes around the peak
-   annotate genomic region of the peak

#### <i class="fa fa-angle-double-right"></i> Comparison

-   estimate the significance of overlap among ChIP peak data sets
-   incorporate GEO database for users to compare their own dataset with
    those deposited in the database

#### <i class="fa fa-angle-double-right"></i> Visualization

-   summarize the coverage of the peak experiment
-   average profile and heatmap of peaks binding to TSS regions
-   genomic annotation
-   distance to TSS
-   overlap of peaks or genes

<!--

## <i class="fa fa-code-fork"></i> Projects that depend on _ChIPseeker_

-->
<i class="fa fa-comment"></i> Feedback
--------------------------------------

<ul class="fa-ul">
    <li><i class="fa-li fa fa-hand-o-right"></i> Please make sure you have followed <a href="https://guangchuangyu.github.io/2016/07/how-to-bug-author/"><strong>the important guide</strong></a> before posting any issue/question</li>
    <li><i class="fa-li fa fa-bug"></i> For bugs or feature requests, please post to <i class="fa fa-github-alt"></i> <a href="https://github.com/GuangchuangYu/ChIPseeker/issues">github issue</a></li>
    <li><i class="fa-li fa fa-question"></i> For user questions, please post to <i class="fa fa-support"></i> <a href="https://support.bioconductor.org">Bioconductor support site</a> or <a href="https://www.biostars.org">Biostars</a></li>
    <li><i class="fa-li fa fa-commenting"></i> Join the group chat in <a href="https://twitter.com/hashtag/ChIPseeker"><i class="fa fa-twitter fa-lg"></i></a> and <a href="http://huati.weibo.com/k/ChIPseeker"><i class="fa fa-weibo fa-lg"></i></a></li>

</ul>
