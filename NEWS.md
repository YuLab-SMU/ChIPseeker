# ChIPseeker 1.15.4

+ plotAnnoBar now visualize barplot according to the order of input list
  (y-axis) (2018-02-27, Tue)
    - <https://github.com/GuangchuangYu/ChIPseeker/issues/73>
+ follow renaming of RangesList class -> IntegerRangesList in IRanges v2.13.12
    - <https://github.com/GuangchuangYu/ChIPseeker/commit/b62d7922fb61e58620bbb685e4def4fb863c8e81>

# ChIPseeker 1.15.3

+ options to ignore '1st exon', '1st intron', 'downstream' and promoter
  subcategory when summarizing result and visualization (2018-01-09, Tue)
    - <https://support.bioconductor.org/p/104676/#104689>
+ throw msg of 'file not found and skip' when requested url is not available
  when downloading BED file from GEO (2017-12-28, Thu)
    - <https://support.bioconductor.org/p/104491/#104507>
+ bug fixed of getGene (2017-12-27, Wed)
