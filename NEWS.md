# ChIPseeker 1.49.3

+ `annotatePeak()` keeps `annotation`, `geneChr/geneStart/geneEnd`, `geneId` and
  `transcriptId` consistent with each other at transcript level (issue #252). The
  transcript id of an exon/intron/UTR hit was taken from
  `names(genomicRegion)[subjectIndex]`, indexing the unlisted ranges with the
  names of the `GRangesList`; this mostly returned NA and occasionally a wrong
  transcript, so a peak could be reported with the metadata of an unrelated
  transcript (even on another chromosome). The ids are now expanded before
  indexing, and the `Promoter` branch writes the feature id reset into `anno`
  instead of a stale local variable, so promoter peaks keep the nearest-TSS
  transcript. The same alignment is now applied at `level = "gene"`: the gene of
  the overlapping exon/intron/UTR is resolved with `TXID2EG()` and used instead
  of the nearest gene, so `geneId`/`geneStart`/`geneEnd`/`distanceToTSS` follow
  the annotation there as well. (2026-10-01, Thu)
+ `seq2gene()` no longer fails with `$ operator is invalid for atomic vectors`
  when none of the queried regions overlaps an exon
  (`getGenomicAnnotation.internal()` returns `NA` then); host-gene extraction is
  skipped and the nearest/flanking genes are still reported (issue #248).
  (2026-10-01, Thu)
+ `getTagMatrix()`/`plotPeakProf2()` no longer fail with
  `Error in cursor:(cursor + seq - 1) : result would be too long a vector` when
  `type = "body"` is used with a flank extension shorter than 1kb (issue #250).
  The share of bins of a flank is derived from its actual length (500bp now
  contributes 5% of the bins instead of being rounded down to zero), a
  non-empty flank always gets at least one column, and the binning loops were
  replaced by edge-based averaging (which also fixes an off-by-one divisor in
  the last bin of the upstream/body sections). The x-axis breaks of
  `plotPeakProf2()` are derived from the same column layout. (2026-10-01, Thu)
+ `downloadGEObedFiles()`/`downloadGSMbedFiles()` rewrite the `ftp://` urls of
  `gsminfo$supplementary_file` to `https://` before downloading and report the
  underlying download error when a file cannot be fetched (issue #254).
  (2026-10-01, Thu)
+ Documented the distance definition of `flank_gene_distances` reported by
  `annotatePeak(..., addFlankGeneInfo = TRUE)`: a distance of 0 means that the
  peak overlaps the feature range - at `level = "transcript"` the feature is the
  whole transcript, which is why most entries can be 0 - while non-overlapping
  peaks get the signed distance to the feature TSS (issue #235).
  (2026-10-01, Thu)

# ChIPseeker 1.49.2

+ `annotatePeak()` now reports `geneChr` and `geneStrand` as characters instead of
  factor codes. `as.data.frame()` returns 'seqnames'/'strand' as factors and
  assigning a factor into `mcols()` dropped the class, so both columns came out
  as integers (e.g. `geneStrand` = 1/2 instead of +/-, and a wrong `geneChr`
  whenever the seqlevels were not in numeric order) (issues #233, #247).
  (2026-09-15, Tue)
+ `annotatePeak()` now warns when peaks are dropped because no feature of `TxDb`
  can be found for them (e.g. peaks on contigs/scaffolds without genes, or a
  seqlevels style mismatch). Previously they disappeared silently (issues
  #251, #258). (2026-09-15, Tue)
+ `annotatePeak(..., sameStrand = TRUE)` no longer assigns a peak to a feature on
  the opposite strand. Overlap detection in `getNearestFeatureIndicesAndDistances()`
  was calling `findOverlaps()` with `unstrand(features)`, which bypassed `sameStrand`
  and overrode the strand-aware nearest-feature result (issues #257, #258). Peaks
  with ambiguous strand (`*`) are unaffected and still match features on any strand.
  (2026-09-14, Mon)
+ `plotAnnoBar()` no longer uses the deprecated `ggplot2::aes_string()` (issue #268).
  It follows the tidy evaluation idiom already used by `plotDistToTSS()`.
  (2026-09-14, Mon)
+ `upsetplot()` no longer uses the deprecated `ggplot2::aes_()`; the x aesthetic is
  now mapped with `.data$anno`, the same tidy evaluation idiom used by
  `plotDistToTSS()`. (2026-09-17, Thu)
+ dropped unused `@importFrom` directives that have no call site in the package
  (`ggplot2::geom_text`, `ggplot2::geom_segment`, `ggplot2::scale_fill_brewer`,
  `ggplot2::scale_fill_hue`), so the generated NAMESPACE imports less.
  (2026-09-17, Thu)

# ChIPseeker 1.49.1

+ Fixed bug in `getNearestFeatureIndicesAndDistances()` where results for
  `overlap == "all"` were silently overridden by the `overlap == "TSS"` branch,
  causing the two modes to behave identically. (2026-09-07, Mon)
+ Comprehensive documentation updates for core functions (`annotatePeak()`,
  `getGenomicAnnotation()`, `getTagMatrix()` and related, `readPeakFile()`,
  `seq2gene()`, dplyr verb extensions, plotting functions, etc.).
+ Simplified redundant logic in `seq2gene()` by merging the promoter and
  flanking-gene extraction into a single condition (no behavior change).

# ChIPseeker 1.48.0

+ Bioconductor RELEASE_3_23 (2026-04-29, Wed)

# ChIPseeker 1.47.1

+ fixed issue in 'test-txdb.R' as 'TxDb.Hsapiens.UCSC.hg19.knownGene' changes its transcript ID from UCSC (e.g., uc002qsd.4) to Ensembl (e.g., ENST00000487630.1_3) (2025-11-04, Tue)

# ChIPseeker 1.46.0

+ Bioconductor RELEASE_3_22 (2025-11-01, Sat)

# ChIPseeker 1.45.2

+ new cache mechanism from 'yulab.utils' (2025-10-15, Wed)

# ChIPseeker 1.44.0

+ Bioconductor RELEASE_3_21 (2025-04-17, Thu)

# ChIPseeker 1.42.0

+ Bioconductor RELEASE_3_20 (2024-10-30, Wed)

# ChIPseeker 1.41.3

+ Better `covplot()`. Support universal chromosome names, and keep the default order of multiple peaks when plot a list of `GRanges` object.
+ Robust `generate_colors()`. Edit the logical of decision, and can validate color code automatically.
+ Extend dplyr verbs (`filter()`, `mutate()`, `arrange()`, `rename()`) to peak (`GRanges` object or `data.frame`), see #242.

# ChIPseeker 1.41.2

+ Enhancement of `plotDistToTSS()`, see #241.

# ChIPseeker 1.41.1

+ use `yulab.utils::yulab_msg()` for startup message (2024-07-26, Fri)

# ChIPseeker 1.40.0

+ Bioconductor RELEASE_3_19 (2024-05-15, Wed)

# ChIPseeker 1.38.0

+ Bioconductor RELEASE_3_18 (2023-10-25, Wed)

# ChIPseeker 1.36.0

+ Bioconductor RELEASE_3_17 (2023-05-03, Wed)

# ChIPseeker 1.35.3

+ fixed R check by removing calling `BiocStyle::Biocpkg()` in vignette, instead we use `yulab.utils::Biocpkg()` (2023-04-11, Tue)

# ChIPseeker 1.35.2

+ fixed R check by adding 'prettydoc' to Suggests (2023-04-04, Tue)

# ChIPseeker 1.35.1

+ use `ggplot` to plot heatmap (2022-12-30, Fri, #203)
+ update startup message to display the 'Current Protocols (2022)' paper. 

# ChIPseeker 1.34.0

+ Bioconductor RELEASE_3_16 (2022-11-02, Wed)


# ChIPseeker 1.33.4

+ add citation Q. Wang (2022) (2022-10-29, Sat)

# ChIPseeker 1.33.3

+ allows passing user defined color to `vennpie()` (2022-10-20, Thu, #202, #207)
+ add `columns` paramter to `annotatePeak()` to better support passing `EnsDb` to `annoDb` (#193, #205)
+ export `getAnnoStat()` (#200, #204)

# ChIPseeker 1.33.2

+ supports `by = "ggVennDiagram"` in `vennplot` function (2022-09-13, Tue)

# ChIPseeker 1.33.1

+ `plotPeakProf()` allows passing GRanges object or a list of GRanges objects to TxDb parameter (2022-06-04, Sat)
+ add test files for `getTagMatrix()` and `plotTagMatrix()`
+ `getBioRegion()` supports UTR regions (3'UTR + 5'UTR)
+ `makeBioRegionFromGranges()` supports generating windoes from self-made GRanges object
+ allow specify colors in `covplot()` (2022-05-09, Mon, #185, #188)

# ChIPseeker 1.32.0

+ Bioconductor 3.15 release

# ChIPseeker 1.31.4

+ `readPeakFile` now supports `.broadPeak` and `.gappedPeak` files (2021-12-17, Fri, #173) 

# ChIPseeker 1.31.3

+ bug fixed of determining promoter region in minus strand (2021-12-16, Thu, #172)

# ChIPseeker 1.31.2

+ update vignette

# ChIPseeker 1.31.1

+ bug fixed to take strand information (2021-11-10, Wed, #167)

# ChIPseeker 1.30.0

+ Bioconductor 3.14 release

# ChIPseeker 1.29.2

+ extend functions for plotting peak profiles to support other types of bioregions (2021-10-15, Fri, @MingLi-929, #156, #160, #162, #163)

# ChIPseeker 1.29.1

+ add example for `seq2gene` function (2021-05-21, Fri)

# ChIPseeker 1.28.0

+ Bioconductor 3.13 release (2021-05-20, Thu)

# ChIPseeker 1.27.5

+ update GEO data (103398/1973025 GSM) (2021-05-14, Fri)

# ChIPseeker 1.27.4

+ bug fixed in determine downstream gene (2021-04-27, Thu)
  - <https://github.com/YuLab-SMU/ChIPseeker/pull/148>
+ `getBioRegion` now supports '3UTR' and '5UTR' (2021-03-30, Tue)
  - <https://github.com/YuLab-SMU/ChIPseeker/pull/146>

# ChIPseeker 1.27.3

+ add two parameter, cex and radius, to `plotAnnoPie` (2021-03-12, Fri)
  - <https://github.com/YuLab-SMU/ChIPseeker/pull/144>

# ChIPseeker 1.27.2

+ bug fixed of `getGenomicAnnotation` (2021-03-03, Wed)
  - <https://github.com/YuLab-SMU/ChIPseeker/issues/142>

# ChIPseeker 1.27.1

+ Add support for `EnsDb` annotation databases in `annotatePeak`. 
  - <https://github.com/YuLab-SMU/ChIPseeker/pull/120>

# ChIPseeker 1.26.0

+ Bioconductor 3.12 release (2020-10-28, Wed)


# ChIPseeker 1.23.1

+ update GEO data (51079/762820 GSM) (2019-12-20, Fri)

# ChIPseeker 1.22.0

+ Bioconductor 3.10 release
 
# ChIPseeker 1.21.1

+ new implementation of `upsetplot` (2019-08-29, Thu)
  - use `ggupset`, `ggimage` and `ggplotify`
+ `subset` method for `csAnno` object (2019-08-27, Tue)

# ChIPseeker 1.20.0

+ Bioconductor 3.9 release

# ChIPseeker 1.19.1

+ add `origin_label = "TSS"` parameter to `plotAvgProf` (2018-12-12, Wed)
  - <https://github.com/GuangchuangYu/ChIPseeker/issues/91>
  
# ChIPseeker 1.18.0

+ Bioconductor 3.8 release

# ChIPseeker 1.17.2

+ add `flip_minor_strand` parameter in `getTagMatrix` (2018-08-10, Fri)
  - should set to FALSE if windows if not symetric
  
# ChIPseeker 1.17.1

+ fixed issue of `vennpie` by adding pseudo-count +1 (2018-07-21, Sat)
  - <https://www.biostars.org/p/326456/>

# ChIPseeker 1.16.0

+ Bioconductor 3.7 release

# ChIPseeker 1.15.4

+ If the required input is a named list and user input a list without name,
  set the name automatically and throw warning msg instead of error <2018-03-14,
  Wed>
    - <https://support.bioconductor.org/p/106903/#106936>
+ change `plotAvgProf`'s default y label <2018-03-14, Wed>
    - <https://github.com/GuangchuangYu/ChIPseeker/issues/76>
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
