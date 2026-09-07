##' ChIPseeker: ChIP peak Annotation, Comparison, and Visualization
##'
##' @description
##' ChIPseeker is an R/Bioconductor package designed for comprehensive analysis
##' of ChIP-seq (Chromatin Immunoprecipitation followed by sequencing) data.
##' The package provides a complete toolkit for annotating ChIP peaks, comparing
##' datasets, and visualizing results to facilitate biological interpretation.
##'
##' @details
##' ChIPseeker offers four main categories of functionality:
##'
##' \strong{1. Peak Annotation:}
##' \itemize{
##'   \item \code{\link{annotatePeak}}: Comprehensive peak annotation including
##'         nearest genes, genomic features (promoter, exon, intron, UTR, etc.),
##'         distance to TSS, and optional gene annotation from databases
##'   \item \code{\link{seq2gene}}: Many-to-many mapping of genomic regions to
##'         genes for functional enrichment analysis, considering host genes,
##'         promoter regions, and flanking genes
##'   \item \code{\link{readPeakFile}}: Read peak files in various formats
##'         (BED, narrowPeak, broadPeak, etc.)
##' }
##'
##' \strong{2. Visualization:}
##' \itemize{
##'   \item \code{\link{plotAnnoBar}}, \code{\link{plotAnnoPie}}: Visualize
##'         genomic annotation distribution
##'   \item \code{\link{plotDistToTSS}}: Distribution of peaks relative to
##'         transcription start sites
##'   \item \code{\link{covplot}}: Coverage plot showing peak distribution
##'         across chromosomes
##'   \item \code{\link{plotAvgProf}}, \code{\link{plotPeakProf}}: Average
##'         profiles and peak profiles around TSS or gene body
##'   \item \code{\link{tagHeatmap}}, \code{\link{peakHeatmap}}: Heatmaps of
##'         peak binding patterns
##'   \item \code{\link{vennplot}}, \code{\link{upsetplot}}: Visualize overlap
##'         among multiple peak sets
##' }
##'
##' \strong{3. Peak Comparison and Overlap Analysis:}
##' \itemize{
##'   \item \code{\link{enrichAnnoOverlap}}: Test significance of gene overlap
##'         between ChIP experiments based on nearest gene annotation
##'   \item \code{\link{enrichPeakOverlap}}: Test significance of genomic
##'         coordinate overlap using permutation testing
##'   \item \code{\link{vennplot}}, \code{\link{vennpie}}: Visualize overlaps
##'         among peak sets
##'   \item \code{\link{shuffle}}: Shuffle peak positions for permutation testing
##' }
##'
##' \strong{4. GEO Database Integration:}
##' \itemize{
##'   \item \code{\link{downloadGEObedFiles}}, \code{\link{downloadGSMbedFiles}}:
##'         Download ChIP-seq data from GEO database
##'   \item \code{\link{getGEOInfo}}, \code{\link{getGEOspecies}},
##'         \code{\link{getGEOgenomeVersion}}: Query GEO database metadata
##'   \item Access to over 17,000 ChIP-seq datasets from GEO for comparison
##' }
##'
##' \strong{Additional Features:}
##' \itemize{
##'   \item \code{\link{getTagMatrix}}: Build tag matrices for profile plotting
##'   \item \code{\link{getPromoters}}: Extract promoter regions
##'   \item \code{\link{getBioRegion}}: Define custom genomic regions
##'   \item dplyr verb extensions for GRanges objects: \code{filter},
##'         \code{mutate}, \code{rename}, \code{arrange}
##'   \item Support for liftOver conversion between genome assemblies
##' }
##'
##' @section Main Workflow:
##' The typical ChIPseeker workflow:
##' \enumerate{
##'   \item Read peak files using \code{\link{readPeakFile}}
##'   \item Annotate peaks using \code{\link{annotatePeak}} to get nearest genes,
##'         genomic features, and distances
##'   \item Visualize results using \code{\link{plotAnnoBar}},
##'         \code{\link{plotDistToTSS}}, etc.
##'   \item (Optional) Compare with other datasets using \code{\link{enrichPeakOverlap}}
##'         or download from GEO
##'   \item (Optional) Map regions to genes using \code{\link{seq2gene}} for
##'         functional enrichment analysis
##' }
##'
##' @section Key Classes:
##' \itemize{
##'   \item \code{\link{csAnno-class}}: Container class for annotated peaks
##'         with metadata and statistics
##' }
##'
##' @section Citation:
##' If you use ChIPseeker in published research, please cite:
##' \itemize{
##'   \item Q Wang, M Li, T Wu, L Zhan, L Li, M Chen, W Xie, Z Xie, E Hu,
##'         S Xu, G Yu. Exploring epigenomic datasets by ChIPseeker.
##'         \emph{Current Protocols}, 2022, 2(10): e585.
##'   \item G Yu, LG Wang, QY He. ChIPseeker: an R/Bioconductor package for ChIP
##'         peak annotation, comparison and visualization.
##'         \emph{Bioinformatics}, 2015, 31(14):2382-2383.
##' }
##'
##' @author
##' \strong{Maintainer:} Guangchuang Yu \email{guangchuangyu@gmail.com}
##'
##' Contributors: Ming Li, Qianwen Wang, Yun Yan, Hervé Pagès, Michael Kluge,
##' Thomas Schwarzl, Zhougeng Xu, Chun-Hui Gao
##'
##' @seealso
##' Useful links:
##' \itemize{
##'   \item \url{https://yulab-smu.top/contribution-knowledge-mining/}
##'   \item Report bugs at \url{https://github.com/YuLab-SMU/ChIPseeker/issues}
##'   \item Bioconductor page: \url{https://bioconductor.org/packages/ChIPseeker}
##' }
##'
##' @keywords package
##' @docType package
##' @name ChIPseeker-package
##' @aliases ChIPseeker ChIPseeker-package
"_PACKAGE"



##' Information Datasets
##'
##' Pre-calculated datasets included in the ChIPseeker package:
##' \itemize{
##'   \item \code{ucsc_release}: UCSC genome version information
##'   \item \code{gsminfo}: GEO Sample (GSM) information for ChIP-seq datasets
##'   \item \code{tagMatrixList}: Pre-calculated tag matrices for example datasets
##' }
##'
##' @name info
##' @aliases ucsc_release gsminfo tagMatrixList
##' @docType data
##' @keywords datasets
NULL

##' Name of the ChIPseeker cache environment (internal static variable)
##'
##' This variable stores the name of the cache environment used internally
##' by ChIPseeker to store annotation data (exons, introns, UTRs, etc.) for
##' efficient repeated access.
##'
##' @format character vector of length 1
##' @keywords internal
ChIPseekerCache <- "ChIPseekerEnv"

