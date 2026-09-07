##' Annotate genomic regions to genes in many-to-many mapping
##'
##' This function associates genomic regions (e.g., ChIP-seq peaks, enhancers,
##' regulatory elements) with coding genes in a many-to-many mapping. It is
##' designed to link both coding and non-coding genomic regions to coding genes
##' to facilitate functional enrichment analysis.
##'
##' @description
##' The function maps genomic regions to genes through three mechanisms:
##' \enumerate{
##'   \item \strong{Host genes}: Regions that overlap with exons or introns of
##'         genes. These are genes that directly contain the genomic regions
##'   \item \strong{Proximal genes}: Regions located within the promoter region
##'         (TSS region) of genes. These are genes whose transcription start sites
##'         are near the genomic regions
##'   \item \strong{Flanking genes}: Regions located within a specified distance
##'         from genes but outside the promoter region. These are nearby genes that
##'         may be regulated via cis-regulatory mechanisms
##' }
##'
##' The function returns a vector of unique gene IDs that can be used for downstream
##' functional enrichment analysis (e.g., GO, KEGG, Reactome pathway analysis).
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Identifies host genes by finding overlaps with exons and introns using
##'         \code{getGenomicAnnotation.internal()}
##'   \item Extracts gene IDs from overlapping exons and introns, converting
##'         transcript-based IDs to gene IDs
##'   \item Finds the nearest gene to each region using
##'         \code{getNearestFeatureIndicesAndDistances()}
##'   \item Categorizes nearest genes into:
##'     \itemize{
##'       \item Promoter genes: genes where the distance falls within \code{tssRegion}
##'       \item Flanking genes: genes where the absolute distance is less than
##'             \code{flankDistance} but outside \code{tssRegion}
##'     }
##'   \item Combines all gene IDs (host, promoter, and flanking) and returns unique
##'         gene IDs
##' }
##'
##' This many-to-many mapping allows a single genomic region to be associated with
##' multiple genes (e.g., overlapping multiple genes or having multiple nearby genes),
##' and a single gene can be associated with multiple genomic regions.
##'
##' @param seq GRanges object containing genomic regions to be mapped to genes.
##'   Typically ChIP-seq peaks, enhancers, or other regulatory elements
##' @param tssRegion numeric vector of length 2 specifying the promoter region
##'   around the transcription start site (TSS). Default example: c(-1000, 1000)
##'   means 1kb upstream and 1kb downstream of TSS. Negative values indicate
##'   upstream, positive values indicate downstream. Regions within this range
##'   are considered to be in promoter regions
##' @param flankDistance numeric, distance in base pairs to search for flanking
##'   genes. Regions that are within this distance from a gene (but outside the
##'   \code{tssRegion}) are considered flanking genes. Default example: 3000 (3kb)
##' @param TxDb TxDb or EnsDb annotation object containing gene/transcript
##'   structure information
##' @param sameStrand logical, whether to only consider genes on the same strand
##'   as the genomic regions when finding overlaps and nearest features. If FALSE,
##'   strand is ignored. Default is FALSE
##' @return Character vector of unique gene IDs. The vector contains gene IDs from:
##'   \itemize{
##'     \item Host genes (regions overlapping exons/introns)
##'     \item Proximal genes (regions within TSS region)
##'     \item Flanking genes (regions within \code{flankDistance} but outside TSS region)
##'   }
##'   The gene IDs are typically Entrez IDs (for TxDb) or gene_id (for EnsDb),
##'   and can be used directly for functional enrichment analysis with packages
##'   like \code{clusterProfiler}, \code{DOSE}, or \code{ReactomePA}
##' @export
##' @examples
##' \dontrun{
##' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' file <- getSampleFiles()[[1]] # a bed file
##' gr <- readPeakFile(file)
##' genes <- seq2gene(gr, tssRegion=c(-1000, 1000), flankDistance = 3000, TxDb)
##'
##' ## Use the gene IDs for enrichment analysis
##' ## library(clusterProfiler)
##' ## ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "BP")
##' }
##' @importFrom yulab.utils get_cache_element
##' @importFrom yulab.utils update_cache_item
##' @seealso \code{\link{annotatePeak}} for detailed peak annotation,
##'   \code{\link{getNearestFeatureIndicesAndDistances}} for finding nearest features
##' @author Guangchuang Yu
seq2gene <- function(seq, tssRegion, flankDistance, TxDb, sameStrand=FALSE) {
    .ChIPseekerEnv(TxDb, item = ChIPseekerCache)
    # ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

    ## get genes whose exons or introns overlap with the genomic regions
    ## Exons
    exonList <- get_cache_element(item = ChIPseekerCache, elements = "exonList")
    if(is.null(exonList)){
        exonList <- exonsBy(TxDb)
        update_cache_item(item = ChIPseekerCache, list("exonList" = exonList))
    }

    # if ( exists("exonList", envir=ChIPseekerEnv, inherits=FALSE) ) {
    #     exonList <- get("exonList", envir=ChIPseekerEnv)
    # } else {
    #     exonList <- exonsBy(TxDb)
    #     assign("exonList", exonList, envir=ChIPseekerEnv)
    # }
    exons <- getGenomicAnnotation.internal(seq, exonList, type = "Exon", sameStrand=sameStrand)

    ## Introns
    intronList <- get_cache_element(item = ChIPseekerCache, elements = "intronList")

    if(is.null(intronList)){
        intronList <- intronsByTranscript(TxDb)
        update_cache_item(item = ChIPseekerCache, list("intronList" = intronList))
    }

    # if ( exists("intronList", envir=ChIPseekerEnv, inherits=FALSE) ) {
    #     intronList <- get("intronList", envir=ChIPseekerEnv)
    # } else {
    #     intronList <- intronsByTranscript(TxDb)
    #     assign("intronList", intronList, envir=ChIPseekerEnv)
    # }
    introns <- getGenomicAnnotation.internal(seq, intronList, type="Intron", sameStrand=sameStrand)

    genes <- c(exons$gene, introns$gene)
    ## > head(genes)
    ## [1] "uc001aed.3/126789"    "uc001aka.3/440556"    "uc001ako.3/49856"
    ## [4] "uc001alg.3/100133612" "uc009vly.2/390992"    "uc001awv.2/79814"
    genes <- gsub("\\w+\\.*\\d*/(\\d+)", "\\1", genes)
    ## > head(genes)
    ## [1] "126789"    "440556"    "49856"     "100133612" "390992"    "79814"

    ### find the nearest gene to the peaks
    features <- getGene(TxDb, by="gene")
    idx.dist <- getNearestFeatureIndicesAndDistances(seq, features, sameStrand=sameStrand)
    nearestFeatures <- features[idx.dist$index]

    distance <- idx.dist$distance
    # pi <- distance > tssRegion[1] & distance < tssRegion[2]
    # promoters <- mcols(nearestFeatures[pi])[["gene_id"]]
    # nearest_genes <- mcols(nearestFeatures[!pi][abs(distance[!pi]) < flankDistance])[["gene_id"]]
    # genes <- c(genes, promoters, nearest_genes)

    # The separate promoter extraction is redundant. Since flankDistance is typically larger
    # than the TSS region (e.g., tssRegion = c(-1000, 1000) = 2kb, flankDistance = 3000bp),
    # all promoters are already within flankDistance.
    # get genes within flankDistance (includes both promoter and flanking genes)
    nearest_genes <- mcols(nearestFeatures[abs(distance) < flankDistance |
                          (distance > tssRegion[1] & distance < tssRegion[2])])[["gene_id"]]
    genes <- c(genes, nearest_genes)
    return(unique(genes))
}
