##' annotate genomic regions to genes in many-to-many mapping
##'
##' 
##' @title seq2gene
##' @param seq genomic regions in GRanges object
##' @param tssRegion TSS region
##' @param flankDistance flanking search radius
##' @param TxDb TranscriptDb object
##' @param sameStrand logical whether find nearest/overlap gene in the same strand
##' @return gene vector
##' @export
##' @examples
##' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' file <- getSampleFiles()[[1]] # a bed file
##' gr <- readPeakFile(file)
##' genes <- seq2gene(gr, tssRegion=c(-1000, 1000), flankDistance = 3000, TxDb) 
##' @author Guangchuang Yu
seq2gene <- function(seq, tssRegion, flankDistance, TxDb, sameStrand=FALSE) {
    .ChIPseekerEnv(TxDb)
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    
    ## Exons
    if ( exists("exonList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        exonList <- get("exonList", envir=ChIPseekerEnv)
    } else {
        exonList <- exonsBy(TxDb)
        assign("exonList", exonList, envir=ChIPseekerEnv)
    }
    exons <- getGenomicAnnotation.internal(seq, exonList, type = "Exon", sameStrand=sameStrand)
    
    ## Introns
    if ( exists("intronList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        intronList <- get("intronList", envir=ChIPseekerEnv)
    } else {
        intronList <- intronsByTranscript(TxDb)
        assign("intronList", intronList, envir=ChIPseekerEnv)
    }
    introns <- getGenomicAnnotation.internal(seq, intronList, type="Intron", sameStrand=sameStrand)
    
    genes <- c(exons$gene, introns$gene)
    ## > head(genes)
    ## [1] "uc001aed.3/126789"    "uc001aka.3/440556"    "uc001ako.3/49856"    
    ## [4] "uc001alg.3/100133612" "uc009vly.2/390992"    "uc001awv.2/79814"   
    genes <- gsub("\\w+\\.*\\d*/(\\d+)", "\\1", genes)
    ## > head(genes)
    ## [1] "126789"    "440556"    "49856"     "100133612" "390992"    "79814"   

    features <- getGene(TxDb, by="gene")
    idx.dist <- getNearestFeatureIndicesAndDistances(seq, features, sameStrand=sameStrand)
    nearestFeatures <- features[idx.dist$index] 
    
    distance <- idx.dist$distance

    pi <- distance > tssRegion[1] & distance < tssRegion[2]
    promoters <- mcols(nearestFeatures[pi])[["gene_id"]]

    nearest_genes <- mcols(nearestFeatures[!pi][abs(distance[!pi]) < flankDistance])[["gene_id"]]

    genes <- c(genes, promoters, nearest_genes)
    return(unique(genes))
}
