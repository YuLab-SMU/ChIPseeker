##' Annotate peaks
##'
##' 
##' @title annotatePeak
##' @param peak peak file or GRanges object
##' @param tssRegion Region Range of TSS 
##' @param TxDb TxDb object
##' @param level one of transcript and gene
##' @param assignGenomicAnnotation logical, assign peak genomic annotation or not
##' @param annoDb annotation package
##' @param addFlankGeneInfo logical, add flanking gene information from the peaks 
##' @param flankDistance distance of flanking sequence
##' @param verbose print message or not
##' @return data.frame or GRanges object with columns of:
##' 
##' all columns provided by input.
##' 
##' annotation: genomic feature of the peak, for instance if the peak is
##' located in 5'UTR, it will annotated by 5'UTR. Possible annotation is
##' Promoter-TSS, Exon, 5' UTR, 3' UTR, Intron, and Intergenic.
##' 
##' geneChr: Chromosome of the nearest gene
##' 
##' geneStart: gene start
##' 
##' geneEnd: gene end
##' 
##' geneLength: gene length
##' 
##' geneStrand: gene strand
##' 
##' geneId: entrezgene ID
##' 
##' distanceToTSS: distance from peak to gene TSS
##' 
##' if annoDb is provided, extra column will be included:
##' 
##' ENSEMBL: ensembl ID of the nearest gene
##' 
##' SYMBOL: gene symbol
##' 
##' GENENAME: full gene name
##' @importFrom S4Vectors metadata
##' @importFrom S4Vectors mcols
##' @importFrom S4Vectors mcols<-
##' @importFrom GenomeInfoDb seqlengths
##' @importFrom GenomeInfoDb seqinfo
##' @importFrom GenomeInfoDb seqinfo<-
## @importFrom GenomicFeatures getChromInfoFromUCSC
##' @importMethodsFrom BiocGenerics as.data.frame
##' @examples
##' require(TxDb.Hsapiens.UCSC.hg38.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 3000), TxDb=txdb)
##' peakAnno
##' @seealso \code{\link{plotAnnoBar}} \code{\link{plotAnnoPie}} \code{\link{plotDistToTSS}}
##' @export
##' @author G Yu 
annotatePeak <- function(peak,
                         tssRegion=c(-3000, 3000),
                         TxDb=NULL,
                         level = "transcript",
                         assignGenomicAnnotation=TRUE,
                         annoDb=NULL,
                         addFlankGeneInfo=FALSE,
                         flankDistance=5000,
                         verbose=TRUE) {
    
    level <- match.arg(level, c("transcript", "gene"))
    
    if ( is(peak, "GRanges") ){
        ## this test will be TRUE
        ## when peak is an instance of class/subclass of "GRanges"
        input <- "gr"
        peak.gr <- peak
    } else {
        input <- "file"
        peak.gr <- loadPeak(peak, verbose)
    }

    peakNum <- length(peak.gr)
    
    if (verbose)
        cat(">> preparing features information...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")

    TxDb <- loadTxDb(TxDb)

    if (level=="transcript") {
        features <- getGene(TxDb, by="transcript")
    } else {
         features <- getGene(TxDb, by="gene")
    }
    if (verbose)
        cat(">> identifying nearest features...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    ## nearest features
    idx.dist <- getNearestFeatureIndicesAndDistances(peak.gr, features)
    nearestFeatures <- features[idx.dist$index]
    if (verbose)
        cat(">> calculating distance from peak to TSS...\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    ## distance
    distance <- idx.dist$distance

    ## update peak, remove un-map peak if exists.
    peak.gr <- idx.dist$peak

    if (verbose)
        cat(">> assigning genomic annotation...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    ## annotation
    if (assignGenomicAnnotation == TRUE) {
        anno <- getGenomicAnnotation(peak.gr, distance, tssRegion, TxDb, level)
        annotation <- anno[["annotation"]]
        detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]
    } else {
        annotation <- NULL
        detailGenomicAnnotation <- NULL
    }
    ## duplicated names since more than 1 peak may annotated by only 1 gene
    names(nearestFeatures) <- NULL
    nearestFeatures.df <- as.data.frame(nearestFeatures)
    if (level == "transcript") {
        colnames(nearestFeatures.df) <- c("geneChr", "geneStart", "geneEnd",
                                          "geneLength", "geneStrand", "geneId", "transcriptId")
        nearestFeatures.df$geneId <- TXID2EG(as.character(nearestFeatures.df$geneId), geneIdOnly=TRUE)
    } else {
        colnames(nearestFeatures.df) <- c("geneChr", "geneStart", "geneEnd",
                                          "geneLength", "geneStrand", "geneId")
    }

    ## append annotation to peak.gr
    if (!is.null(annotation))
        ## elementMetadata(peak.gr)[["annotation"]] <- annotation
	mcols(peak.gr)[["annotation"]] <- annotation

    for(cn in colnames(nearestFeatures.df)) {
        ## elementMetadata(peak.gr)[[cn]] <- unlist(nearestFeatures.df[, cn])
    	mcols(peak.gr)[[cn]] <- unlist(nearestFeatures.df[, cn])
    }

    ## elementMetadata(peak.gr)[["distanceToTSS"]] <- distance
    mcols(peak.gr)[["distanceToTSS"]] <- distance
 
    if (!is.null(annoDb)) {
        if (verbose)
            cat(">> adding gene annotation...\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
        IDType <- metadata(TxDb)[8,2]     
        geneAnno <- addGeneAnno(annoDb, peak.gr$geneId, type=IDType)
        if (! all(is.na(geneAnno))) {
            for(cn in colnames(geneAnno)[-1]) {
                ## elementMetadata(peak.gr)[[cn]] <- geneAnno[, cn]
		mcols(peak.gr)[[cn]] <- geneAnno[, cn]
	    }
        }
    }

    if (addFlankGeneInfo == TRUE) {
        if (verbose)
            cat(">> adding flank feature information from peaks...\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
 
        flankInfo <- getAllFlankingGene(peak.gr, features, flankDistance)
        ## elementMetadata(peak.gr)[["flank_txIds"]] <- NA
        ## elementMetadata(peak.gr)[["flank_geneIds"]] <- NA
        ## elementMetadata(peak.gr)[["flank_gene_distances"]] <- NA
        
        ## elementMetadata(peak.gr)[["flank_txIds"]][flankInfo$peakIdx] <- flankInfo$flank_txIds
        ## elementMetadata(peak.gr)[["flank_geneIds"]][flankInfo$peakIdx] <- flankInfo$flank_geneIds
        ## elementMetadata(peak.gr)[["flank_gene_distances"]][flankInfo$peakIdx] <- flankInfo$flank_gene_distances

        mcols(peak.gr)[["flank_txIds"]] <- NA
        mcols(peak.gr)[["flank_geneIds"]] <- NA
        mcols(peak.gr)[["flank_gene_distances"]] <- NA
        
        mcols(peak.gr)[["flank_txIds"]][flankInfo$peakIdx] <- flankInfo$flank_txIds
        mcols(peak.gr)[["flank_geneIds"]][flankInfo$peakIdx] <- flankInfo$flank_geneIds
        mcols(peak.gr)[["flank_gene_distances"]][flankInfo$peakIdx] <- flankInfo$flank_gene_distances

    }
    
    if(verbose)
        cat(">> assigning chromosome lengths\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")

    seqinfo(peak.gr) <- seqinfo(TxDb)[names(seqlengths(peak.gr))]
    
    if(verbose)
        cat(">> done...\t\t\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")

    res <- new("csAnno",
               anno = peak.gr,
               tssRegion = tssRegion,
               level=level,
               detailGenomicAnnotation=detailGenomicAnnotation,
               annoStat=getGenomicAnnoStat(peak.gr),
               peakNum=peakNum
               )
    return(res)
}


