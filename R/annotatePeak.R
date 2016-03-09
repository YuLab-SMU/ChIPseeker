##' Annotate peaks
##'
##' 
##' @title annotatePeak
##' @param peak peak file or GRanges object
##' @param tssRegion Region Range of TSS 
##' @param TxDb TxDb object
##' @param level one of transcript and gene
##' @param assignGenomicAnnotation logical, assign peak genomic annotation or not
##' @param genomicAnnotationPriority genomic annotation priority
##' @param annoDb annotation package
##' @param addFlankGeneInfo logical, add flanking gene information from the peaks 
##' @param flankDistance distance of flanking sequence
##' @param sameStrand logical, whether find nearest/overlap gene in the same strand
##' @param ignoreOverlap logical, whether ignore overlap of TSS with peak
##' @param ignoreUpstream logical, if True only annotate gene at the 3' of the peak.
##' @param ignoreDownstream logical, if True only annotate gene at the 5' of the peak.
##' @param overlap one of 'TSS' or 'all', if overlap="all", then gene overlap with peak will be reported as nearest gene, no matter the overlap is at TSS region or not.
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
##' @import BiocGenerics S4Vectors GenomeInfoDb 
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 3000), TxDb=txdb)
##' peakAnno
##' }
##' @seealso \code{\link{plotAnnoBar}} \code{\link{plotAnnoPie}} \code{\link{plotDistToTSS}}
##' @export
##' @author G Yu 
annotatePeak <- function(peak,
                         tssRegion=c(-3000, 3000),
                         TxDb=NULL,
                         level = "transcript",
                         assignGenomicAnnotation=TRUE,
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         annoDb=NULL,
                         addFlankGeneInfo=FALSE,
                         flankDistance=5000,
                         sameStrand = FALSE,
                         ignoreOverlap=FALSE,
                         ignoreUpstream=FALSE,
                         ignoreDownstream=FALSE,
                         overlap = "TSS",
                         verbose=TRUE) {

    is_GRanges_of_TxDb <- FALSE
    if (is(TxDb, "GRanges")) {
        is_GRanges_of_TxDb <- TRUE
        assignGenomicAnnotation <- FALSE
        annoDb <- NULL
        addFlankGeneInfo <- FALSE
        message("#\n#.. 'TxDb' is a self-defined 'GRanges' object...\n#")
        message("#.. Some parameters of 'annotatePeak' will be disable,")
        message("#.. including:")
        message("#..\tlevel, assignGenomicAnnotation, genomicAnnotationPriority,")
        message("#..\tannoDb, addFlankGeneInfo and flankDistance.")
        message("#\n#.. Some plotting functions are designed for visualizing genomic annotation")
        message("#.. and will not be available for the output object.\n#")
    }

    if (is_GRanges_of_TxDb) {
        level <- "USER_DEFINED"
    } else {
        level <- match.arg(level, c("transcript", "gene"))
    }

    if (assignGenomicAnnotation && all(genomicAnnotationPriority %in% c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")) == FALSE) {
        stop('genomicAnnotationPriority should be any order of c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")')
    }
    
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

    if (is_GRanges_of_TxDb) {
        features <- TxDb
    } else {
        TxDb <- loadTxDb(TxDb)
        
        if (level=="transcript") {
            features <- getGene(TxDb, by="transcript")
        } else {
            features <- getGene(TxDb, by="gene")
        }
    } 
    if (verbose)
        cat(">> identifying nearest features...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    
    ## nearest features
    idx.dist <- getNearestFeatureIndicesAndDistances(peak.gr, features,
                                                     sameStrand, ignoreOverlap,
                                                     ignoreUpstream,ignoreDownstream,
                                                     overlap=overlap)
    nearestFeatures <- features[idx.dist$index]
    if (verbose)
        cat(">> calculating distance from peak to TSS...\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    ## distance
    distance <- idx.dist$distance

    ## update peak, remove un-map peak if exists.
    peak.gr <- idx.dist$peak
    
    ## annotation
    if (assignGenomicAnnotation == TRUE) {
        if (verbose)
            cat(">> assigning genomic annotation...\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
        
        anno <- getGenomicAnnotation(peak.gr, distance, tssRegion, TxDb, level, genomicAnnotationPriority, sameStrand=sameStrand)
        annotation <- anno[["annotation"]]
        detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]
    } else {
        annotation <- NULL
        detailGenomicAnnotation <- NULL
    }
    ## duplicated names since more than 1 peak may annotated by only 1 gene
    names(nearestFeatures) <- NULL
    nearestFeatures.df <- as.data.frame(nearestFeatures)
    if (is_GRanges_of_TxDb) {
        colnames(nearestFeatures.df)[1:5] <- c("geneChr", "geneStart", "geneEnd",
                                          "geneLength", "geneStrand")
    } else if (level == "transcript") {
        colnames(nearestFeatures.df) <- c("geneChr", "geneStart", "geneEnd",
                                          "geneLength", "geneStrand", "geneId", "transcriptId")
        nearestFeatures.df$geneId <- TXID2EG(as.character(nearestFeatures.df$geneId), geneIdOnly=TRUE)
    } else {
        colnames(nearestFeatures.df) <- c("geneChr", "geneStart", "geneEnd",
                                          "geneLength", "geneStrand", "geneId")
    }

    ## append annotation to peak.gr
    if (!is.null(annotation))
        mcols(peak.gr)[["annotation"]] <- annotation

    for(cn in colnames(nearestFeatures.df)) {
        mcols(peak.gr)[[cn]] <- unlist(nearestFeatures.df[, cn])
    }

    mcols(peak.gr)[["distanceToTSS"]] <- distance
 
    if (!is.null(annoDb)) {
        if (verbose)
            cat(">> adding gene annotation...\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
        .idtype <- IDType(TxDb)
        if (length(.idtype) == 0 || is.na(.idtype) || is.null(.idtype)) {
            if (grepl('ENSG', peak.gr$geneId[1])) {
                .idtype <- "Ensemble Gene ID"
            } else if (grepl('^\\d+$',, peak.gr$geneId[1])) {
                .idtype <- "Entrez Gene ID"                
            } else {
                warning("Unknown ID type, gene annotation will not be added...")
                .idtype <- NA
            }
        }

        if (!is.na(.idtype)) {
            peak.gr %<>% addGeneAnno(annoDb, IDType(TxDb))
        }
    }
    
    if (addFlankGeneInfo == TRUE) {
        if (verbose)
            cat(">> adding flank feature information from peaks...\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
 
        flankInfo <- getAllFlankingGene(peak.gr, features, level, flankDistance)

        if (level == "transcript") {
            mcols(peak.gr)[["flank_txIds"]] <- NA
            mcols(peak.gr)[["flank_txIds"]][flankInfo$peakIdx] <- flankInfo$flank_txIds
        }
        
        mcols(peak.gr)[["flank_geneIds"]] <- NA
        mcols(peak.gr)[["flank_gene_distances"]] <- NA
        
        mcols(peak.gr)[["flank_geneIds"]][flankInfo$peakIdx] <- flankInfo$flank_geneIds
        mcols(peak.gr)[["flank_gene_distances"]][flankInfo$peakIdx] <- flankInfo$flank_gene_distances

    }

    if (!is_GRanges_of_TxDb) {
        if(verbose)
            cat(">> assigning chromosome lengths\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
        
        peak.gr@seqinfo <- seqinfo(TxDb)[names(seqlengths(peak.gr))]
    }
    
    if(verbose)
        cat(">> done...\t\t\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")

    if (assignGenomicAnnotation) {
        res <- new("csAnno",
                   anno = peak.gr,
                   tssRegion = tssRegion,
                   level=level,
                   hasGenomicAnnotation = TRUE,
                   detailGenomicAnnotation=detailGenomicAnnotation,
                   annoStat=getGenomicAnnoStat(peak.gr),
                   peakNum=peakNum
                   )
    } else {
        res <- new("csAnno",
                   anno = peak.gr,
                   tssRegion = tssRegion,
                   level=level,
                   hasGenomicAnnotation = FALSE,
                   peakNum=peakNum
                   )  
    }
    return(res)
}


##' dropAnno
##'
##' drop annotation exceeding distanceToTSS_cutoff
##' @title dropAnno
##' @param csAnno output of annotatePeak
##' @param distanceToTSS_cutoff distance to TSS cutoff
##' @return csAnno object
##' @export
##' @author Guangchuang Yu
dropAnno <- function(csAnno, distanceToTSS_cutoff=10000) {
    idx <- which(abs(mcols(csAnno@anno)[["distanceToTSS"]]) < distanceToTSS_cutoff)
    csAnno@anno <- csAnno@anno[idx]
    csAnno@peakNum <- length(idx)
    if (csAnno@hasGenomicAnnotation) {
        csAnno@annoStat <- getGenomicAnnoStat(csAnno@anno)
        detailGenomicAnnotation = detailGenomicAnnotation[idx,]
    }
    csAnno
}
