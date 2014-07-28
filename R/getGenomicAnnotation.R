updateGenomicAnnotation <- function(peaks, genomicRegion, type, annotation) {
    hits <- getGenomicAnnotation.internal(peaks, genomicRegion, type)
    if (length(hits) == 2) {
        annotation[hits$queryIndex] <- hits$annotation
    }
    return(annotation)
 }         


##' get Genomic Annotation of peaks
##'
##' 
##' @title getGenomicAnnotation
##' @param peaks peaks in GRanges object
##' @param distance distance of peak to TSS
##' @param tssRegion tssRegion, default is -3kb to +3kb
##' @param TxDb TxDb object
##' @importFrom GenomicFeatures intronsByTranscript
##' @importFrom GenomicFeatures threeUTRsByTranscript
##' @importFrom GenomicFeatures fiveUTRsByTranscript
##' @importFrom GenomicFeatures exonsBy
##' @return character vector
##' @author G Yu
getGenomicAnnotation <- function(peaks,
                                 distance,
                                 tssRegion=c(-3000, 3000),
                                 TxDb
                                 ) {
    
    ##
    ## since some annotation overlap,
    ## a priority is assign based on the following:
    ## 1. Promoter
    ## 2. 5' UTR
    ## 3. 3' UTR
    ## 4. Exon
    ## 5. Intron
    ## 6. Downstream
    ## 7. Intergenic
    ##


    .ChIPseekerEnv(TxDb)
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    
            
    annotation <- rep(NA, length(distance))
    ## Intergenic
    annotation[is.na(annotation)] <- "Intergenic"
    
    ## Introns
    if ( exists("intronList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        intronList <- get("intronList", envir=ChIPseekerEnv)
    } else {
        intronList <- intronsByTranscript(TxDb)
        assign("intronList", intronList, envir=ChIPseekerEnv)
    }
    annotation <- updateGenomicAnnotation(peaks, intronList, "Intron", annotation)

    ## Exon
    if ( exists("exonList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        exonList <- get("exonList", envir=ChIPseekerEnv)
    } else {
        exonList <- exonsBy(TxDb)
        assign("exonList", exonList, envir=ChIPseekerEnv)
    }
    annotation <- updateGenomicAnnotation(peaks, exonList, "Exon", annotation)

    ## 3' UTR Exons
    if ( exists("threeUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        threeUTRList <- get("threeUTRList", envir=ChIPseekerEnv)
    } else {
        threeUTRList <- threeUTRsByTranscript(TxDb)
        assign("threeUTRList", threeUTRList, envir=ChIPseekerEnv)
    }
    annotation <- updateGenomicAnnotation(peaks, threeUTRList, "3' UTR", annotation)
    
    ## 5' UTR Exons
    if ( exists("fiveUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        fiveUTRList <- get("fiveUTRList", envir=ChIPseekerEnv)
    } else {
        fiveUTRList <- fiveUTRsByTranscript(TxDb)
        assign("fiveUTRList", fiveUTRList, envir=ChIPseekerEnv)
    }
    annotation <- updateGenomicAnnotation(peaks, fiveUTRList, "5' UTR", annotation)
    
    ## TSS
    annotation[distance >= tssRegion[1] &
               distance <= tssRegion[2] ] <- "Promoter"

    pm <- max(abs(tssRegion))
    if (pm/1000 >= 2) {
        dd <- seq(1:ceiling(pm/1000))*1000
        for (i in 1:length(dd)) {
            if (i == 1) {
                lbs <- paste("Promoter", " (<=", dd[i]/1000, "kb)", sep="")
                annotation[abs(distance) <= dd[i] &
                           annotation == "Promoter"] <- lbs
            } else {
                lbs <- paste("Promoter", " (", dd[i-1]/1000, "-", dd[i]/1000, "kb)", sep="")
                annotation[abs(distance) <= dd[i] &
                           abs(distance) > dd[i-1] &
                           annotation == "Promoter"] <- lbs
            }
        }
    }

    features <- getGene(TxDb, by="gene")

    ## nearest from gene end
    idx <- follow(peaks, features)
    peF <- features[idx]
    dd <- ifelse(strand(peF) == "+",
                 start(peaks) - end(peF),
                 end(peaks) - start(peF))
    for (i in 1:3) {
        j <- which(annotation == "Intergenic" & abs(dd) <= i*1000)
        if (length(j) > 0) {
            if (i == 1) {
                lbs <- "Downstream (<1kb)"
            } else {
                lbs <- paste("Downstream (", i-1, "-", i, "kb)", sep="")
            }
            annotation[j] <- lbs
        }
    }
    annotation[annotation == "Intergenic"] = "Distal Intergenic"
    return(annotation)
}


##' @importFrom IRanges elementLengths
##' @importFrom IRanges findOverlaps
## @importFrom IRanges queryHits
## @importFrom IRanges subjectHits
##' @importFrom S4Vectors queryHits
##' @importFrom S4Vectors subjectHits
##' @importMethodsFrom BiocGenerics unlist
getGenomicAnnotation.internal <- function(peaks, genomicRegion, type){
    GRegion <- unlist(genomicRegion)
    GRegionLen <- elementLengths(genomicRegion)
    if (type == "Intron" || type =="Exon") {
        nn <- TXID2EG(names(genomicRegion))
        names(GRegionLen) <- nn
        GRegion$gene_id <- rep(nn, times=GRegionLen)
    } else {
        names(GRegionLen) <- names(genomicRegion)
        GRegion$gene_id <- rep(names(genomicRegion), times=GRegionLen)
    }

    if (type == "Intron") {
        intron_rank <- unlist(sapply(GRegionLen, function(i) seq(0, i)))
        intron_rank <- intron_rank[intron_rank != 0]
        GRegion$intron_rank <- intron_rank
    }
    ## find overlap
    GRegionHit <- findOverlaps(peaks, GRegion)
    if (length(GRegionHit) == 0) {
        return(NA)
    }
    qh <- queryHits(GRegionHit)
    hit.idx <- getFirstHitIndex(qh)
    GRegionHit <- GRegionHit[hit.idx]
    queryIndex <- queryHits(GRegionHit)
    subjectIndex <- subjectHits(GRegionHit)

    hits <- GRegion[subjectIndex]
    geneID <- hits$gene_id

    if (type == "Intron") {
        anno <- paste(type, " (", geneID, ", intron ", hits$intron_rank,
                      " of ", GRegionLen[geneID], ")", sep="")
    } else if (type == "Exon") {
        anno <- paste(type, " (", geneID, ", exon ", hits$exon_rank,
                      " of ", GRegionLen[geneID], ")", sep="")
    } else {
        anno <- type
    }
    res <- list(queryIndex=queryIndex, annotation=anno)
    return(res)
}
