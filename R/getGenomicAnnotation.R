updateGenomicAnnotation <- function(peaks, genomicRegion, type, anno) {
    hits <- getGenomicAnnotation.internal(peaks, genomicRegion, type)
    if (length(hits) == 2) {
        hitIndex <- hits$queryIndex
        anno[["annotation"]][hitIndex] <- hits$annotation
        anno[["detailGenomicAnnotation"]][hitIndex, type] <- TRUE
    }
    return(anno)
 }         


##' get Genomic Annotation of peaks
##'
##' 
##' @title getGenomicAnnotation
##' @param peaks peaks in GRanges object
##' @param distance distance of peak to TSS
##' @param tssRegion tssRegion, default is -3kb to +3kb
##' @param TxDb TxDb object
##' @param level one of gene or transcript
##' @importFrom GenomicFeatures intronsByTranscript
##' @importFrom GenomicFeatures threeUTRsByTranscript
##' @importFrom GenomicFeatures fiveUTRsByTranscript
##' @importFrom GenomicFeatures exonsBy
##' @return character vector
##' @author G Yu
getGenomicAnnotation <- function(peaks,
                                 distance,
                                 tssRegion=c(-3000, 3000),
                                 TxDb,
                                 level
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

    flag <- rep(FALSE, length(distance))
    detailGenomicAnnotation <- data.frame(
        genic=flag,
        Intergenic=flag,
        Promoter=flag,
        fiveUTR=flag,
        threeUTR=flag,
        Exon=flag,
        Intron=flag,
        downstream=flag,
        distal_intergenic=flag)
        
    ## Intergenic
    annotation[is.na(annotation)] <- "Intergenic"

    anno <- list(annotation=annotation,
                 detailGenomicAnnotation=detailGenomicAnnotation)
    
    ## Introns
    if ( exists("intronList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        intronList <- get("intronList", envir=ChIPseekerEnv)
    } else {
        intronList <- intronsByTranscript(TxDb)
        assign("intronList", intronList, envir=ChIPseekerEnv)
    }
    anno <- updateGenomicAnnotation(peaks, intronList, "Intron", anno)

    ## Exon
    if ( exists("exonList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        exonList <- get("exonList", envir=ChIPseekerEnv)
    } else {
        exonList <- exonsBy(TxDb)
        assign("exonList", exonList, envir=ChIPseekerEnv)
    }
    anno <- updateGenomicAnnotation(peaks, exonList, "Exon", anno)

    intergenicIndex <- which(anno[["annotation"]] == "Intergenic")
    anno[["detailGenomicAnnotation"]][intergenicIndex, "Intergenic"] <- TRUE
    anno[["detailGenomicAnnotation"]][-intergenicIndex, "genic"] <- TRUE
    
    ## 3' UTR Exons
    if ( exists("threeUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        threeUTRList <- get("threeUTRList", envir=ChIPseekerEnv)
    } else {
        threeUTRList <- threeUTRsByTranscript(TxDb)
        assign("threeUTRList", threeUTRList, envir=ChIPseekerEnv)
    }
    anno <- updateGenomicAnnotation(peaks, threeUTRList, "threeUTR", anno)
    
    ## 5' UTR Exons
    if ( exists("fiveUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        fiveUTRList <- get("fiveUTRList", envir=ChIPseekerEnv)
    } else {
        fiveUTRList <- fiveUTRsByTranscript(TxDb)
        assign("fiveUTRList", fiveUTRList, envir=ChIPseekerEnv)
    }
    anno <- updateGenomicAnnotation(peaks, fiveUTRList, "fiveUTR", anno)

    annotation <- anno[["annotation"]]
    detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]

    ## TSS
    tssIndex <- distance >= tssRegion[1] & distance <= tssRegion[2] 
    annotation[tssIndex] <- "Promoter"
    detailGenomicAnnotation[tssIndex, "Promoter"] <- TRUE
    
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

 
    features <- getGene(TxDb, by=level)

    ## nearest from gene end
    idx <- follow(peaks, features)
    na.idx <- which(is.na(idx))
    peF <- features[idx[-na.idx]]
    dd <- ifelse(strand(peF) == "+",
                 start(peaks[-na.idx]) - end(peF),
                 end(peaks[-na.idx]) - start(peF))
    dd2 <- numeric(length(idx))
    dd2[-na.idx] <- dd
    for (i in 1:3) { ## downstream within 3k
        j <- which(annotation == "Intergenic" & abs(dd2) <= i*1000 & dd2 != 0)
        if (length(j) > 0) {
            if (i == 1) {
                lbs <- "Downstream (<1kb)"
            } else {
                lbs <- paste("Downstream (", i-1, "-", i, "kb)", sep="")
            }
            annotation[j] <- lbs
        }
    }
    annotation[which(annotation == "Intergenic")] = "Distal Intergenic"

    downstreamIndex <- dd2 > 0 & dd2 < 3000 ## downstream 3k
    detailGenomicAnnotation[downstreamIndex, "downstream"] <- TRUE
    detailGenomicAnnotation[which(annotation == "Distal Intergenic"), "distal_intergenic"] <- TRUE
    return(list(annotation=annotation, detailGenomicAnnotation=detailGenomicAnnotation))
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
    } else if (type == "fiveUTR") {
        anno <- "5' UTR"
    } else if (type == "threeUTR") {
        anno <- "3' UTR"
    } else {
        anno <- type
    }
    res <- list(queryIndex=queryIndex, annotation=anno)
    return(res)
}
