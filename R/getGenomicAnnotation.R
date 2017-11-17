updateGenomicAnnotation <- function(peaks, genomicRegion, type, anno, sameStrand=FALSE) {
    hits <- getGenomicAnnotation.internal(peaks, genomicRegion, type, sameStrand=sameStrand)
    if (length(hits) > 1) {
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
##' @param genomicAnnotationPriority genomic Annotation Priority
##' @param sameStrand whether annotate gene in same strand
##' @importFrom GenomicFeatures threeUTRsByTranscript
##' @importFrom GenomicFeatures fiveUTRsByTranscript
##' @return character vector
##' @author G Yu
getGenomicAnnotation <- function(peaks,
                                 distance,
                                 tssRegion=c(-3000, 3000),
                                 TxDb,
                                 level,
                                 genomicAnnotationPriority,
                                 sameStrand = FALSE
                                 ) {

    ##
    ## since some annotation overlap,
    ## a priority is assign based on *genomicAnnotationPriority*
    ## use the following priority by default:
    ##
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

    anno <- list(annotation=annotation,
                 detailGenomicAnnotation=detailGenomicAnnotation)

    genomicAnnotationPriority <- rev(genomicAnnotationPriority)
    for (AP in genomicAnnotationPriority) {
        if (AP == "Intergenic") {
            ## Intergenic
            annotation[is.na(annotation)] <- "Intergenic"
            anno[["annotation"]] <- annotation
        } else if (AP == "Intron") {
            ## Introns
            intronList <- get_intronList(ChIPseekerEnv)
            anno <- updateGenomicAnnotation(peaks, intronList, "Intron", anno, sameStrand=sameStrand)
        } else if (AP == "Exon") {
            ## Exons
            exonList <- get_exonList(ChIPseekerEnv)
            anno <- updateGenomicAnnotation(peaks, exonList, "Exon", anno, sameStrand=sameStrand)
        } else if (AP == "3UTR") {
            ## 3' UTR Exons
            if ( exists("threeUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
                threeUTRList <- get("threeUTRList", envir=ChIPseekerEnv)
            } else {
                threeUTRList <- threeUTRsByTranscript(TxDb)
                assign("threeUTRList", threeUTRList, envir=ChIPseekerEnv)
            }
            anno <- updateGenomicAnnotation(peaks, threeUTRList, "threeUTR", anno, sameStrand=sameStrand)
        } else if (AP == "5UTR") {
            ## 5' UTR Exons
            if ( exists("fiveUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
                fiveUTRList <- get("fiveUTRList", envir=ChIPseekerEnv)
            } else {
                fiveUTRList <- fiveUTRsByTranscript(TxDb)
                assign("fiveUTRList", fiveUTRList, envir=ChIPseekerEnv)
            }
            anno <- updateGenomicAnnotation(peaks, fiveUTRList, "fiveUTR", anno, sameStrand=sameStrand)
        }

        annotation <- anno[["annotation"]]
        detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]

        if (AP == "Promoter") {
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
        }

    }


    genicIndex <- which(apply(detailGenomicAnnotation[, c("Exon", "Intron")], 1, any))
    detailGenomicAnnotation[-genicIndex, "Intergenic"] <- TRUE
    detailGenomicAnnotation[genicIndex, "genic"] <- TRUE

    ## intergenicIndex <- anno[["annotation"]] == "Intergenic"
    ## anno[["detailGenomicAnnotation"]][intergenicIndex, "Intergenic"] <- TRUE
    ## anno[["detailGenomicAnnotation"]][!intergenicIndex, "genic"] <- TRUE


    features <- getGene(TxDb, by=level)

    ## nearest from gene end
    if (sameStrand) {
        idx <- precede(peaks, features)
    } else {
        idx <- precede(peaks, unstrand(features))
    }
    na.idx <- which(is.na(idx))
    if (length(na.idx)) {
        idx <- idx[-na.idx]
        peaks <- peaks[-na.idx]
    }
    peF <- features[idx]
    dd <- ifelse(strand(peF) == "+",
                 start(peaks) - end(peF),
                 end(peaks) - start(peF))
    if (length(na.idx)) {
        dd2 <- numeric(length(idx) + length(na.idx))
        dd2[-na.idx] <- dd
    } else {
        dd2 <- dd
    }

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

    dsd <- getOption("ChIPseeker.downstreamDistance")
    if (is.null(dsd))
        dsd <- 3000 ## downstream 3k by default

    downstreamIndex <- dd2 > 0 & dd2 < dsd
    detailGenomicAnnotation[downstreamIndex, "downstream"] <- TRUE
    detailGenomicAnnotation[which(annotation == "Distal Intergenic"), "distal_intergenic"] <- TRUE
    return(list(annotation=annotation, detailGenomicAnnotation=detailGenomicAnnotation))
}


##' @import BiocGenerics S4Vectors IRanges
getGenomicAnnotation.internal <- function(peaks, genomicRegion, type, sameStrand=FALSE){
    GRegion <- unlist(genomicRegion)
    GRegionLen <- elementNROWS(genomicRegion)

    names(GRegionLen) <- names(genomicRegion)
    GRegion$gene_id <- rep(names(genomicRegion), times=GRegionLen)


    if (type == "Intron") {
        gr2 <- GRegion[!duplicated(GRegion$gene_id)]
        strd <- as.character(strand(gr2))
        len <- GRegionLen[GRegionLen != 0]

        GRegion$intron_rank <- lapply(seq_along(strd), function(i) {
            rank <- seq(1, len[i])
            if (strd[i] == '-')
                rank <- rev(rank)
            return(rank)
        }) %>% unlist
    }

    if (type == "Intron" || type =="Exon") {
        nn <- TXID2EG(names(genomicRegion))
        names(GRegionLen) <- nn
        GRegion$gene_id <- rep(nn, times=GRegionLen)
    }

    ## find overlap
    if (sameStrand) {
        GRegionHit <- findOverlaps(peaks, GRegion)
    } else {
        GRegionHit <- findOverlaps(peaks, unstrand(GRegion))
    }

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
    res <- list(queryIndex=queryIndex, annotation=anno, gene=geneID)
    return(res)
}
