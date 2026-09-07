updateGenomicAnnotation <- function(peaks, genomicRegion, type,
                                    anno, sameStrand=FALSE) {
    hits <- getGenomicAnnotation.internal(peaks, genomicRegion,
                                          type, sameStrand=sameStrand)
    if (length(hits) > 1) {
        hitIndex <- hits$queryIndex
        anno[["annotation"]][hitIndex] <- hits$annotation
        anno[["detailGenomicAnnotation"]][hitIndex, type] <- TRUE
    }
    return(anno)
}


##' Get genomic annotation of peaks
##'
##' This function assigns genomic feature categories to peaks based on their
##' location relative to gene structures (promoters, UTRs, exons, introns, etc.).
##'
##' @description
##' The function categorizes each peak into one of several genomic feature types
##' by checking for overlaps with different genomic regions. When a peak overlaps
##' multiple features, a priority system is used to assign the most specific
##' annotation. The function also provides detailed annotation information
##' indicating which features each peak overlaps.
##'
##' @details
##' The annotation process follows these steps:
##' \enumerate{
##'   \item Initializes annotation vectors and a detailed annotation data frame
##'   \item Processes genomic regions in reverse priority order (from lowest to
##'         highest priority) to ensure higher priority annotations overwrite
##'         lower priority ones
##'   \item Checks for overlaps with: Introns, Exons, 3' UTR, 5' UTR, Promoter
##'         (based on distance to TSS), and Intergenic regions
##'   \item For promoter regions, creates distance-based subcategories (e.g.,
##'         "Promoter (<=1kb)", "Promoter (1-2kb)") when the TSS region is >= 2kb
##'   \item Identifies downstream regions (within a configurable distance from
##'         gene end) and categorizes remaining intergenic peaks as "Distal Intergenic"
##'   \item Marks peaks as "genic" if they overlap exons or introns, otherwise
##'         marks as "Intergenic"
##' }
##'
##' The priority order (from highest to lowest) is typically: Promoter > 5' UTR >
##' 3' UTR > Exon > Intron > Downstream > Intergenic. This can be customized via
##' \code{genomicAnnotationPriority}.
##'
##' @param peaks GRanges object containing genomic ranges of peaks to be annotated
##' @param distance numeric vector of distances from each peak to the TSS of the
##'   nearest gene. Positive values indicate downstream, negative values indicate
##'   upstream. Used to identify promoter regions
##' @param tssRegion numeric vector of length 2 specifying the TSS region for
##'   promoter annotation. Default is c(-3000, 3000), meaning 3kb upstream and
##'   3kb downstream of TSS. Peaks within this region are annotated as "Promoter"
##' @param TxDb TxDb or EnsDb annotation object containing gene/transcript
##'   structure information
##' @param level character, one of "gene" or "transcript". Determines whether to
##'   annotate at gene level or transcript level
##' @param genomicAnnotationPriority character vector specifying the priority
##'   order of genomic annotations. Must be a permutation of c("Promoter", "5UTR",
##'   "3UTR", "Exon", "Intron", "Downstream", "Intergenic"). Higher priority
##'   annotations (appearing earlier in the vector) will overwrite lower priority
##'   ones when a peak overlaps multiple features
##' @param sameStrand logical, whether to only consider overlaps with features on
##'   the same strand as the peak. If FALSE, strand is ignored. Default is FALSE
##' @return A list with two components:
##'   \itemize{
##'     \item \code{annotation}: Character vector of length equal to input peaks,
##'       containing the assigned genomic annotation for each peak. Possible
##'       values include: "Promoter" (with distance subcategories if applicable),
##'       "5' UTR", "3' UTR", "Exon (geneID, exon N of M)", "Intron (geneID,
##'       intron N of M)", "Downstream (distance ranges)", "Distal Intergenic"
##'     \item \code{detailGenomicAnnotation}: Data frame with logical columns
##'       indicating which features each peak overlaps: \code{genic},
##'       \code{Intergenic}, \code{Promoter}, \code{fiveUTR}, \code{threeUTR},
##'       \code{Exon}, \code{Intron}, \code{downstream}, \code{distal_intergenic}
##'   }
##' @importFrom GenomicFeatures threeUTRsByTranscript
##' @importFrom GenomicFeatures fiveUTRsByTranscript
##' @importFrom yulab.utils get_cache_element
##' @importFrom yulab.utils update_cache_item
##' @author G Yu
##' @noRd
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



    .ChIPseekerEnv(TxDb, item = ChIPseekerCache)
    # ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

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
        if (AP == "Intron") {
            ## Introns
            # intronList <- get_intronList(ChIPseekerEnv)
            intronList <- get_intronList(item = ChIPseekerCache)
            anno <- updateGenomicAnnotation(peaks, intronList, "Intron", anno, sameStrand=sameStrand)
        } else if (AP == "Exon") {
            ## Exons
            # exonList <- get_exonList(ChIPseekerEnv)
            exonList <- get_exonList(item = ChIPseekerCache)
            anno <- updateGenomicAnnotation(peaks, exonList, "Exon", anno, sameStrand=sameStrand)
        } else if (AP == "3UTR") {
            ## 3' UTR Exons
            threeUTRList <- get_cache_element(item = ChIPseekerCache, elements = "threeUTRList")

            if(is.null(threeUTRList)){
                threeUTRList <- threeUTRsByTranscript(TxDb)
                update_cache_item(item = ChIPseekerCache, list("threeUTRList" = threeUTRList))
            }

            # if ( exists("threeUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
            #     threeUTRList <- get("threeUTRList", envir=ChIPseekerEnv)
            # } else {
            #     threeUTRList <- threeUTRsByTranscript(TxDb)
            #     assign("threeUTRList", threeUTRList, envir=ChIPseekerEnv)
            # }
            anno <- updateGenomicAnnotation(peaks, threeUTRList, "threeUTR", anno, sameStrand=sameStrand)
        } else if (AP == "5UTR") {
            ## 5' UTR Exons
            fiveUTRList <- get_cache_element(item = ChIPseekerCache, elements = "fiveUTRList")

            if(is.null(fiveUTRList)){
                fiveUTRList <- fiveUTRsByTranscript(TxDb)
                update_cache_item(item = ChIPseekerCache, list("fiveUTRList" = fiveUTRList))
            }

            # if ( exists("fiveUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
            #     fiveUTRList <- get("fiveUTRList", envir=ChIPseekerEnv)
            # } else {
            #     fiveUTRList <- fiveUTRsByTranscript(TxDb)
            #     assign("fiveUTRList", fiveUTRList, envir=ChIPseekerEnv)
            # }
            anno <- updateGenomicAnnotation(peaks, fiveUTRList, "fiveUTR", anno, sameStrand=sameStrand)
        } else if (AP == "Promoter") {
            annotation <- anno[["annotation"]]
            ## detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]

            ## TSS
            tssIndex <- distance >= tssRegion[1] & distance <= tssRegion[2]
            annotation[tssIndex] <- "Promoter"
            anno$detailGenomicAnnotation[tssIndex, "Promoter"] <- TRUE

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
            anno[["annotation"]] <- annotation
        } else {
            ## Intergenic
            annotation[is.na(annotation)] <- "Intergenic"
            anno[["annotation"]] <- annotation
        }
    }

    annotation <- anno[["annotation"]]
    detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]
    genicIndex <- which(apply(detailGenomicAnnotation[, c("Exon", "Intron")], 1, any))
    detailGenomicAnnotation[-genicIndex, "Intergenic"] <- TRUE
    detailGenomicAnnotation[genicIndex, "genic"] <- TRUE

    ## intergenicIndex <- anno[["annotation"]] == "Intergenic"
    ## anno[["detailGenomicAnnotation"]][intergenicIndex, "Intergenic"] <- TRUE
    ## anno[["detailGenomicAnnotation"]][!intergenicIndex, "genic"] <- TRUE


    features <- getGene(TxDb, by=level)

    ## nearest from gene end
    if (sameStrand) {
        idx <- follow(peaks, features)
    } else {
        idx <- follow(peaks, unstrand(features))
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

    dsd <- getOption("ChIPseeker.downstreamDistance")
    if (is.null(dsd))
	    dsd <- 3000 ## downstream 3k by default

    ## downstream within dsd
    if(dsd/1000<=1){
        j <- which(annotation == "Intergenic" & abs(dd2) <= dsd & dd2 != 0)
        if(length(j)>0){
            lbs <- paste("Downstream (<=", dsd, "bp)", sep="")
            annotation[j] <- lbs
        }
    }else{

        ## downstream within 0-dsd/1000 kb
        for(i in 1:(dsd/1000)){
            j <- which(annotation == "Intergenic" & abs(dd2) <= i*1000 & dd2 != 0)
            if (length(j) > 0){
                if (i == 1){
                    lbs <- "Downstream (<1kb)"
                }else{
                    lbs <- paste("Downstream (", i-1, "-", i, "kb)", sep="")
                }
		annotation[j] <- lbs
            }
        }

        ## downstream (dsd/1000) kb - dsd bp
        z <- which(annotation == "Intergenic" & abs(dd2) <= dsd & dd2 != 0)
        if(length(z)>0){
            lbs <- paste("Downstream (",dsd/1000,"kb-", dsd, "bp)", sep="")
            annotation[z] <- lbs
        }
    }
    annotation[which(annotation == "Intergenic")] = "Distal Intergenic"

    downstreamIndex <- dd2 > 0 & dd2 < dsd
    detailGenomicAnnotation[downstreamIndex, "downstream"] <- TRUE
    detailGenomicAnnotation[which(annotation == "Distal Intergenic"), "distal_intergenic"] <- TRUE
    return(list(annotation=annotation, detailGenomicAnnotation=detailGenomicAnnotation))
}


##' @import BiocGenerics S4Vectors IRanges
##' @noRd
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
