##' Annotate peaks
##'
##' 
##' @title annotatePeak
##' @param peak peak file or GRanges object
##' @param tssRegion Region Range of TSS 
##' @param as one of "data.frame", "GRanges" and "txt"
##' @param TranscriptDb TranscriptDb object
##' @param assignGenomicAnnotation logical, assign peak genomic annotation or not
##' @param annoDb annotation package
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
##' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
##' @importFrom GenomicFeatures genes
##' @importFrom AnnotationDbi get
##' @importMethodsFrom BiocGenerics as.data.frame
##' @examples
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 100), as="GRanges", TranscriptDb=txdb)
##' head(peakAnno)
##' @seealso \code{\link{plotAnnoBar}} \code{\link{plotAnnoPie}} \code{\link{plotDistToTSS}}
##' @export
##' @author G Yu
annotatePeak <- function(peak,
                         tssRegion=c(-3000, 3000),
                         as="GRanges",
                         TranscriptDb=NULL,
                         assignGenomicAnnotation=TRUE,
                         annoDb=NULL,
                         verbose=TRUE) {
    
    output <- match.arg(as, c("data.frame", "GRanges", "txt"))
    
    if ( is(peak, "GRanges") ){
        ## this test will be TRUE
        ## when peak is an instance of class/subclass of "GRanges"
        input <- "gr"
        peak.gr <- peak
    } else {
        if (file.exists(peak)) {
            if (verbose)
                cat(">> loading peak file...\t\t\t\t",
                    format(Sys.time(), "%Y-%m-%d %X"), "\n")
            ## peak.df <- peak2DF(peak)
            ## peak.gr <- peakDF2GRanges(peak.df)
            input <- "file"
            peak.gr <- readPeakFile(peak, as="GRanges")
        } else {
            stop("peak should be a GRanges object or a peak file...")
        }
    }
    
    if (verbose)
        cat(">> preparing features information...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    if ( is.null(TranscriptDb) ) {
        TranscriptDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }

    if (!exists("ChIPseekerEnv"))
        .ChIPseekerEnv(TranscriptDb)
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    
    if ( exists("features", envir=ChIPseekerEnv, inherits=FALSE) ) {
        features <- get("features", envir=ChIPseekerEnv)
    } else {
        features <- genes(TranscriptDb)
        assign("features", features, envir=ChIPseekerEnv)
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

    if (verbose)
        cat(">> assigning genomic annotation...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    ## annotation
    if (assignGenomicAnnotation == TRUE) {
        annotation <- getGenomicAnnotation(peak.gr, distance, tssRegion, TranscriptDb)
    } else {
        annotation <- NULL
    }
    ## duplicated names since more than 1 peak may annotated by only 1 gene
    names(nearestFeatures) <- NULL
    nearestFeatures.df <- as.data.frame(nearestFeatures)
    colnames(nearestFeatures.df) <- c("geneChr", "geneStart", "geneEnd",
                                      "geneLength", "geneStrand", "geneId")

    ## append annotation to peak.gr
    if (!is.null(annotation))
        elementMetadata(peak.gr)[["annotation"]] <- annotation

    for(cn in colnames(nearestFeatures.df)) {
        elementMetadata(peak.gr)[[cn]] <- unlist(nearestFeatures.df[, cn])
    }

    elementMetadata(peak.gr)[["distanceToTSS"]] <- distance
    ## res <- cbind(peak.df,
    ##              annotation=annotation,
    ##              nearestFeatures.df,
    ##              distanceToTSS=distance)
    if (!is.null(annoDb)) {
        if (verbose)
            cat(">> adding gene annotation...\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
        geneAnno <- addGeneAnno(annoDb, peak.gr$geneId)
        for(cn in colnames(geneAnno)[-1])
            elementMetadata(peak.gr)[[cn]] <- geneAnno[, cn]
        
        ## out <- cbind(res, geneAnno[,-1])
    }

    if (output == "txt") {
        if (input == "gr") {
            outfile <- "peak_anno.txt"
        } else {
            outfile <- sub("\\.\\w+", "_anno.txt", peak)
        }
        write.table(as.data.frame(peak.gr), file=outfile,
                    sep="\t", row.names=FALSE, quote=FALSE)
    }
    
    if(verbose)
        cat(">> done...\t\t\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")

    if (output == "df") {
        return(as.data.frame(peak.gr))
    }
    return(peak.gr)
}

##' add gene annotation, symbol, gene name etc.
##'
##' 
##' @title addGeneAnno 
##' @param annoDb annotation package
##' @param geneID query geneID
##' @return data.frame
##' @importFrom AnnotationDbi select
##' @author G Yu
addGeneAnno <- function(annoDb, geneID){
    kk <- unlist(geneID)
    require(annoDb, character.only = TRUE)
    annoDb <- eval(parse(text=annoDb))
    ann <- suppressWarnings(select(annoDb,
                                   keys=kk,
                                   keytype="ENTREZID",
                                   columns=c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")))
    idx <- getFirstHitIndex(ann$ENTREZID)
    ann <- ann[idx,]

    idx <- unlist(sapply(kk, function(x) which(x==ann$ENTREZID)))
    ann <- ann[idx,]
    return(ann)
}



## ## the function works but too slow
## getNearestIndexAndDistance <- function(peak, features) {
##     ## peak only conatin one peak record, in GRanges object
##     ## feature is the annotation in GRanges object
       
##     ## only keep start position based on strand
##     start(features) <- end(features) <- ifelse(strand(features) == "+", start(features), end(features))
##     ## feature selected
##     ## match seqnames
##     fs <- features[seqnames(features) == as.character(seqnames(peak))]

##     ## peak start
##     ps <- rep(peak, length(fs))
##     start(ps) <- end(ps) <- ifelse(strand(fs) == "+", start(ps), end(ps))

##     dd <- start(ps) - start(fs)
##     ii <- which.min(abs(dd))

##     idx <- which(fs[ii] == features)
##     d <- dd[ii]
##     result <- c(idx,d)
##     names(result) <- c("index", "dist")
##     return(result)
## }

##' get index of features that closest to peak and calculate distance
##'
##' 
##' @title getNearestFeatureIndicesAndDistances 
##' @param peaks peak in GRanges 
##' @param features features in GRanges
##' @return data.frame
##' @importFrom IRanges precede
##' @importFrom IRanges follow
##' @importFrom IRanges start
##' @importFrom IRanges end
##' @importFrom BiocGenerics strand
## @importMethodsFrom GenomicRanges strand
##' @author G Yu
getNearestFeatureIndicesAndDistances <- function(peaks, features) {
    ## peaks only conatin all peak records, in GRanges object
    ## feature is the annotation in GRanges object
    
    ## only keep start position based on strand
    start(features) <- end(features) <- ifelse(strand(features) == "+", start(features), end(features))
    
    ## nearest from peak start
    ps.idx <- precede(peaks, features)
    ## nearest from peak end
    pe.idx <- follow(peaks, features)
    
    ## features from nearest peak start
    psF <- features[ps.idx]
    ## feature distances from peak start
    psD <- ifelse(strand(psF) == "+",
                 start(peaks) - start(psF),
                 end(peaks)-end(psF))

    ## features from nearest peak end
    peF <- features[pe.idx]
    ## feature distances from peak end
    peD <- ifelse(strand(peF) == "+",
                  start(peF) - start(peaks),
                  end(peF)-end(peaks))

    pse <- data.frame(ps=psD, pe=peD)
    j <- apply(pse, 1, function(i) which.min(abs(i)))

    ## index
    idx <- ps.idx
    idx[j==2] <- pe.idx[j==2]
    
    ## distance
    dd <- psD
    dd[j==2] <- peD[j==2]
    
    res <- data.frame(index=idx, distance=dd)
    return(res)
}


##' get Genomic Annotation of peaks
##'
##' 
##' @title getGenomicAnnotation
##' @param peaks peaks in GRanges object
##' @param distance distance of peak to TSS
##' @param tssRegion tssRegion, default is -3kb to +100bp
##' @param TranscriptDb TranscriptDb object
##' @importFrom GenomicFeatures intronsByTranscript
##' @importFrom GenomicFeatures threeUTRsByTranscript
##' @importFrom GenomicFeatures fiveUTRsByTranscript
##' @importFrom GenomicFeatures exonsBy
##' @return character vector
##' @author G Yu
getGenomicAnnotation <- function(peaks,
                                 distance,
                                 tssRegion=c(-3000, 100),
                                 TranscriptDb
                                 ) {
    
    ##
    ## since some annotation overlap,
    ## a priority is assign based on the following:
    ## 1. TSS
    ## 2. Exon
    ## 3. 5' UTR
    ## 4. 3' UTR
    ## 5. Intron
    ## 6. Intergenic
    ##

    if (!exists("ChIPseekerEnv"))
        .ChIPseekerEnv(TranscriptDb)
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    
            
    annotation <- rep(NA, length(distance))
    ## Intergenic
    annotation[is.na(annotation)] <- "Intergenic"
    
    ## Introns
    if ( exists("intronList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        intronList <- get("intronList", envir=ChIPseekerEnv)
    } else {
        intronList <- intronsByTranscript(TranscriptDb)
        assign("intronList", intronList, envir=ChIPseekerEnv)
    }
    intronHits <- getGenomicAnnotation.internal(peaks, intronList, "Intron")
    annotation[intronHits$queryIndex] <- intronHits$annotation
      
    ## 3' UTR Exons
    if ( exists("threeUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        threeUTRList <- get("threeUTRList", envir=ChIPseekerEnv)
    } else {
        threeUTRList <- threeUTRsByTranscript(TranscriptDb)
        assign("threeUTRList", threeUTRList, envir=ChIPseekerEnv)
    }
    threeUTRHits <- getGenomicAnnotation.internal(peaks, threeUTRList, "3' UTR")
    annotation[threeUTRHits$queryIndex] <- threeUTRHits$annotation
    
    ## 5' UTR Exons
    if ( exists("fiveUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        fiveUTRList <- get("fiveUTRList", envir=ChIPseekerEnv)
    } else {
        fiveUTRList <- fiveUTRsByTranscript(TranscriptDb)
        assign("fiveUTRList", fiveUTRList, envir=ChIPseekerEnv)
    }
    fiveUTRHits <- getGenomicAnnotation.internal(peaks, fiveUTRList, "5' UTR")
    annotation[fiveUTRHits$queryIndex] <- fiveUTRHits$annotation
    
    ## Exon
    if ( exists("exonList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        exonList <- get("exonList", envir=ChIPseekerEnv)
    } else {
        exonList <- exonsBy(TranscriptDb)
        assign("exonList", exonList, envir=ChIPseekerEnv)
    }
    exonHits <- getGenomicAnnotation.internal(peaks, exonList, "Exon")
    annotation[exonHits$queryIndex] <- exonHits$annotation
    
    ## TSS
    annotation[distance >= tssRegion[1] &
               distance <= tssRegion[2]] <- "Promoter-TSS"
    
    return(annotation)
}
