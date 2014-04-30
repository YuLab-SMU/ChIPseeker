##' @importFrom AnnotationDbi get
.ChIPseekerEnv <- function(TranscriptDb) {
    if (!exists("ChIPseekerEnv", envir=.GlobalEnv)) {
        assign("ChIPseekerEnv", new.env(), .GlobalEnv)
    }
    
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    if (!exists("TXDB", envir=ChIPseekerEnv, inherits=FALSE)) {
        ## first run
        assign("TXDB", TranscriptDb, envir=ChIPseekerEnv)
    } else {
        TXDB <- get("TXDB", envir=ChIPseekerEnv)
        m1 <- unlist(metadata(TXDB))
        m2 <- unlist(metadata(TranscriptDb))
        m1 <- m1[!is.na(m1)]
        m2 <- m2[!is.na(m2)]

        if ( any(m1 != m2) ) {
            rm(ChIPseekerEnv)
            assign("ChIPseekerEnv", new.env(), .GlobalEnv)
            ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
            assign("TXDB", TranscriptDb, envir=ChIPseekerEnv)
        }
    }
    
}

getCols <- function(n) {
    colorRampPalette(brewer.pal(12, "Set3"))(n)  
}


getTagCount <- function(tagMatrix, xlim) {
    ss <- colSums(tagMatrix)
    ss <- ss/sum(ss)
    ## plot(1:length(ss), ss, type="l", xlab=xlab, ylab=ylab)
    pos <- value <- NULL
    dd <- data.frame(pos=c(xlim[1]:xlim[2]), value=ss)
    return(dd)
}

updateGenomicAnnotation <- function(peaks, genomicRegion, type, annotation) {
    hits <- getGenomicAnnotation.internal(peaks, genomicRegion, type)
    if (length(hits) == 2) {
        annotation[hits$queryIndex] <- hits$annotation
    }
    return(annotation)
 }         
 

##' @importFrom IRanges elementLengths
##' @importFrom IRanges findOverlaps
##' @importFrom IRanges queryHits
##' @importFrom IRanges subjectHits
##' @importMethodsFrom BiocGenerics unlist
getGenomicAnnotation.internal <- function(peaks, genomicRegion, type){
    GRegion <- unlist(genomicRegion)
    GRegionLen <- elementLengths(genomicRegion)
    names(GRegionLen) <- names(genomicRegion)
    GRegion$gene_id <- rep(names(genomicRegion), times=GRegionLen)

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
        anno <- paste(type, " (", geneID, " intron ", hits$intron_rank,
                      " of ", GRegionLen[geneID], ")", sep="")
    } else if (type == "Exon") {
        anno <- paste(type, " (", geneID, " exon ", hits$exon_rank,
                      " of ", GRegionLen[geneID], ")", sep="")
    } else {
        anno <- type
    }
    res <- list(queryIndex=queryIndex, annotation=anno)
    return(res)
}

getFirstHitIndex <- function(x) {
    sapply(unique(x), function(i) which(x == i)[1])
}

##' calculate the overlap matrix, which is useful for vennplot
##'
##' 
##' @title overlap 
##' @param Sets a list of objects
##' @return data.frame
##' @importFrom gtools permutations
##' @export
##' @author G Yu
overlap <- function(Sets) {
    ## this function is very generic.
    ## it call the getIntersectLength function to calculate
    ## the number of the intersection.
    ## if it fail, take a look at the object type were supported by getIntersectLength function.
    
    nn <- names(Sets)
    w <- t(apply(permutations(2,length(Sets),0:1, repeats.allowed=TRUE), 1 , rev))
    rs <- rowSums(w)
    wd <- as.data.frame(w)
    wd$n <- NA
    for (i in length(nn):0) {
        idx <- which(rs == i)
        if (i == length(nn)) {
            len <- getIntersectLength(Sets, as.logical(w[idx,]))
            wd$n[idx] <- len
        } else if (i == 0) {
            wd$n[idx] <- 0
        } else {
            for (ii in idx) {
                ##print(ii)
                len <- getIntersectLength(Sets, as.logical(w[ii,]))
                ww = w[ii,]
                jj <- which(ww == 0)
                pp <- permutations(2, length(jj), 0:1, repeats.allowed=TRUE)
                
                for (aa in 2:nrow(pp)) {
                    ## 1st row is all 0, abondoned
                    xx <- jj[as.logical(pp[aa,])]
                    ww[xx] =ww[xx] +1
                    bb <-  t(apply(w, 1, function(i) i == ww))
                    wd$n[rowSums(bb) == length(ww) ]
                         ww <- w[ii,]
                    len <- len - wd$n[rowSums(bb) == length(ww) ]
                    ww <- w[ii,]
                }
                wd$n[ii] <- len
            }
        }
    }
    colnames(wd) = c(names(Sets), "Weight")
    return(wd)
}


getIntersectLength <- function(Sets, idx) {
    ## only use intersect and length methods in this function
    ## works fine with GRanges object
    ## and easy to extend to other objects.
    ss= Sets[idx]
    ol <- ss[[1]]
    
    if (sum(idx) == 1) {
        return(length(ol))
    }
    
    for (j in 2:length(ss)) {
        ol <-  intersect(ol, ss[[j]])
    }
    return(length(ol))
}

shuffle <- function(peak, chrLens) {
    peak.gr <- loadPeak(peak)

    nn <- as.vector(seqnames(peak.gr))
    ii <- order(nn)
    w <- width(peak.gr)
    nnt <- table(nn)
    jj <- order(names(nnt))
    nnt <- nnt[jj]
    chrLens <- chrLens[jj]
    ss <- unlist(sapply(1:length(nnt), function(i) sample(chrLens[i],nnt[i])))

    res <- GRanges(seqnames=nn[ii], ranges=IRanges(ss, ss+w[ii]), strand="*")
    return(res)   
}

loadPeak <- function(peak, verbose=FALSE) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        if (verbose)
            cat(">> loading peak file...\t\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
        peak.gr <- readPeakFile(peak, as="GRanges")
    } else {
        stop("peak should be GRanges object or a peak file...")
    }
    return(peak.gr)
}

##' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
loadTxDb <- function(TranscriptDb) {
    if ( is.null(TranscriptDb) ) {
        TranscriptDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    return(TranscriptDb)
}

##' @importFrom AnnotationDbi get
##' @importFrom GenomicFeatures genes
##' @importFrom GenomicFeatures transcriptsBy
getGene <- function(TranscriptDb, by="gene") {
    .ChIPseekerEnv(TranscriptDb)
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

    if (by == "gene") {
        if ( exists("features", envir=ChIPseekerEnv, inherits=FALSE) ) {
            features <- get("features", envir=ChIPseekerEnv)
        } else {
            features <- genes(TranscriptDb)
            assign("features", features, envir=ChIPseekerEnv)
        }
    } else {
        if ( exists("Transcripts", envir=ChIPseekerEnv, inherits=FALSE) ) {
            features <- get("Transcripts", envir=ChIPseekerEnv)
        } else {
            features <- transcriptsBy(TranscriptDb)
            features <- unlist(features)
            assign("Transcripts", features, envir=ChIPseekerEnv)
        }
    }
    
    return(features)
}


##' get filenames of sample files
##'
##' 
##' @title getSampleFiles 
##' @return list of file names
##' @export
##' @author G Yu
getSampleFiles <- function() {
    dir <- system.file("extdata", "GSE40740", package="ChIPseeker")
    files <- list.files(dir)
    protein <- sub("GSM\\d+_", "", files)
    protein <- sub("_.+", "", protein)
    res <- paste(dir, files, sep="/")
    res <- as.list(res)
    names(res) <- protein
    return(res)
}
## @importFrom RCurl getURL
## getDirListing <- function (url) {
##     ## from GEOquery
##     print(url)
##     a <- getURL(url)
##     b <- textConnection(a)
##     d <- read.table(b, header = FALSE)
##     close(b)
##     return(d)
## }
