##' @importFrom AnnotationDbi get
.ChIPseekerEnv <- function(TxDb) {
    if (!exists("ChIPseekerEnv", envir=.GlobalEnv)) {
        assign("ChIPseekerEnv", new.env(), .GlobalEnv)
    }
    
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    if (!exists("TXDB", envir=ChIPseekerEnv, inherits=FALSE)) {
        ## first run
        assign("TXDB", TxDb, envir=ChIPseekerEnv)
    } else {
        TXDB <- get("TXDB", envir=ChIPseekerEnv)
        m1 <- unlist(metadata(TXDB))
        m2 <- unlist(metadata(TxDb))
        m1 <- m1[!is.na(m1)]
        m2 <- m2[!is.na(m2)]

        if ( length(m1) != length(m2) || any(m1 != m2) ) {
            rm(ChIPseekerEnv)
            assign("ChIPseekerEnv", new.env(), .GlobalEnv)
            ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
            assign("TXDB", TxDb, envir=ChIPseekerEnv)
        }
    }
    
}

getCols <- function(n) {
    col <- c("#8dd3c7", "#ffffb3", "#bebada",
             "#fb8072", "#80b1d3", "#fdb462",
             "#b3de69", "#fccde5", "#d9d9d9",
             "#bc80bd", "#ccebc5", "#ffed6f")
             
    col2 <- c("#1f78b4", "#ffff33", "#c2a5cf",
             "#ff7f00", "#810f7c", "#a6cee3",
             "#006d2c", "#4d4d4d", "#8c510a",
             "#d73027", "#78c679", "#7f0000",
             "#41b6c4", "#e7298a", "#54278f")
    
    col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
              "#33a02c", "#fb9a99", "#e31a1c",
              "#fdbf6f", "#ff7f00", "#cab2d6",
              "#6a3d9a", "#ffff99", "#b15928")
    
    ## colorRampPalette(brewer.pal(12, "Set3"))(n)
    colorRampPalette(col3)(n)  
}


getTagCount <- function(tagMatrix, xlim, conf) {
    ss <- colSums(tagMatrix)
    ss <- ss/sum(ss)
    ## plot(1:length(ss), ss, type="l", xlab=xlab, ylab=ylab)
    pos <- value <- NULL
    dd <- data.frame(pos=c(xlim[1]:xlim[2]), value=ss)
    return(dd)
}


##
## estimate CI using bootstraping
##
getSgn <- function(data, idx){
    d <- data[idx, ]
    ss <- colSums(d)
    ss <- ss / sum(ss)
    return(ss)
}
parseBootCiPerc <- function(bootCiPerc){
    bootCiPerc <- bootCiPerc$percent
    tmp <- length(bootCiPerc)
    ciLo <- bootCiPerc[tmp - 1]
    ciUp <- bootCiPerc[tmp]
    return(c(ciLo, ciUp))
}
##' @importFrom boot boot
##' @importFrom boot boot.ci
getTagCountCI <- function(tagMatrix, xlim, conf = 0.95){
    RESAMPLE_TIME <- 500
    trackLen <- ncol(tagMatrix)
    tagMxBoot <- boot(data = tagMatrix, statistic = getSgn, R = RESAMPLE_TIME)
    tagMxBootCi <- sapply(seq_len(trackLen), function(i) {
                        bootCiToken <- boot.ci(tagMxBoot, type = "perc", index = i) 
                        ## parse boot.ci results
                        return(parseBootCiPerc(bootCiToken))
                        }
                    )
    row.names(tagMxBootCi) <- c("Lower", "Upper")
    # pos <- value <- NULL
    # dd <- data.frame(pos = c(xlim[1]:xlim[2]), 
    #                 lw = tagMxBootCi[1, ], up = tagMxBootCi[2, ])
    # return(dd)
    return(tagMxBootCi)
}


TXID2EG <- function(txid, geneIdOnly=FALSE) {
    txid <- as.character(txid)
    if (geneIdOnly == TRUE) {
        res <- TXID2EGID(txid)
    } else {
        res <- TXID2TXEG(txid)
    }
    return(res)
}

##' @importFrom GenomicFeatures transcripts
TXID2TXEG <- function(txid) {
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    
    if (exists("txid2geneid", envir=ChIPseekerEnv, inherits=FALSE)) {
        txid2geneid <- get("txid2geneid", envir=ChIPseekerEnv)
    } else {
        txdb <- get("TXDB", envir=ChIPseekerEnv)
        txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
        idx <- which(sapply(txidinfo$gene_id, length) == 0)
        txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
        txid2geneid <- paste(mcols(txidinfo)[["tx_name"]],
                             mcols(txidinfo)[["gene_id"]],
                             sep="/")
        txid2geneid <- sub("/NA", "", txid2geneid)
        
        names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
        assign("txid2geneid", txid2geneid, envir=ChIPseekerEnv)
    }
    return(as.character(txid2geneid[txid]))
}

TXID2EGID <- function(txid) {
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
    
    if (exists("txid2eg", envir=ChIPseekerEnv, inherits=FALSE)) {
        txid2geneid <- get("txid2eg", envir=ChIPseekerEnv)
    } else {
        txdb <- get("TXDB", envir=ChIPseekerEnv)
        txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
        idx <- which(sapply(txidinfo$gene_id, length) == 0)
        txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
        txid2geneid <- as.character(mcols(txidinfo)[["gene_id"]])
                
        names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
        assign("txid2eg", txid2geneid, envir=ChIPseekerEnv)
    }
    return(as.character(txid2geneid[txid]))
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
loadTxDb <- function(TxDb) {
    if ( is.null(TxDb) ) {
        TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    return(TxDb)
}

##' @importFrom AnnotationDbi get
##' @importFrom GenomicFeatures genes
##' @importFrom GenomicFeatures transcriptsBy
getGene <- function(TxDb, by="gene") {
    .ChIPseekerEnv(TxDb)
    ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)

    by <- match.arg(by, c("gene", "transcript"))
    
    if (by == "gene") {
        if ( exists("features", envir=ChIPseekerEnv, inherits=FALSE) ) {
            features <- get("features", envir=ChIPseekerEnv)
        } else {
            features <- genes(TxDb)
            assign("features", features, envir=ChIPseekerEnv)
        }
    } else {
        if ( exists("Transcripts", envir=ChIPseekerEnv, inherits=FALSE) ) {
            features <- get("Transcripts", envir=ChIPseekerEnv)
        } else {
            features <- transcriptsBy(TxDb)
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
    dir <- system.file("extdata", "GEO_sample_data", package="ChIPseeker")
    files <- list.files(dir)
    ## protein <- sub("GSM\\d+_", "", files)
    ## protein <- sub("_.+", "", protein)
    protein <- gsub(pattern='GSM\\d+_(\\w+_\\w+)_.*', replacement='\\1',files)
    protein <- sub("_Chip.+", "", protein)
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


is.dir <- function(dir) {
    if (file.exists(dir) == FALSE)
        return(FALSE)
    return(file.info(dir)$isdir)
}


parse_targetPeak_Param <- function(targetPeak) {
    if (length(targetPeak) == 1) {
        if (is.dir(targetPeak)) {
            files <- list.files(path=targetPeak)
            idx <- unlist(sapply(c("bed", "bedGraph", "Peak"), grep, x=files))
            idx <- sort(unique(idx))
            files <- files[idx]
            targetPeak <- sub("/$", "", targetPeak)
            res <- paste(targetPeak, files, sep="/")
        } else {
            if (!file.exists(targetPeak)) {
                stop("bed file is not exists...")
            } else {
                res <- targetPeak
            }
        }
    } else {
        if (is.dir(targetPeak[1])) {
            stop("targetPeak should be a vector of bed file names or a folder containing bed files...")
        } else {
            res <- targetPeak[file.exists(targetPeak)]
            if (length(res) == 0) {
                stop("targetPeak file not exists...")
            }
        }
    }
    return(res)
}
