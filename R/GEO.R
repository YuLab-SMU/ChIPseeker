########################################
##                                    ##
## data last update: Mar 03, 2015     ##
##                                    ##
########################################



##' accessing species statistics collecting from GEO database
##'
##'
##' @title getGEOspecies
##' @return data.frame
##' @author G Yu
##' @export
getGEOspecies <- function() {
    gsminfo <- get_gsminfo()
    species <- gsminfo$organism
    res <- as.data.frame(table(species))
    return(res)
}

##' get genome version statistics collecting from GEO ChIPseq data
##'
##'
##' @title getGEOgenomeVersion
##' @return data.frame
##' @author G Yu
##' @export
getGEOgenomeVersion <- function() {
    gsminfo <- get_gsminfo()
    gv <- gsminfo[, c("organism",
                      "genomeVersion")]
    genomeVersion <- gv$genomeVersion
    res <- as.data.frame(table(genomeVersion))
    gv <- unique(gv)

    res <- merge(gv, res, by.x="genomeVersion", by.y="genomeVersion", all.y=TRUE)
    res <- res[, c("organism", "genomeVersion", "Freq")]
    return(res)
}

##' get subset of GEO information by genome version keyword
##'
##'
##' @title getGEOInfo
##' @param genome genome version
##' @param simplify simplify result or not
##' @return data.frame
##' @author G Yu
##' @export
getGEOInfo <- function(genome, simplify =TRUE) {
    gsminfo <- get_gsminfo()
    genomeVersion <- NULL ## to satisfy codetools
    res <- subset(gsminfo, subset = genomeVersion == genome)
    if (simplify) {
        res <- res[,c("series_id", "gsm", "organism", "title", "supplementary_file", "genomeVersion", "pubmed_id")]
    }
    return(res)
}

##' download all BED files of a particular genome version
##'
##'
##' @title downloadGEObedFiles
##' @param genome genome version
##' @param destDir destination folder
##' @return NULL
##' @author G Yu
##' @export
downloadGEObedFiles <- function(genome, destDir=getwd()) {
    info <- getGEOInfo(genome)
    downloadGEO.internal(info, destDir)
}

##' download BED supplementary files of a list of GSM accession numbers
##'
##'
##' @title downloadGSMbedFiles
##' @param GSM GSM accession numbers
##' @param destDir destination folder
##' @return NULL
##' @author G Yu
##' @export
downloadGSMbedFiles <- function(GSM, destDir=getwd()) {
    gsminfo <- get_gsminfo()
    info <- gsminfo[gsminfo$gsm %in% GSM,]
    downloadGEO.internal(info, destDir)
}

##' @importFrom utils download.file
downloadGEO.internal <- function(info, destDir) {
    fnames <- as.character(info$supplementary_file)
    destfiles <- sub(".*\\/", paste(destDir, "/", sep=""), fnames)
    names(destfiles) <- NULL

    for (i in seq_along(fnames)) {
        if ( ! file.exists(destfiles[i]) )
            tryCatch(download.file(fnames[i],
                          destfile=destfiles[i],
                          mode="wb"),
                     error = function(e) message(fnames[i], ': file not found and skip'))
    }
}

##' @importFrom utils data
## @importFrom GEOmetadb
## @importFrom RSQLite dbConnect
## @importFrom RSQLite dbGetQuery
prepareGSMInfo <- function() {
    pkg <- "GEOmetadb"
    require(pkg, character.only=TRUE)
    getSQLiteFile <- eval(parse(text="getSQLiteFile"))
    ## get the latest version of sql file
    getSQLiteFile()

    GEOmetadbFile="GEOmetadb.sqlite"
    file.info(GEOmetadbFile)

    sqlpkg <- "RSQLite"
    require(sqlpkg, character.only=TRUE)
    dbConnect <- eval(parse(text="dbConnect"))
    dbGetQuery <- eval(parse(text="dbGetQuery"))
    SQLite <- eval(parse(text="SQLite"))

    con <- dbConnect(SQLite(),GEOmetadbFile)
    ## dbListTables(con)


    pkg <- "GEOquery"
    require(pkg, character.only=TRUE)
    getGEO <- eval(parse(text="getGEO"))
    Meta <- eval(parse(text="Meta"))

    ## get all GPL IDs
    ## download soft using gpl = getGEO("GPLXXX")
    ## using Meta(gpl) find the technology match sequencing
    ## get all gsm IDs
    ## parse it

    gpl <- dbGetQuery(con, 'select gpl, technology from gpl')
    gpl <- gpl[gpl[,2] == "high-throughput sequencing",1]
    gpl <- gpl[!is.na(gpl)]

    ## save the processedGSM vector that contain all the GSM that have been processed.
    ## next time when preparing GSMInfo, filter those have been processed before.
    load(system.file("extdata/processedGSM.rda", package="ChIPseeker"))
    processedGSM <- get("processedGSM")
    newGSM <- c()

    gpldir <- "GPL"
    if (!file.exists(gpldir)) {
        dir.create(gpldir)
    }

    for (gid in gpl) {
        gg <- tryCatch(getGEO(gid, destdir=gpldir), error=function(e) NULL)
        if (is.null(gg)) {
            next
        }
        gsm <- Meta(gg)$sample_id
        gsm <- gsm[! (gsm %in% processedGSM) ]
        if (length(gsm) == 0) {
            next
        }
        newGSM <- c(newGSM, gsm)

        sf <- batchGetGSMsuppFile(gsm)
        if (!is.null(sf)) {
            save(sf, file=paste(gid, "_sf.rda", sep=""))
        }
    }

    processedGSM <- c(processedGSM, newGSM)
    processedGSM <- unique(processedGSM)
    save(processedGSM, file="../processedGSM.rda", compress="xz")


    sfiles <- list.files(pattern="_sf.rda")
    res <- data.frame(gsm=NULL, remoteFile=NULL)
    for (ff in sfiles) {
        load(ff)
        if (!is.null(sf)) {
            res <- rbind(res, sf)
        }
    }
    colnames(res)[2] <- "supplementary_file"


    GSMInfo <- lapply(unique(as.character(res$gsm)), function(i) {
        dbGetQuery(con,paste("select gsm,series_id,gpl,organism_ch1,title,characteristics_ch1,source_name_ch1,extract_protocol_ch1,description,data_processing,submission_date ",
                             "from gsm where gsm='", i, "'", sep=""))
    })

    GSMInfo <- do.call("rbind", GSMInfo)

    colnames(GSMInfo) <- sub("_ch1", "", colnames(GSMInfo))

    gsminfo <- merge(GSMInfo, res, by.x="gsm", by.y="gsm")

    tryCatch(utils::data("ucsc_release", package="ChIPseeker"))
    ucsc_release <- get("ucsc_release")

    genVer <- lapply(1:nrow(gsminfo), function(i)
                     getGenomicVersion(ucsc_release,
                                       gsminfo[i, "data_processing"],
                                       gsminfo[i, "organism"],
                                       gsminfo[i, "supplementary_file"])
                     )

    gsminfo$genomeVersion <- unlist(genVer)

    gse <- as.character(gsminfo$series_id)
    pubmed <- lapply(gse, function(i) {
        dbGetQuery(con,paste("select gse,pubmed_id ",
                             "from gse where gse='", i, "'", sep=""))
    })
    pm <- do.call(rbind, pubmed)
    pm <- unique(pm)
    gsminfo <- merge(gsminfo, pm, by.x="series_id", by.y="gse", all.x=TRUE)

    ## remove non-ASCII characters
    for(i in 1:ncol(gsminfo)) {
        gsminfo[,i] = iconv(gsminfo[,i], "latin1", "ASCII", sub="")
    }
    gsminfo2 <- gsminfo
    rm(gsminfo)

    data("gsminfo", package="ChIPseeker")
    gsminfo <- get("gsminfo")
    gsminfo <- rbind(gsminfo, gsminfo2)
    gsminfo <- unique(gsminfo)

    save(gsminfo, file="../gsminfo.rda", compress="xz")
}


getGenomicVersion <- function(ucsc_release, data_processing, organism, supplementary_file) {
    data_processing <- as.character(data_processing)
    organism <- as.character(organism)
    supplementary_file <- as.character(supplementary_file)

    species <- NULL
    gs <- subset(ucsc_release, subset = species == organism)
    if (nrow(gs) == 0) return(NA)

    genMatch <- unlist(sapply(gs$ucsc_version, grep, data_processing))
    if (length(genMatch) == 0) {
        genMatch <- unlist(sapply(gs$ucsc_version, grep, supplementary_file))
        if (length(genMatch) == 0) {
            return(NA)
        }
    }

    genVer <- names(genMatch)
    if (length(genVer) > 1) {
        genVer <- genVer[1]
    }

    return(genVer)
}

## getGSE_ENCODE <- function() {
##     encode=readLines("http://www.ncbi.nlm.nih.gov/geo/info/ENCODE.html")
##     encode.chipseq <- encode[grep("ChIP-Seq", encode)]
##     ## require(gsubfn)
##     ## gse <- sapply(encode.chipseq, function(i) {
##     ##     res <- strapply(i, "(GSE\\d+)")
##     ##     unique(unlist(res))
##     ## })
##     gse <- sapply(encode.chipseq, gsub, pattern='.*(GSE\\d+).*', replacement='\\1')
##     names(gse) <- NULL
##     return(gse)
## }

## GSE2GSM <- function(GSE) {
##     info <- getGEO(GSE, GSEMatrix=FALSE)
##     metaInfo <- Meta(info)
##     gsm <- metaInfo$sample_id
##     return(gsm)
## }

##' @importFrom parallel mclapply
##' @importFrom parallel detectCores
batchGetGSMsuppFile <- function(gsm) {
    suppfiles <- mclapply(seq_along(gsm), function(i) {
        cat("processing ", gsm[i], "\t",  i , " of ", length(gsm), "\n")
        tryCatch(getGSMsuppFile(gsm[i]), error=function(e) NULL)
    }, mc.cores=detectCores())

    suppfiles <- suppfiles[!unlist(lapply(suppfiles, is.null))]

    sf <- do.call("rbind", suppfiles)
    return(sf)
}

getGSMsuppFile <- function(GSM) {
    pkg <- "GEOquery"
    require(pkg, character.only=TRUE)
    getGEO <- eval(parse(text="getGEO"))
    Meta <- eval(parse(text="Meta"))

    destdir="geo_soft"
    if (!file.exists(destdir)) {
        dir.create(destdir)
    }
    info <- getGEO(GSM, GSEMatrix=FALSE, destdir=destdir)
    ## http://www.ncbi.nlm.nih.gov/geo/info/soft2.html
    metaInfo <- Meta(info)

    ## suppmentary file names
    fnames <- unlist(metaInfo[grep("supplementary", names(metaInfo))])
    names(fnames) <- NULL
    i <- c(grep("bed.gz", fnames), grep("Peak.gz", fnames), grep("bedGraph.gz", fnames))

    if (length(i) == 0) {
        message("No bed files found")
        return(NULL)
    }
    fnames <- fnames[i]
    res <- data.frame(gsm=GSM, remoteFile = fnames)
    return(res)
}

get_gsminfo <- function() {
    tryCatch(utils::data("gsminfo", package="ChIPseeker"))
    gsminfo <- get("gsminfo")
    return(gsminfo)
}
