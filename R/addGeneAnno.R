##' get gene annotation, symbol, gene name etc.
##'
##'
##' @title getGeneAnno
##' @param annoDb annotation package
##' @param geneID query geneID
##' @param type gene ID type
##' @param columns names of columns to be obtained from database
##' @return data.frame
##' @importFrom AnnotationDbi select
##' @author G Yu
getGeneAnno <- function(annoDb, geneID, type, columns){
    kk <- unlist(geneID)
    require(annoDb, character.only = TRUE)
    annoDb <- eval(parse(text=annoDb))

    if (type == "Entrez Gene ID") {
        kt <- "ENTREZID"
    } else if (type =="Ensembl gene ID" || type == "Ensembl Gene ID") {
        kt <- "ENSEMBL"
    } else {
        message("geneID type is not supported...\tPlease report it to developer...\n")
        return(NA)
    }

    i <- which(!is.na(kk))
    kk <- gsub("\\.\\d+$", "", kk)
    ann <- tryCatch(
        suppressWarnings(select(annoDb,
                                keys=unique(kk[i]),
                                keytype=kt,
                                columns=columns)),
        error = function(e) NULL)

    if (is.null(ann)) {
        warning("ID type not matched, gene annotation will not be added...")
        return(NA)
    }
    idx <- getFirstHitIndex(ann[,kt])
    ann <- ann[idx,]

    ## idx <- unlist(sapply(kk, function(x) which(x==ann[,kt])))
    ## res <- matrix(NA, ncol=ncol(ann), nrow=length(kk)) %>% as.data.frame
    ## colnames(res) <- colnames(ann)
    ## res[i,] <- ann[idx,]

    rownames(ann) <- ann[, kt]
    res <- ann[as.character(kk),]

    return(res)
}


addGeneAnno <- function(peak.gr, annoDb, type, columns) {
    geneAnno <- getGeneAnno(annoDb, peak.gr$geneId, type, columns)
    if (! all(is.na(geneAnno))) {
        for(cn in colnames(geneAnno)[-1]) {
            mcols(peak.gr)[[cn]] <- geneAnno[, cn]
        }
    }
    return(peak.gr)
}

