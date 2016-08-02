
##' get gene annotation, symbol, gene name etc.
##'
##' 
##' @title getGeneAnno 
##' @param annoDb annotation package
##' @param geneID query geneID
##' @param type gene ID type
##' @return data.frame
##' @importFrom AnnotationDbi select
##' @author G Yu
getGeneAnno <- function(annoDb, geneID, type){
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
                                keys=kk[i],
                                keytype=kt,
                                columns=c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))),
        error = function(e) NULL)

    if (is.null(ann)) {
        warning("ID type not matched, gene annotation will not be added...")
        return(NA)
    }
    idx <- getFirstHitIndex(ann[,kt])
    ann <- ann[idx,]

    idx <- unlist(sapply(kk, function(x) which(x==ann[,kt])))

    res <- matrix(NA, ncol=ncol(ann), nrow=length(kk)) %>% as.data.frame
    colnames(res) <- colnames(ann)
    res[i,] <- ann[idx,]
    return(res)
}


addGeneAnno <- function(peak.gr, annoDb, type) {
    geneAnno <- getGeneAnno(annoDb, peak.gr$geneId, type)
    if (! all(is.na(geneAnno))) {
        for(cn in colnames(geneAnno)[-1]) {
            mcols(peak.gr)[[cn]] <- geneAnno[, cn]
        }
    }
    return(peak.gr)
}

