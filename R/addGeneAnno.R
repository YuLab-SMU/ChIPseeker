
##' add gene annotation, symbol, gene name etc.
##'
##' 
##' @title addGeneAnno 
##' @param annoDb annotation package
##' @param geneID query geneID
##' @param type gene ID type
##' @return data.frame
##' @importFrom AnnotationDbi select
##' @author G Yu
addGeneAnno <- function(annoDb, geneID, type){
    kk <- unlist(geneID)
    require(annoDb, character.only = TRUE)
    annoDb <- eval(parse(text=annoDb))
    if (type == "Entrez Gene ID") {
        kt <- "ENTREZID"
    } else if (type =="Ensembl gene ID" || type == "Ensembl Gene ID") {
        kt <- "ENSEMBL"
    } else {
        warnings("geneID type is not supported...\tPlease report it to developer...\n")
        return(NA)
    }

    ann <- suppressWarnings(select(annoDb,
                                   keys=kk,
                                   keytype=kt,
                                   columns=c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")))
    idx <- getFirstHitIndex(ann[,kt])
    ann <- ann[idx,]

    idx <- unlist(sapply(kk, function(x) which(x==ann[,kt])))
    ann <- ann[idx,]
    return(ann)
}

