##' Get gene annotation from annotation database
##'
##' This function retrieves gene annotations (symbol, gene name, etc.) from
##' annotation databases (e.g., org.Hs.eg.db) based on gene IDs.
##'
##' @description
##' The function queries an annotation database to retrieve gene information such
##' as gene symbols, full gene names, Ensembl IDs, and other identifiers. It
##' handles both Entrez Gene IDs and Ensembl Gene IDs, and returns a data frame
##' with the requested annotation columns.
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Validates and converts the gene ID type to the appropriate keytype
##'         for the annotation database (ENTREZID or ENSEMBL)
##'   \item Removes version numbers from Ensembl IDs (e.g., "ENSG00000123456.1"
##'         becomes "ENSG00000123456")
##'   \item Queries the annotation database using \code{AnnotationDbi::select()}
##'   \item Handles duplicate gene IDs by keeping only the first occurrence
##'   \item Returns a data frame with rows matching the input gene IDs
##' }
##'
##' @param annoDb character, name of the annotation database package (e.g.,
##'   "org.Hs.eg.db" for human). The package must be installed and loaded
##' @param geneID character vector of gene IDs to query. Can be Entrez Gene IDs
##'   or Ensembl Gene IDs depending on the \code{type} parameter
##' @param type character, gene ID type. Must be one of "Entrez Gene ID" or
##'   "Ensembl gene ID" / "Ensembl Gene ID"
##' @param columns character vector, names of columns to retrieve from the
##'   annotation database. Common values include "SYMBOL", "GENENAME", "ENSEMBL",
##'   "ENTREZID", "GO", "PATH", etc. See the annotation package documentation
##'   for available columns
##' @return A data.frame with:
##'   \itemize{
##'     \item Rows corresponding to the input gene IDs (in the same order)
##'     \item Columns as specified in the \code{columns} parameter, plus the
##'           keytype column (ENTREZID or ENSEMBL)
##'     \item NA values for gene IDs that could not be found in the database
##'   }
##'   Returns NA if the ID type is not supported or if the database query fails
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


##' Add gene annotation to GRanges object
##'
##' This function adds gene annotation columns (symbol, gene name, etc.) to a
##' GRanges object containing peaks with gene IDs.
##'
##' @description
##' The function retrieves gene annotations from an annotation database based on
##' the gene IDs stored in the GRanges object, then adds the annotation columns
##' as metadata columns to the GRanges object.
##'
##' @details
##' The function:
##' \enumerate{
##'   \item Extracts gene IDs from the \code{geneId} metadata column of the
##'         GRanges object
##'   \item Calls \code{getGeneAnno()} to retrieve annotations from the database
##'   \item Adds each annotation column as a new metadata column to the GRanges
##'         object (excluding the keytype column which already exists as geneId)
##' }
##'
##' @param peak.gr GRanges object containing peaks with a \code{geneId} metadata
##'   column containing gene IDs
##' @param annoDb character, name of the annotation database package (e.g.,
##'   "org.Hs.eg.db")
##' @param type character, gene ID type. Must be one of "Entrez Gene ID" or
##'   "Ensembl gene ID" / "Ensembl Gene ID"
##' @param columns character vector, names of columns to retrieve from the
##'   annotation database (e.g., "SYMBOL", "GENENAME", "ENSEMBL")
##' @return GRanges object with additional metadata columns containing gene
##'   annotations. The original GRanges object is returned if annotation retrieval
##'   fails
##' @seealso \code{\link{getGeneAnno}} for the underlying annotation retrieval
##'   function
##' @author G Yu
##' @noRd
addGeneAnno <- function(peak.gr, annoDb, type, columns) {
    geneAnno <- getGeneAnno(annoDb, peak.gr$geneId, type, columns)
    if (! all(is.na(geneAnno))) {
        for(cn in colnames(geneAnno)[-1]) {
            mcols(peak.gr)[[cn]] <- geneAnno[, cn]
        }
    }
    return(peak.gr)
}

