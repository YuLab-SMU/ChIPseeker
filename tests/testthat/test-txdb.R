library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(yulab.utils)

context("TXDB")

test_that("Update txdb", {
    hg19_txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    ChIPseeker:::.ChIPseekerEnv(hg19_txdb)
    expect_equal(ChIPseeker:::get_env_genome(), "hg19")

    hg38_txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    ChIPseeker:::.ChIPseekerEnv(hg38_txdb)
    expect_equal(ChIPseeker:::get_env_genome(), "hg38")
})


test_that("txdb", {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    ChIPseeker:::.ChIPseekerEnv(txdb)
    expect_equal(ChIPseeker:::IDType(txdb), "Entrez Gene ID")

    ## The tx_id -> (tx_name, gene_id) mapping is annotation data that is
    ## rebuilt for every Bioconductor release. Pinning a snapshot made this
    ## test fail on the Bioc 3.24 builders: tx 70455 used to be
    ## "ENST00000487630.1_3" (no gene id, so TXID2EG() fell back to the
    ## transcript name) and is now "ENST00000959756.2_2/3554". Take the ids
    ## and derive the expected values from the installed TxDb, so the test
    ## follows the annotation instead of a particular release. Cover both a
    ## transcript with a gene id and one without (the tx_name fallback).
    tx <- GenomicFeatures::transcripts(txdb, columns = c("tx_id", "gene_id"))
    ngene <- lengths(tx$gene_id)
    with_gene <- which(ngene == 1)
    without_gene <- which(ngene == 0)
    txids <- character()
    if (length(with_gene))
        txids <- c(txids, as.character(tx$tx_id[with_gene[1]]))
    if (length(without_gene))
        txids <- c(txids, as.character(tx$tx_id[without_gene[1]]))

    map <- suppressMessages(AnnotationDbi::select(
        txdb, keys = txids, columns = c("TXNAME", "GENEID"), keytype = "TXID"
    ))
    i <- match(txids, as.character(map$TXID))
    txname <- as.character(map$TXNAME[i])
    geneid <- as.character(map$GENEID[i])
    ## TXID2EG() uses the transcript name when no gene id is available
    fallback <- is.na(geneid)
    geneid[fallback] <- txname[fallback]

    expect_equal(ChIPseeker:::TXID2EG(txids),
                 paste(txname, geneid, sep = "/"))
    expect_equal(ChIPseeker:::TXID2EG(txids, geneIdOnly = TRUE),
                 geneid)
})
