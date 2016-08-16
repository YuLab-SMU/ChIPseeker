library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)

context("TXDB")

test_that("txdb", {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    expect_equal(ChIPseeker:::IDType(txdb), "Entrez Gene ID")
    expect_equal(ChIPseeker:::TXID2EG("70455"), "uc002qsd.4/1")
    expect_equal(ChIPseeker:::TXID2EG("70455", geneOnly=TRUE), "1")
}

