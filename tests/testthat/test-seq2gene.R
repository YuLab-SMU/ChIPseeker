library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

context("test function for seq2gene")

test_that("seq2gene handles regions without exon overlap", {
    ## issue #248: exons is NA when no region overlaps an exon, and
    ## `exons$gene` used to fail with
    ## "$ operator is invalid for atomic vectors"
    ## (chr1:100000-100100 overlaps an intron but no exon)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    peak <- GRanges("chr1", IRanges(100000, 100100))

    genes <- seq2gene(peak, tssRegion = c(-3000, 3000),
                      flankDistance = 5000, TxDb = txdb)

    expect_true(is.character(genes))
    expect_false(any(is.na(genes)))
    ## the intron host gene is reported
    expect_true(length(genes) > 0)
})

test_that("seq2gene returns an empty vector when nothing is nearby", {
    ## a peak in the chr1 centromere without any exon/intron overlap and
    ## far away from every gene
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    peak <- GRanges("chr1", IRanges(122000000, 122000100))

    genes <- seq2gene(peak, tssRegion = c(0, 0),
                      flankDistance = 1, TxDb = txdb)

    expect_equal(genes, character(0))
})
