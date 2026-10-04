library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

context("test function for enrichOverlap")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

make_peaks <- function(pos, width = 300) {
    GRanges("chr1", IRanges(pos, width = width), strand = "*")
}

test_that("enrichPeakOverlap accepts a single GRanges as target", {
    ## a bare GRanges was not wrapped into a list, so the permutation test
    ## failed with "GRanges objects don't support [[, as.list(), lapply()"
    set.seed(1)
    q <- make_peaks(sample(1000000:1400000, 5))
    t <- make_peaks(sample(1000000:1400000, 5))

    res <- enrichPeakOverlap(q, t, TxDb = txdb, nShuffle = 20,
                             mc.cores = 1, verbose = FALSE)

    expect_true(is.data.frame(res))
    expect_equal(nrow(res), 1)
    expect_equal(res$qLen, 5)
    expect_equal(res$tLen, 5)
    expect_true(res$N_OL >= 0 && res$N_OL <= 5)

    ## a list of GRanges keeps working as well
    res2 <- enrichPeakOverlap(q, list(t), TxDb = txdb, nShuffle = 20,
                              mc.cores = 1, verbose = FALSE)
    expect_equal(res2$N_OL, res$N_OL)
})

test_that("the overlap count is symmetric while the tested ratio is not", {
    ## issue #84: exchanging queryPeak and targetPeak reports the same N_OL but
    ## the tested quantity differs, because the observed ratio is normalised by
    ## the target size and the target is the set that gets shuffled
    set.seed(11)
    hot <- sample(1000000:1400000, 3)
    A <- make_peaks(c(sample(hot, 2), sample(2000000:2400000, 18)))
    B <- make_peaks(c(sample(hot, 3), sample(3000000:3400000, 97)))

    set.seed(1)
    ab <- enrichPeakOverlap(A, B, TxDb = txdb, nShuffle = 100, mc.cores = 1,
                            verbose = FALSE)
    set.seed(1)
    ba <- enrichPeakOverlap(B, A, TxDb = txdb, nShuffle = 100, mc.cores = 1,
                            verbose = FALSE)

    ## the number of overlapping peaks is direction free ...
    expect_equal(ab$N_OL, ba$N_OL)
    expect_equal(ab$N_OL, length(intersect(A, B)))
    expect_equal(ab$qLen, ba$tLen)
    expect_equal(ab$tLen, ba$qLen)

    ## ... while the ratio the p-value is based on is target normalised, so it
    ## changes when the two peak sets differ in size
    expect_false(isTRUE(all.equal(ab$N_OL / ab$tLen, ba$N_OL / ba$tLen)))

    ## the p-values stay valid permutation p-values
    expect_true(all(c(ab$pvalue, ba$pvalue) > 0))
    expect_true(all(c(ab$pvalue, ba$pvalue) <= 1))
})