library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

context("test function for getAllFlankingGene")

test_that("flank distances are 0 for overlaps and TSS based otherwise", {
    ## issue #235: a distance of 0 means that the peak overlaps the feature
    ## range, while non-overlapping peaks are assigned the (signed) distance
    ## to the feature TSS.
    features <- GRanges("chr1",
                        IRanges(start = c(1000, 20000), width = 3000),
                        strand = "+",
                        gene_id = c("G1", "G2"))

    ## peak 1 overlaps the body of G1, peak 2 is 950bp upstream of G2
    peaks <- GRanges("chr1", IRanges(start = c(1500, 19000),
                                     end = c(1600, 19050)))

    res <- getAllFlankingGene(peaks, features, level = "gene",
                              distance = 2000)

    res <- res[order(res$peakIdx), ]
    expect_equal(res$peakIdx, c(1L, 2L))
    expect_equal(res$flank_geneIds, c("G1", "G2"))
    expect_equal(res$flank_gene_distances, c("0", "-950"))
})

test_that("annotatePeak reports flank distances following the same rule", {
    ## a peak inside a transcript body gets distance 0 even though its
    ## distance to the TSS is large
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    peak <- GRanges("chr1", IRanges(243832958, 243833008))

    pa <- annotatePeak(peak, TxDb = txdb, level = "transcript",
                       addFlankGeneInfo = TRUE, flankDistance = 5000,
                       verbose = FALSE)
    df <- as.data.frame(pa)

    expect_true(is.character(df$flank_geneIds))
    expect_true(is.character(df$flank_gene_distances))

    distances <- as.numeric(strsplit(df$flank_gene_distances, ";")[[1]])
    expect_true(any(distances == 0))
    expect_true(all(!is.na(distances)))
})
