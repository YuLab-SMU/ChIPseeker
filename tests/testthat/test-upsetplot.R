library(ChIPseeker)
library(GenomicRanges)

context("test function for upsetplot")

## a minimal csAnno object, so that the test does not depend on a TxDb
make_csAnno <- function(n = 60) {
    detail <- data.frame(
        genic = rep(c(TRUE, FALSE), length.out = n),
        Intergenic = rep(c(FALSE, TRUE), length.out = n),
        Promoter = rep(c(TRUE, FALSE), length.out = n),
        fiveUTR = FALSE,
        threeUTR = FALSE,
        Exon = rep(c(TRUE, FALSE), length.out = n),
        Intron = rep(c(FALSE, TRUE), length.out = n),
        downstream = FALSE,
        distal_intergenic = FALSE
    )
    anno <- GRanges("chr1", IRanges(seq_len(n) * 1000, width = 500),
                    strand = "*")
    mcols(anno)$distanceToTSS <- seq_len(n) * 10

    methods::new(
        "csAnno",
        anno = anno,
        tssRegion = c(-3000, 3000),
        level = "transcript",
        hasGenomicAnnotation = TRUE,
        detailGenomicAnnotation = detail,
        annoStat = as.data.frame(table(detail$Exon)),
        peakNum = n
    )
}

test_that("upsetplot can be drawn, also with the vennpie sub-view", {
    skip_if_not_installed("ggupset")
    skip_if_not_installed("ggplot2")

    cs <- make_csAnno()
    fn <- tempfile(fileext = ".png")
    on.exit(unlink(fn), add = TRUE)

    ## the plain UpSet plot has to be drawable ...
    expect_error(
        suppressWarnings(ggplot2::ggsave(
            fn, upsetplot(cs), width = 7, height = 5, dpi = 72
        )),
        NA
    )

    ## ... and so has the one with the vennpie sub-view, which is embedded as an
    ## annotation_custom() layer and therefore only works below
    ## coord_cartesian() with ggplot2 >= 4.0
    skip_if_not_installed("ggimage")
    expect_error(
        suppressWarnings(ggplot2::ggsave(
            fn, upsetplot(cs, vennpie = TRUE), width = 7, height = 5, dpi = 72
        )),
        NA
    )
})