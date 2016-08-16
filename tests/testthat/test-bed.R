library(ChIPseeker)

context("bed file")

test_that("parse bed file", {
    files <- getSampleFiles()
    for (i in seq_along(files)) {
        expect_true(is(readPeakFile(files[[i]]), "GRanges"))
    }
})

