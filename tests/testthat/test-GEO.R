library(ChIPseeker)

context("test function for GEO data mining")

test_that("GEO supplementary file urls are rewritten to https", {
    ## issue #254: GEO does not serve these files over ftp:// any more
    url <- c("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1nnn/GSM1/suppl/a.bed.gz",
             "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2nnn/GSM2/suppl/b.bed.gz")

    expect_equal(geoHttpsUrl(url),
                 c("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1nnn/GSM1/suppl/a.bed.gz",
                   "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2nnn/GSM2/suppl/b.bed.gz"))
})

test_that("downloadGEO.internal skips files that are already downloaded", {
    ## offline test: the destination file exists, so nothing is downloaded
    destDir <- tempdir()
    fname <- "GSM288348_Smad1.bed.gz"
    destfile <- file.path(destDir, fname)
    file.create(destfile)
    on.exit(unlink(destfile), add = TRUE)

    info <- data.frame(
        supplementary_file = paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/",
                                    "GSM288nnn/GSM288348/suppl/", fname))

    expect_silent(downloadGEO.internal(info, destDir))
})
