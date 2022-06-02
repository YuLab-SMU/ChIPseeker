library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

context("test plotTagMatrix() and related functions")

test_that("test plotPeakProf2 function use txdb",{
  
  peak <- getSampleFiles()[[4]]
  peak_list <- getSampleFiles()[4:5]
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  # test single peak file
  p1_1 <- plotPeakProf(peak = peak,
                       upstream = 1000,
                       downstream = 1000,
                       by = "gene",
                       type = "start_site",
                       TxDb = txdb,
                       nbin = 800)
  
  expect_is(p1_1, "gg")
  
  # test a list of peak files
  p1_2 <- plotPeakProf(peak = peak_list,
                       upstream = 1000,
                       downstream = 1000,
                       by = "gene",
                       type = "start_site",
                       TxDb = txdb,
                       nbin = 800)
  
  expect_is(p1_2, "gg")
  
  # test body region
  # without extension
  p2_1 <- plotPeakProf(peak = peak_list,
                       by = "gene",
                       type = "body",
                       TxDb = txdb,
                       nbin = 800)
  
  expect_is(p2_1, "gg")
  
  # extend with rel object
  p2_2 <- plotPeakProf(peak = peak_list,
                       by = "gene",
                       type = "body",
                       TxDb = txdb,
                       upstream = rel(0.2),
                       downstream = rel(0.2),
                       nbin = 800)
  
  expect_is(p2_2, "gg")
  
  # extend with actual number
  p2_3 <- plotPeakProf(peak = peak_list,
                       by = "gene",
                       type = "body",
                       TxDb = txdb,
                       upstream = 1000,
                       downstream = 1000,
                       nbin = 800)
  
  expect_is(p2_3, "gg")
  
})


test_that("test plotPeakProf2 function use self-made granges",{
  
  peak <- getSampleFiles()[[4]]
  peak_list <- getSampleFiles()[4:5]
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  # we consider transcript region as enhancer region
  # and make self-made granges object
  # they can be the same in the form of granges object
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  enhancer <- transcripts(txdb)[1:5000,]
  
  # test single peak file
  p1_1 <- plotPeakProf(peak = peak,
                       upstream = 1000,
                       downstream = 1000,
                       by = "gene",
                       type = "start_site",
                       TxDb = enhancer,
                       nbin = 800)
  
  expect_is(p1_1, "gg")
  
  # test a list of peak files
  p1_2 <- plotPeakProf(peak = peak_list,
                       upstream = 1000,
                       downstream = 1000,
                       by = "gene",
                       type = "start_site",
                       TxDb = enhancer,
                       nbin = 800)
  
  expect_is(p1_2, "gg")
  
  # test body region
  # without extension
  p2_1 <- plotPeakProf(peak = peak_list,
                       by = "gene",
                       type = "body",
                       TxDb = enhancer,
                       nbin = 800)
  
  expect_is(p2_1, "gg")
  
  # extend with rel object
  p2_2 <- plotPeakProf(peak = peak_list,
                       by = "gene",
                       type = "body",
                       TxDb = enhancer,
                       upstream = rel(0.2),
                       downstream = rel(0.2),
                       nbin = 800)
  
  expect_is(p2_2, "gg")
  
  # extend with actual number
  p2_3 <- plotPeakProf(peak = peak_list,
                       by = "gene",
                       type = "body",
                       TxDb = enhancer,
                       upstream = 1000,
                       downstream = 1000,
                       nbin = 800)
  
  expect_is(p2_3, "gg")
  
})