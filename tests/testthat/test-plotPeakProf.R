library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

context("test plotPeakProf() for a list of windows")

peak <- getSampleFiles()[[4]]
peak_list <- getSampleFiles()[4:5]
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## self-made enhancer region in the form of granges object
enhancer <- transcripts(txdb)[1:5000,]

## self-made non-enhancer region in the form of granges object
non_enhancer <- unlist(fiveUTRsByTranscript(txdb))[1:5000]

gr <- list(enhancer,non_enhancer)

test_that("input two self-made granges object",{
  
  p <- plotPeakProf(peak = peak,
                    conf = 0.95,
                    by = c("enhancer","non-enhancer"),
                    windows_name = c("enhancer","non-enhancer"),
                    weightCol = "V5",
                    type = "start_site",
                    upstream = 1000,
                    downstream = 1000,
                    TxDb = list(enhancer,non_enhancer))
  
  expect_is(p,"gg")
  
})

test_that("input a list of peaks",{
  
  p <- plotPeakProf(peak = peak_list,
                    TxDb = list(enhancer,non_enhancer),
                    conf = 0.95,
                    by = c("enhancer","non-enhancer"),
                    windows_name = c("enhancer","non-enhancer"),
                    weightCol = "V5",
                    type = "start_site",
                    upstream = 1000,
                    downstream = 1000)
  
  expect_is(p,"gg")
})


test_that("input gr and txdb input",{
  
  p <- plotPeakProf(peak = peak,
                    TxDb = list(enhancer,txdb),
                    conf = 0.95,
                    by = c("enhancer","gene"),
                    windows_name = c("enhancer","gene"),
                    weightCol = "V5",
                    type = "start_site",
                    upstream = 1000,
                    downstream = 1000)
  
  expect_is(p,"gg")
})


test_that("check body region",{
  
  p <- plotPeakProf(peak = peak,
                    TxDb = list(enhancer,txdb),
                    conf = 0.95,
                    by = c("enhancer","gene"),
                    windows_name = c("enhancer","gene"),
                    weightCol = "V5",
                    type = "body",
                    upstream = 1000,
                    downstream = 1000,
                    nbin = 800)
  
  expect_is(p,"gg")
  
  p <- plotPeakProf(peak = peak,
                    TxDb = list(enhancer,txdb),
                    conf = 0.95,
                    by = c("enhancer","gene"),
                    windows_name = c("enhancer","gene"),
                    weightCol = "V5",
                    type = "body",
                    nbin = 800)
  
  expect_is(p,"gg")
  
  p <- plotPeakProf(peak = peak,
                    TxDb = list(enhancer,txdb),
                    conf = 0.95,
                    by = c("enhancer","gene"),
                    windows_name = c("enhancer","gene"),
                    weightCol = "V5",
                    type = "body",
                    upstream = rel(0.2),
                    downstream = rel(0.2),
                    nbin = 800)
  
  expect_is(p,"gg")
  
  p <- plotPeakProf(peak = peak_list,
                    TxDb = list(enhancer,txdb),
                    conf = 0.95,
                    by = c("enhancer","gene"),
                    windows_name = c("enhancer","gene"),
                    weightCol = "V5",
                    type = "body",
                    upstream = rel(0.2),
                    downstream = rel(0.2),
                    nbin = 800)
  
  expect_is(p,"gg")
})