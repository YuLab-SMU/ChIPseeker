# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# context("test getTagMatrix() and related functions")

# test_that("getBioRegion function", {
  
#   # test three kinds of regions derived from getBioRegion()
#   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
#   gene_start <- getBioRegion(TxDb = txdb,
#                              upstream = 1000,
#                              downstream = 1000,
#                              by = 'gene',
#                              type = "start_site")
#   expect_is(gene_start,"GRanges")
  
#   gene_end <- getBioRegion(TxDb = txdb,
#                            upstream = 1000,
#                            downstream = 1000,
#                            by = 'gene',
#                            type = "end_site")
#   expect_is(gene_end,"GRanges")
  
#   gene_body <- getBioRegion(TxDb = txdb,
#                             upstream = 1000,
#                             downstream = 1000,
#                             by = 'gene',
#                             type = "body")
#   expect_is(gene_body,"GRanges")
  
#   })

# test_that("getPromoters functions",{
  
#   # test two kinds of regions derived from getPromoters
#   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
#   gene <- getPromoters(TxDb=txdb,
#                        upstream=1000,
#                        downstream=1000,
#                        by = "gene")
  
#   transcript <- getPromoters(TxDb=txdb,
#                              upstream=1000,
#                              downstream=1000,
#                              by = "transcript")
  
#   expect_is(gene,"GRanges")
#   expect_is(transcript,"GRanges")
  
# })

# test_that("makeBioRegionFromGranges function",{
  
#   # we consider transcript region as enhancer region
#   # and make self-made granges object
#   # they can be the same in the form of granges object
#   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#   enhancer <- transcripts(txdb)[1:5000,]
  
#   ## we test three kinds of region, start_site, end_site and body
#   enhancer_body <- makeBioRegionFromGranges(gr = enhancer,
#                                             by = "enhancer",
#                                             type = "body")
  
#   enhancer_start <- makeBioRegionFromGranges(gr = enhancer,
#                                              by = "enhancer",
#                                              type = "start_site",
#                                              upstream = 1000,
#                                              downstream = 1000)
  
#   enhancer_end <- makeBioRegionFromGranges(gr = enhancer,
#                                            by = "enhancer",
#                                            type = "end_site",
#                                            upstream = 1000,
#                                            downstream = 1000)
  
#   expect_is(enhancer_body,"GRanges")
#   expect_is(enhancer_start,"GRanges")
#   expect_is(enhancer_end,"GRanges")
  
#   ## test the label
#   expect_equal(attr(enhancer_body,'label'),c("enhancer_SS","enhancer_TS"))
#   expect_equal(attr(enhancer_start,'label'),"enhancer_SS")
#   expect_equal(attr(enhancer_end,'label'),"enhancer_TS")
  
# })

# test_that("getTagMatrix function for single peak file",{
  
#   peak <- getSampleFiles()[[4]]
#   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
#   # make the window by getBioRegion()
#   gene_start <- getBioRegion(TxDb = txdb,
#                              upstream = 1000,
#                              downstream = 1000,
#                              by = 'gene',
#                              type = "start_site")
  
#   # make the window by makeBioRegionFromGranges()
#   enhancer <- transcripts(txdb)[1:5000,]
  
#   enhancer_body <- makeBioRegionFromGranges(gr = enhancer,
#                                             by = "enhancer",
#                                             type = "body")
  
#   # test input window parameter
#   mt1 <- getTagMatrix(peak = peak,
#                       windows = gene_start,
#                       weightCol = "V5")
  
#   expect_is(mt1, "matrix")
  
#   # without extending flank
#   mt2_1 <- getTagMatrix(peak = peak,
#                         windows = enhancer_body,
#                         weightCol = "V5",
#                         nbin = 800)
  
#   expect_is(mt2_1, "matrix")
  
#   # extend flank by rel object
#   mt2_2 <- getTagMatrix(peak = peak,
#                         windows = enhancer_body,
#                         weightCol = "V5",
#                         upstream = rel(0.2),
#                         downstream = rel(0.2),
#                         nbin = 800)
  
#   expect_is(mt2_2, "matrix")
  
#   # extend flank by actual number
#   mt2_3 <- getTagMatrix(peak = peak,
#                         windows = enhancer_body,
#                         weightCol = "V5",
#                         upstream = 1000,
#                         downstream = 1000,
#                         nbin = 800)
  
#   expect_is(mt2_3, "matrix")
  
#   # test input without window parameter 
  
#   # make window through txdb object
#   mt3 <- getTagMatrix(peak = peak,
#                       weightCol = "V5",
#                       TxDb = txdb,
#                       by = "gene",
#                       type = "start_site",
#                       upstream = 3000,
#                       downstream = 3000)
  
#   expect_is(mt3, "matrix")
  
#   # make window through self-made grange object
#   mt4 <- getTagMatrix(peak = peak,
#                       weightCol = "V5",
#                       TxDb = enhancer,
#                       by = "gene",
#                       type = "start_site",
#                       upstream = 1000,
#                       downstream = 1000)
  
#   expect_is(mt4, "matrix")
  
#   # without extending flank
#   mt5_1 <- getTagMatrix(peak = peak,
#                         weightCol = "V5",
#                         TxDb = txdb,
#                         by = "gene",
#                         type = "body",
#                         nbin = 800)
  
#   expect_is(mt5_1, "matrix")
  
#   # extend flank by rel object
#   mt5_2 <- getTagMatrix(peak = peak,
#                         TxDb = enhancer,
#                         weightCol = "V5",
#                         by = "enhancer",
#                         type = "body",
#                         upstream = rel(0.2),
#                         downstream = rel(0.2),
#                         nbin = 800)
  
#   expect_is(mt5_2, "matrix")
  
#   # extend flank by actual number
#   mt5_3 <- getTagMatrix(peak = peak,
#                         TxDb = txdb,
#                         weightCol = "V5",
#                         by = "gene",
#                         type = "body",
#                         upstream = 1000,
#                         downstream = 1000,
#                         nbin = 800)
  
#   expect_is(mt5_3, "matrix")
  
# })
