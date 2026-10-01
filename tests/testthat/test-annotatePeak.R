library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peaks <- GRanges(
    c("chr1", "chr2", "chr20"),
    IRanges(start = c(1000000, 3000000, 3000000), width = 1000)
)

test_that("geneChr / geneStrand are characters, not factor codes", {
    ## issues #233, #247:
    ## as.data.frame() turns 'seqnames' / 'strand' into factors and assigning
    ## a factor into mcols() dropped the class, so geneChr / geneStrand came out
    ## as integers (e.g. 1/2 instead of chr1/chr2 and +/-)
    pa <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000),
                       level = "transcript", verbose = FALSE)
    m <- mcols(as.GRanges(pa))

    expect_true(is.character(m$geneChr))
    expect_true(is.character(m$geneStrand))
    expect_true(all(m$geneStrand %in% c("+", "-")))
    expect_true(all(grepl("^chr", m$geneChr)))

    ## numeric columns must stay numeric
    expect_true(is.integer(m$geneStart))
    expect_true(is.integer(m$geneEnd))
})

test_that("peaks without any feature in TxDb are reported, not silently dropped", {
    ## issue #251: peaks on contigs/scaffolds that carry no gene used to be
    ## dropped without any message at all
    peaks2 <- GRanges(
        c("chr1", "chr1_gl000191_random", "chrM"),
        IRanges(start = c(1000000, 5000, 5000), width = 1000)
    )

    expect_warning(
        annotatePeak(peaks2, TxDb = txdb, tssRegion = c(-3000, 3000),
                     level = "transcript", verbose = FALSE),
        "peaks were dropped"
    )
})

test_that("transcript annotation metadata follows the overlapping isoform", {
    ## issue #252: a peak can overlap one isoform while another nested
    ## transcript has the closer TSS.  All transcript-level fields must then
    ## describe the isoform that supplied the genomic annotation.
    features <- GRanges(
        "chr1",
        IRanges(start = c(100, 500), end = c(1000, 900)),
        strand = c("+", "-"),
        tx_id = c(101L, 202L),
        gene_id = c("GENE1", "GENE1")
    )
    peak <- GRanges("chr1", IRanges(700, 710), strand = "*")

    aligned <- ChIPseeker:::.alignTranscriptAnnotation(
        peak, features, index = 1L, distance = 999,
        annotationFeatureId = "202"
    )

    expect_equal(aligned$index, 2L)
    expect_equal(as.character(mcols(features)$gene_id[aligned$index]), "GENE1")
    expect_equal(as.integer(mcols(features)$tx_id[aligned$index]), 202L)
    expect_equal(aligned$distance, 190)
})

test_that("gene level annotation is aligned by gene_id", {
    ## issue #252 (gene level): .alignAnnotationFeature() also aligns the
    ## nearest gene with the gene of an exon/intron hit
    features <- GRanges(
        "chr1",
        IRanges(start = c(100, 500), end = c(1000, 900)),
        strand = c("+", "-"),
        gene_id = c("GENE1", "GENE2")
    )
    peak <- GRanges("chr1", IRanges(700, 710), strand = "*")

    aligned <- ChIPseeker:::.alignAnnotationFeature(
        peak, features, index = 1L, distance = 999,
        annotationFeatureId = "GENE2", idColumn = "gene_id"
    )

    expect_equal(aligned$index, 2L)
    expect_equal(aligned$distance, 190)

    ## unknown ids leave the nearest feature untouched
    unchanged <- ChIPseeker:::.alignAnnotationFeature(
        peak, features, index = 1L, distance = 999,
        annotationFeatureId = NA_character_, idColumn = "gene_id"
    )
    expect_equal(unchanged$index, 1L)
    expect_equal(unchanged$distance, 999)
})

test_that("gene-level annotation and geneId describe the same gene", {
    ## issue #252 (gene level): the annotation comes from an exon/intron hit
    ## while the nearest gene is found independently, so both have to follow
    ## the gene of the annotation
    g <- suppressMessages(genes(txdb))
    gs <- g[as.character(seqnames(g)) == "chr17" & width(g) > 20000]
    strandG <- as.character(strand(gs))
    tss <- ifelse(strandG == "+", start(gs), end(gs))
    pos <- round(tss + rep(c(0.15, 0.4, 0.6, 0.85), each = length(gs)) *
                     width(gs) * ifelse(strandG == "+", 1, -1))
    own <- ifelse(strandG == "+", 1, -1) * (pos - tss)

    ## keep the peaks for which the closest TSS belongs to another gene
    d <- abs(outer(pos, tss, "-"))
    diag(d) <- Inf
    cand <- unique(pos[apply(d, 1, min) < abs(own)])
    expect_true(length(cand) > 0)

    peaks <- GRanges("chr17", IRanges(utils::head(cand, 20), width = 200))
    pa <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000),
                       level = "gene", verbose = FALSE)
    df <- as.data.frame(pa)

    ## annotation, gene coordinates and geneId describe the same gene
    m <- match(as.character(df$geneId), names(g))
    expect_false(any(is.na(m)))
    expect_equal(as.character(df$geneId), as.character(names(g)[m]))
    expect_equal(df$geneChr, as.character(seqnames(g))[m])
    expect_equal(df$geneStart, start(g)[m])
    expect_equal(df$geneEnd, end(g)[m])
    expect_equal(as.character(df$geneStrand), as.character(strand(g))[m])

    strandRep <- as.character(strand(g))[m]
    tssRep <- as.numeric(ifelse(strandRep == "+", start(g)[m], end(g)[m]))
    dStart <- ifelse(strandRep == "+", start(peaks) - tssRep,
                     tssRep - start(peaks))
    dEnd <- ifelse(strandRep == "+", end(peaks) - tssRep, tssRep - end(peaks))
    expect_equal(df$distanceToTSS, ifelse(abs(dStart) <= abs(dEnd),
                                          dStart, dEnd))

    ## genic rows belong to the reported gene and match the annotation
    genic <- grepl("^(Exon|Intron)", df$annotation)
    expect_true(any(genic))
    expect_true(all(overlapsAny(peaks[genic], g[m[genic]])))
    annGene <- sub("^[A-Za-z']+ \\([^/]+/([^,]+),.*$", "\\1",
                   df$annotation[genic])
    expect_equal(annGene, as.character(df$geneId[genic]))
})

test_that("annotation and transcript-level columns describe the same transcript", {
    ## issue #252: this peak is located inside an intron of the long BRCA1
    ## isoform while the closest TSS belongs to another isoform.  The reported
    ## gene/transcript columns have to follow the transcript named in the
    ## annotation instead of the nearest one.
    peak <- GRanges("chr1", IRanges(243832958, 243833008))
    pa <- annotatePeak(peak, TxDb = txdb, tssRegion = c(-3000, 3000),
                       level = "transcript", verbose = FALSE)
    df <- as.data.frame(pa)

    ## the annotation text names the transcript that overlaps the peak ...
    expect_match(df$annotation, "^Intron ")
    annotatedTx <- sub("^[A-Za-z']+ \\(([^/]+)/.*$", "\\1", df$annotation)

    ## ... and all other transcript-level fields come from that transcript
    expect_equal(df$transcriptId, annotatedTx)
    expect_equal(df$geneChr, as.character(seqnames(peak)))

    tx <- transcripts(txdb)
    tx <- tx[tx$tx_name == df$transcriptId]
    expect_equal(length(tx), 1L)
    expect_true(overlapsAny(peak, tx))
    expect_equal(df$geneStart, start(tx))
    expect_equal(df$geneEnd, end(tx))
    expect_equal(as.character(df$geneStrand), as.character(strand(tx)))

    ## transcripts() has no gene_id column, build the mapping from
    ## transcriptsBy() whose names are the gene ids
    txg <- transcriptsBy(txdb, by = "gene")
    geneOfTx <- rep(names(txg), elementNROWS(txg))
    names(geneOfTx) <- unlist(txg)$tx_name
    expect_equal(as.character(df$geneId),
                 unname(geneOfTx[df$transcriptId]))

    ## the gene id of the annotation refers to the same gene
    annotatedGene <- sub("^[A-Za-z']+ \\([^/]+/([^,]+),.*$", "\\1", df$annotation)
    expect_equal(as.character(df$geneId), annotatedGene)

    ## distanceToTSS is the strand aware distance to the TSS of that transcript
    strandTx <- as.character(strand(tx))
    tss <- ifelse(strandTx == "+", start(tx), end(tx))
    dStart <- ifelse(strandTx == "+", start(peak) - tss, tss - start(peak))
    dEnd <- ifelse(strandTx == "+", end(peak) - tss, tss - end(peak))
    expect_equal(df$distanceToTSS, ifelse(abs(dStart) <= abs(dEnd),
                                          dStart, dEnd))
})

test_that("transcript-level annotation is self-consistent for many peaks", {
    ## issue #252: geneChr/geneStart/geneEnd/geneId/transcriptId/distanceToTSS
    ## must all be derived from the transcript that the annotation refers to
    tx <- transcripts(txdb)
    tssPeak <- tx[tx$tx_name == "ENST00000336199.9_6"]
    tssPeak <- GRanges(seqnames(tssPeak),
                       IRanges(start(tssPeak), width = 200))

    peaks <- GRanges("chr1",
                     IRanges(start = seq(243660000, 244000000, length.out = 24),
                             width = 200))
    peaks <- c(peaks, tssPeak)

    pa <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000),
                       level = "transcript", verbose = FALSE)
    df <- as.data.frame(pa)
    expect_equal(nrow(df), length(peaks))

    m <- match(df$transcriptId, tx$tx_name)
    expect_false(any(is.na(m)))

    expect_equal(df$geneChr, as.character(seqnames(tx))[m])
    expect_equal(df$geneStart, start(tx)[m])
    expect_equal(df$geneEnd, end(tx)[m])
    expect_equal(as.character(df$geneStrand), as.character(strand(tx))[m])

    ## transcripts() has no gene_id column, build the mapping from
    ## transcriptsBy() whose names are the gene ids
    txg <- transcriptsBy(txdb, by = "gene")
    geneOfTx <- rep(names(txg), elementNROWS(txg))
    names(geneOfTx) <- unlist(txg)$tx_name
    expect_equal(as.character(df$geneId),
                 unname(geneOfTx[df$transcriptId]))

    ## distances are the strand aware distances to the TSS of reported tx
    strandTx <- as.character(strand(tx))[m]
    tss <- as.numeric(ifelse(strandTx == "+", start(tx)[m], end(tx)[m]))
    dStart <- ifelse(strandTx == "+", start(peaks) - tss, tss - start(peaks))
    dEnd <- ifelse(strandTx == "+", end(peaks) - tss, tss - end(peaks))
    expect_equal(df$distanceToTSS, ifelse(abs(dStart) <= abs(dEnd),
                                          dStart, dEnd))

    ## genic annotations must overlap the reported transcript
    genic <- grepl("^(Exon|Intron|5' UTR|3' UTR)", df$annotation)
    expect_true(any(genic))
    expect_true(all(overlapsAny(peaks[genic], tx[m[genic]])))

    ## promoter annotations are TSS based, the nearest TSS is in tssRegion
    promoter <- grepl("^Promoter", df$annotation)
    expect_true(any(promoter))
    expect_true(all(abs(df$distanceToTSS[promoter]) <= 3000))
})
