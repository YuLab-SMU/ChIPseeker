\name{vennplot}
\alias{vennplot}
\title{vennplot}
\usage{
vennplot(Sets, by = "gplots")
}
\arguments{
  \item{Sets}{a list of object, can be vector or GRanges
  object}

  \item{by}{one of gplots or Vennerable}
}
\value{
venn plot that summarize the overlap of peaks from
different experiments or gene annotation from different
peak files.
}
\description{
plot the overlap of a list of object
}
\examples{
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
peakAnno <- annotatePeak(peakfile, TranscriptDb=txdb)
peakAnnoList <- lapply(1:3, function(i) peakAnno[sample(1:length(peakAnno), 100),])
names(peakAnnoList) <- paste("peak", 1:3, sep="_")
genes= lapply(peakAnnoList, function(i) unlist(i$geneId))
vennplot(genes)
}
\author{
G Yu
}

