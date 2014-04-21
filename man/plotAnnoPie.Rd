\name{plotAnnoPie}
\alias{plotAnnoPie}
\title{plotAnnoPie}
\usage{
plotAnnoPie(peakAnno, ndigit = 2, cex = 0.9, col = NA, pie3D = FALSE,
  ...)
}
\arguments{
  \item{peakAnno}{peakAnno}

  \item{ndigit}{number of digit to round}

  \item{cex}{label cex}

  \item{col}{color}

  \item{pie3D}{plot in 3D or not}

  \item{...}{extra parameter}
}
\value{
pie plot of peak genomic feature annotation
}
\description{
pieplot from peak genomic annotation
}
\examples{
## example not run
## require(TxDb.Hsapiens.UCSC.hg19.knownGene)
## txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
## peakfile <- system.file("extdata", "sample_peaks.txt", package="chipseeker")
## peakAnno <- annotatePeak(peakfile, TranscriptDb=txdb)
## plotAnnoPie(peakAnno)
}
\author{
G Yu
}
\seealso{
\code{\link{annotatePeak}} \code{\link{plotAnnoBar}}
}

