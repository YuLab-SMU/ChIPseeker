\name{plotAnnoPie}
\alias{plotAnnoPie}
\title{plotAnnoPie}
\usage{
plotAnnoPie(peakAnno, ndigit = 2, radius = 0.8, explode = 0.05,
  labelcex = 1, ...)
}
\arguments{
  \item{peakAnno}{peakAnno}

  \item{ndigit}{ndigit}

  \item{radius}{the radius of the pie}

  \item{explode}{the amount to "explode" the pie}

  \item{labelcex}{label font size}

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

