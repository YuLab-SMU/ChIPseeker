\name{plotPeakProf}
\alias{plotPeakProf}
\title{plotPeakProf}
\usage{
plotPeakProf(peak, TranscriptDb = NULL, upstream = 1000,
  downstream = 1000, xlab = "Genomic Region (5'->3')",
  ylab = "Read count Per Million mapped reads")
}
\arguments{
  \item{peak}{peak file or GRanges object}

  \item{TranscriptDb}{TranscriptDb object}

  \item{upstream}{upstream position}

  \item{downstream}{downstream position}

  \item{xlab}{xlab}

  \item{ylab}{ylab}
}
\value{
ggplot object
}
\description{
plot the profile of peaks that align to flank sequences of
TSS
}
\author{
G Yu
}

