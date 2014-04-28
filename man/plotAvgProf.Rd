\name{plotAvgProf}
\alias{plotAvgProf}
\title{plotAvgProf}
\usage{
plotAvgProf(peak, weightCol = NULL, TranscriptDb = NULL, upstream = 1000,
  downstream = 1000, xlab = "Genomic Region (5'->3')",
  ylab = "Read Count Frequency", verbose = TRUE)
}
\arguments{
  \item{peak}{peak file or GRanges object}

  \item{weightCol}{column name of weight}

  \item{TranscriptDb}{TranscriptDb object}

  \item{upstream}{upstream position}

  \item{downstream}{downstream position}

  \item{xlab}{xlab}

  \item{ylab}{ylab}

  \item{verbose}{print message or not}
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

