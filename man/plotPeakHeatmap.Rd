\name{plotPeakHeatmap}
\alias{plotPeakHeatmap}
\title{plotPeakHeatmap}
\usage{
plotPeakHeatmap(peak, TranscriptDb = NULL, upstream = 1000,
  downstream = 1000, color = "red")
}
\arguments{
  \item{peak}{peak file or GRanges object}

  \item{TranscriptDb}{TranscriptDb object}

  \item{upstream}{upstream position}

  \item{downstream}{downstream position}

  \item{color}{color}
}
\value{
figure
}
\description{
plot the heatmap of peaks align to flank sequences of TSS
}
\author{
G Yu
}

