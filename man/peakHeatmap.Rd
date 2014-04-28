\name{peakHeatmap}
\alias{peakHeatmap}
\title{peakHeatmap}
\usage{
peakHeatmap(peak, weightCol = NULL, TranscriptDb = NULL, upstream = 1000,
  downstream = 1000, xlab = "", ylab = "", title = NULL, color = NULL,
  verbose = TRUE)
}
\arguments{
  \item{peak}{peak file or GRanges object}

  \item{weightCol}{column name of weight}

  \item{TranscriptDb}{TranscriptDb object}

  \item{upstream}{upstream position}

  \item{downstream}{downstream position}

  \item{xlab}{xlab}

  \item{ylab}{ylab}

  \item{title}{title}

  \item{color}{color}

  \item{verbose}{print message or not}
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

