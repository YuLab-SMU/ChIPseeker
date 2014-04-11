\name{getGenomicAnnotation}
\alias{getGenomicAnnotation}
\title{getGenomicAnnotation}
\usage{
getGenomicAnnotation(peaks, distance, tssRegion = c(-3000, 100), TranscriptDb)
}
\arguments{
  \item{peaks}{peaks in GRanges object}

  \item{distance}{distance of peak to TSS}

  \item{tssRegion}{tssRegion, default is -3kb to +100bp}

  \item{TranscriptDb}{TranscriptDb object}
}
\value{
character vector
}
\description{
get Genomic Annotation of peaks
}
\author{
G Yu
}

