\name{getPromoters}
\alias{getPromoters}
\title{getPromoters}
\usage{
getPromoters(TranscriptDb = NULL, upstream = 1000, downstream = 1000)
}
\arguments{
  \item{TranscriptDb}{TranscriptDb}

  \item{upstream}{upstream position}

  \item{downstream}{downstream position}
}
\value{
GRanges object
}
\description{
extract promoter positions from TranscriptDb.
}
\author{
G Yu
}

