\name{bigPreDiv}
\alias{bigPreDiv}
\title{Memory efficient diversity statistic calculations}
\description{Internal function for \code{bigDivPart}}
\usage{
bigPreDiv(prePopList, bs = FALSE, nloci, npops, popSizes, fstat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prePopList}{A custom charcter vector to be parsed}
  \item{bs}{Logical resample indicator}
  \item{nloci}{Integer. Number of loci}
  \item{npops}{Integer. Number of population samples}
  \item{popSizes}{Integer vector. Number of individuals per sample}
  \item{fstat}{Logical WC calculation indicator}
}
\author{Kevin Keenan}
\keyword{internal}