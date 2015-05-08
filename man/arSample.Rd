\name{arSample}
\alias{arSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resample based calculation of allelic richness
}
\description{
Resample based calculation of allelic richness
}
\usage{
arSample(dat, nrep, ci = FALSE, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{List output from diveRsity:::rgp}
  \item{nrep}{Number of bootstrap replicates}
  \item{ci}{A logical indicating if a confidence interval should be calculated}
  \item{alpha}{A numeric argument specifying the alpha/2 and 1-(alpha/2) limits for confidence intervals.}
}
\keyword{internal}