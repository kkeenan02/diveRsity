\name{statCalc}
\alias{statCalc}
\title{Internal differentiation calculator}
\description{Internal differentiation calculator}
\usage{
statCalc(rsDat, idx = NULL, al, fst, bs = TRUE, ci_type = "individuals")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{Character string pointing to input file}
  \item{idx}{an index list for resamples}
  \item{al}{A list of allele frequency matrices}
  \item{fst}{Logical indication of whether WC f-stats should be calculated}
  \item{bs}{Logical indiction of whether \code{rsDat} should be bootstrapped}
  \item{ci_type}{A character string indicating whether bootstrapping should be
  carried out over individuals within samples (\code{``individuals''}.), or
  across loci (\code{``loci''}).}
}
\author{Kevin Keenan}
\keyword{internal}