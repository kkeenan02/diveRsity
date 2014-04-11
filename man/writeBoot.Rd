\name{writeBoot}
\alias{writeBoot}
\title{Generate bootstrap pairwise matrices per locus}
\description{Write all pairwise matrices from a bootstrap run to file}
\usage{
writeBoot(infile = NULL, outfile = NULL, gp = 3, bootstraps = 0, parallel = FALSE)
}
\arguments{
  \item{infile}{Character string of file name}
  \item{outfile}{Character string prefix fr output files}
  \item{gp}{Genepop format}
  \item{bootstraps}{Number of bootstrap iterations}
  \item{parallel}{Logical multicore indication}
}
\author{Kevin Keenan}
\keyword{internal}