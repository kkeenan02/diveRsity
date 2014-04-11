\name{bsFun}
\alias{bsFun}
\title{
divMigrate bootstrap routine
}
\description{
divMigrate bootstrap routine
}
\usage{
bsFun(genos, idx, af, pw)
}

\arguments{
  \item{genos}{A list of arrays containing population genotypes}
  \item{idx}{A list of bootstrap indeces}
  \item{af}{A list of allele frequency matrices}
  \item{pw}{Pairwise combination integer matrix}
}
\value{
  \item{dRel}{A matrix containing relative migration values}
}
\author{
Kevin Keenan
}
\keyword{internal}