\name{genos2mat}
\alias{genos2mat}

\title{
Convert genos loci into presance absence matrices
}
\description{
A matrix for each locus is constructed with n rows and c columns, where n is the number of genotype individuals, and c is the total number of observed alleles at a locus. Individual are surveyed for alleles and assigned 0 for an absence, 0.5 for a single presence (i.e. heterozygote), and 1 for double presence (i.e. homozygote).
}
\usage{
genos2mat(mat, ip, na)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{A template matrix}
  \item{ip}{A character vector of alleles (e.g. \code{genos[,i,j]}) for all genotyped individual.}
  \item{na}{An NA replacement template}
}
\author{
Kevin Keenan
}
\keyword{internal}
