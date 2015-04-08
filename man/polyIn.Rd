\name{polyIn}
\alias{polyIn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for calculating informativeness for the inference of ancestry from loci of any ploidy.
}
\description{
This function will calculate Rosenberg et al's, (2003) In for loci of any ploidy level. The statistic can be calculated across all samples or on a pairwise basis. The function is efficient for a large number of loci.
}
\usage{
polyIn(infile = NULL, pairwise = FALSE, para = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{A character string pointing to an input file containing genotypes in a modified genepop format. Rather that genotypes being coded as either two digit or three digit numbers, polyploid genotypes are coded as \code{AABA, AAAA, ... ABAB} etc. for a tetraploid or \code{ABA, AAA, ... ABB} etc. for a triploid. The overall structure of the input file should be that of a genepop file.}
  \item{pairwise}{A logical argument specifying whether In should be calculated for pairwise population comparisons. Large number of population samples will result in long computation times.}
  \item{para}{A logical argument specifying whether multiple CPU cores should be used. If \code{para = TRUE}, all available system CPUs will be used.}
}
\value{
If \code{pairwise = TRUE}, a list of pairwise matrices for each locus is returned. These list elements are named as per loci names in \code{infile} to allow simple indexing for loci of interest. If \code{pairwise = FALSE} and named numeric vector is returned.
}
\references{
Rosenberg, N., Li, L., Ward, R., and Pritchard, J., (2003) 'Informativeness of genetic markers for inference of of ancestry.', American Journal of Human Genetics, vol. 73, no. 6 pp. 1402-22.
}
\author{
Kevin Keenan 2014
}