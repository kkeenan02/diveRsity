\name{divSimCo}
\alias{divSimCo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Similarity coefficients for co-dominant diploid genenotype data
}
\description{
Calculation of similarity coefficients for co-dominant genotype data, following Kosman and Leonard, (2005).
}
\usage{
divSimCo(infile = NULL, loci = FALSE, boots = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{Specifying the name of the \emph{`genepop'} (Rousset, 2008) file from which the statistics are to be calculated This file can be in either the 3 digit of 2 digit format. See \url{http://genepop.curtin.edu.au/help_input.html} for detail on the genepop file format.}
  
\item{loci}{A logical argument, specifying whether similarity matrices for each locus should be calculated. If \code{FALSE}, only a global similarity matrix is given.}

\item{boots}{An integer, specifying the number of bootstrapped matrices to generate. Loci are resampled with replacement.}

}
\details{
Similarity coefficients are calculated according to Kosman and Leonard, (2005).
Multi-locus values are calculated from all loci without missing data.
}
\value{
A list containing either one (\code{global}) if \code{loci = FALSE}, or two elements (\code{global} and \code{loci}) if \code{loci = TRUE}.
  \item{global}{A pairwise matrix of multilocus similarity coefficients.}
  \item{loci}{A list of pairwise matrices of locus similarity coefficients.}
}
\references{
Kosman E, Leonard K (2005) Similarity coefficients for molecular markers in studies of genetic relationships between individuals for haploid, diploid, and polyploid species. Molecular ecology, 14, 415-424.

}
\author{
Kevin Keenan
}