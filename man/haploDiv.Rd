\name{haploDiv}
\alias{haploDiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function, allowing the calculation of Weir & Cockerham's (1984) \eqn{F_{ST}} from haploid genotypes in the genepop
}
\description{
haploDiv allows users to calculate Weir & Cockerham's \eqn{F_{ST}} from haploid genotypes for locus, global and pairwise population levels.
}
\usage{
haploDiv(infile = NULL, outfile = NULL, pairwise = FALSE, 
         boots = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{A genepop file/data frame containing haploid genotypes. This file/data frame should contain locus information in either the two digit or three digit format. The argument can be a character string indicating the name of a file or a data frame in the R workspace (e.g. see \code{data("Test_data")}).}
  \item{outfile}{A character string specifying a prefix to be added to output files. A character string specifying a directory location will result in the output files being written to the specified location, rather than the current working directory. If \code{outfile = NULL}, no results will be written to disk.}
  \item{pairwise}{Specifies whether a population pairwise matrix containing \eqn{F_{ST}} values should be calculated.}
  \item{boots}{Specifies whether bootstrapped 95\% confidence intervals should be calculated for each pairwise estimate of \eqn{F_{ST}}. If \code{bootstraps = 0} and \code{pairwise = TRUE}, only a pairwise matrix of \eqn{F_{ST}} will be returned.}
}
\details{
This function uses the same fundamental algorithms as \code{divPart} and \code{fastDivPart}, the only difference being that if \emph{diploidizes} haploid genotypes before calculating statistics. The diploidization process has the effect of changing a haploid genotype into a homozygous diploid genotype for all individuals.
}
\value{
\item{locus}{A named vector of locus estimates of Weir & Cockerham's \eqn{F_{ST}} across all population samples.}
\item{overall}{A global estimate of Weir & Cockerham's \eqn{F_{ST}}.}
\item{pairwise}{A diagonal matrix containing pairwise estimates of Weir & Cockerham's \eqn{F_{ST}} across all loci. Returned when \code{pairwise = TRUE}.}
\item{bs_pairwise}{A data frame with three data columns containing bootstrapped mean, lower 95\% confidence limit and upper 95\% confidence limit for each population pair (rows). Returned when \code{bootstraps > 0} and \code{pairwise = TRUE}.}
}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}
