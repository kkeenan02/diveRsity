\name{snp2gp}
\alias{snp2gp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion function for SNP nucleotide genotype matrix to a genepop file.
}
\description{
This function converts SNP nucleotide genotype to genepop file format.
}
\usage{
snp2gp(infile, prefix_length)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{A character string indicating the name of the text file containing SNP genotypes.}
  \item{prefix_length}{This argument specifies the population specific prefix affixed to individual names.}
}
\value{
A genepop file
}
\author{
Kevin Keenan, 2014
}


\examples{
## NOT RUN
# data(SNPs, package = "diveRsity")
# snp2gp(infile = SNPs, prefix_length = 2)
}