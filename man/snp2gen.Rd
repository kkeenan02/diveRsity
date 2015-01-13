\name{snp2gen}
\alias{snp2gen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion function for SNP nucleotide genotype matrix to a genepop file or genepop data object.
}
\description{
This function converts SNP nucleotide genotype to genepop file format. The genepop object can be written to file using the \code{write} argument.
}
\usage{
snp2gen(infile = NULL, prefix_length = 2, write = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{A character string indicating the name of the text file containing SNP genotypes.}
  \item{prefix_length}{This argument specifies the population specific prefix within individual names.}
  \item{write}{A logical indication of whether the genepop object should be written to file.}
}
\value{
A genepop file
}
\author{
Kevin Keenan, 2014
}


\examples{
\dontrun{
data(SNPs, package = "diveRsity")
snp2gen(infile = SNPs, prefix_length = 2, write = FALSE)
}
}
