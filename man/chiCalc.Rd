\name{chiCalc}
\alias{chiCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A functions for the assessment of population genetic heterogeneity
}

\description{
\code{chiCalc} allows the calculation of locus $X^{2}$ as well as $X^{2}$ across loci and populations. P.values and significance indicators are also provided. The function accepts co-dominant genetic data in both 2 and 3 digit genepop formats.
}

\details{
All results will be written to a user defined folder ("working_directory/outfile"), providing an argument is passed for 'outfile'. Results will be written to \code{.txt} files, and an \code{R}, object is also written to the current environment.
}

\usage{
chiCalc(infile = NULL, outfile = NULL, gp = 3, minFreq = NULL)
}
\arguments{

\item{infile}{Specifying the name of the \emph{`genepop'}(Rousset, 2008) file from which the statistics are to be calculated. This file can be in either the 3 digit of 2 digit format, and must contain only one \emph{whitespace} separator (e.g. \dQuote{space} or \dQuote{tab}) between each column including the individual names column. The number of columns must be equal to the number of loci + 1 (the individual names column). If this file is not in the \code{working directory} the file path must be given. The name must be a character string (i.e. enclosed in \dQuote{} or `').}

\item{outfile}{Allows users to specify a prefix for an output folder. Name must a character string enclosed in either ``'' or `'.}

\item{gp}{Specifies the digit format of the \code{infile}. Either 3 (default) or 2.}

\item{minFreq}{This argument allows users to specify a minimum mean allele frequency, below which, allele classes should be excluded from the analysis. The chi-square test for heterogeneity works well for classes which are observed at least five times. As many alleles at microsatellite loci are not observed this frequently, this option allow for a more conservative assessment of population heterogeneity, only considering alleles equal to or above the specified frequency. minFreq can be single nominal frequency or a vector of frequencies.}
}

\value{

\item{chiTable}{A character matrix of locus and overall $X^{2}$ values, as well as respective degrees of freedom, p.values and visual statistical significance.}

}

\references{
Rousset, F., ``genepop'007: a complete re-implementation of the genepop
software for Windows and Linux.,'' Molecular ecology resources, vol. 8,
no. 1, pp. 103-6, (2008).
}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}

\examples{
# To run an example use the following format
#
# test_results <- chiCalc(infile = Test_data, outfile = 'out',
#                       gp = 3, minFreq = 0.05)
}
