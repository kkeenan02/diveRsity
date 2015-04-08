\name{divBasic}
\alias{divBasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate basic population parameters such as allelic richness, observed heterozygosity, as well as expected heterozygosity.
}

\description{
\code{divBasic} allows the calculation of locus and overall basic population parameters. \code{divBasic} will write results to a \emph{.xlsx} workbook. The function accepts co-dominant genetic data in both 2 and 3 digit genepop formats.
}

\details{
All results will be written to a user defined folder ("working_directory/outfile"), providing an argument is passed for 'outfile'. Results will be written to \code{.xlsx} files, and multiple \code{R} objects are also written to the current environment.

HWE tests can be carried out using either a standard Chisq goodness of fit method, or using Fisher's exact method. The standard chisq test behave poorly when there are classes with low numbers of observations (e.g. hypervariable microsatellite loci). In such instance it is advisable to use exact testing. Multi-locus HWE is tested using the standard chisq method by summing chisq difference and degrees of freedom across loci, and using these parameter to derive a pvalue for the test. When using exact testing, the multi-locus pvalue is determined using Fisher's method for combining pvalue from independent tests. This process assumes that loci are unlinked.
}

\usage{
divBasic(infile = NULL, outfile = NULL, gp = 3, bootstraps = NULL,
         HWEexact = FALSE, mcRep = 2000)
}
\arguments{

\item{infile}{Specifying the name of the \emph{`genepop'}(Rousset, 2008) file from which the statistics are to be calculated. This file can be in either the 3 digit of 2 digit format, and must contain only one \emph{whitespace} separator (e.g. \dQuote{space} or \dQuote{tab}) between each column including the individual names column. The number of columns must be equal to the number of loci + 1 (the individual names column). If this file is not in the \code{working directory} the file path must be given. The name must be a character string (i.e. enclosed in \dQuote{} or `').}

\item{outfile}{Allows users to specify a prefix for an output folder. Name must a character string enclosed in either ``'' or `'.}

\item{gp}{Specifies the digit format of the \code{infile}. Either 3 (default) or 2.}

\item{bootstraps}{This argument specifies how many bootstrap iterations should be executed when calculating 95\% confidence intervals for $F_is$. The argument should be an integer greater than 1. Setting \code{bootstrap = NULL} suppresses the calculation of F_is. Users should note that setting this argument to values greater than 1000 may result in longer executions times.}

\item{HWEexact}{A logical argument specifying if HWE testing should be carried out using Fisher's exact tests.}

\item{mcRep}{An integer specifying the number of replicates to use for the Monte Carlo tests if HWEexact is TRUE.}

}

\value{

\item{locus_pop_size}{A matrix containing the number of individuals typed per locus per population sample. Mean values across loci are also given.}

\item{Allele_number}{A matrix containing the number of alleles observed per locus per population sample. Mean values across loci are also given.}

\item{proportion_Alleles}{A matrix containing the percentage of total alleles observed per locus per population sample. Mean values across loci are also given.}

\item{Allelic_richness}{A matrix containing the allelic richness per locus per population sample. Allelic richness is calculated using 1000 re-samples (n = smallest sample in the input data file), with replacement per population sample locus per population sample. Mean values across loci are also given.}

\item{Ho}{A matrix containing observed heterozygosity per locus per population sample. Mean values across loci are also given.}

\item{He}{A matrix containing expected heterozygosity per locus per population sample. Mean values across loci are also given.}

\item{HWE}{A matrix containing uncorrected \emph{p}-values from chi-square test for goodness-of-fit to Hardy-Weinberg equilibrium. Overall \emph{p}-values are also given per population sample.}

\item{fis}{A list of dataframes containing locus and global F_is values for each population sample. In each dataframe the actual F_is is listed in the first column of the matrix, lower and upper 95\% confidence intervals are listed in the next two columns, while bias corrected 95\% CI are listed in the last two columns. This object is only returned when \code{bootstraps} is not \code{NULL}}

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
\dontrun{
# To run an example use the following format

test_results <- divBasic(infile = Test_data, outfile = 'out', gp = 3, bootstraps = 1000)
}
}
