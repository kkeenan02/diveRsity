\name{divRatio}
\alias{divRatio}

\title{
Calculates the standardised diversity ratios relative to a 'yardstick' reference population following Skrbinsek et al., 2012.}

\description{
Diversity ratios derived from allelic richness and expected heterozygosity are calculated from either a genepop file containing raw data for all populations of interest, or from a genepop file containing raw data for only the reference population and a data frame containing relevant information for the populations of interest. See below for more details the structure of this data frame.
}

\details{

All results will be written to a user defined folder, providing an argument is passed for 'outfile'. Results will be written to \code{.xlsx} files if the package \code{xlsx} and its dependencies are installed, or a \code{.txt} file otherwise.

}

\usage{
divRatio(infile = NULL, outfile = NULL, gp = 3, pop_stats = NULL, 
         refPos = NULL, boots = 1000, para = FALSE)
}

\arguments{

\item{infile}{A character string argument specifying the name of either a 3 digit or 2 digit genepop file containing the raw genotypes of at least the reference population sample.}

\item{outfile}{A character string specifying a prefix name for an automatically generated results folder, to which results file will be written.}

\item{gp}{Specifies the digit format of the \code{infile}. Either 3 (default) or 2.}

\item{pop_stats}{A character string indicating the name of the population statistics data frame file. This argument is required if only raw data for the reference population are give in \code{infile}. The data frame should be structured in a specific way. An example can be seen by typing \code{data(pop_stats)} into the console. The validloci column is only required if mean allelic richness and expected heterozygosity for populations of interest have been calculated from loci for which data is not present in the reference population. This column should contain a single character string of common loci between each population sample and the reference population sample.}

\item{refPos}{A numeric argument specifying the position of the reference population in \code{infile}. The argument is only valid when raw genotype data has been provided for the reference population sample and all other populations of interest and \code{pop_stats} is \code{NULL}.}

\item{boots}{Specifies the number of times the reference population should be resampled when calculating the sample size standardised allelic richness and expected heterozygosity for calculating the diversity ratios. The larger the number of bootstraps the longer the analysis will take to run. As an indication of runtime, running \code{divRatio} on the \code{Big_data} data set (type \code{?Big_data} for details), takes 10min 42s on a Toshiba Satellite R830 with 6GB RAM, and an Intel Core i5 - 2435M CPU running Linux.}

\item{para}{A logical argument indicating whether the analysis should make use of all available cores on the users system.}

}

\value{

A data frame containing the following columns:

\item{pop}{The names of each population of interest, including the reference population}
\item{n}{The sample size of each population}
\item{alr}{Mean allelic richness across loci}
\item{alrSE}{The standard error of the allelic richness across loci}
\item{He}{Mean expected heterozygosity across loci}
\item{HeSE}{Standard error of expected heterozygosity across loci}
\item{alrRatio}{The ratio of the allelic richness of the subject population sample and the sample size standardised reference population allelic richness}
\item{alrSEratio}{The standard error of divisions for the allelic richness ratio}
\item{heRatio}{The ratio of expected heterozygosity between the standardised reference population sample and subject population samples}
\item{heSEratio}{The standard error of divisions for the expected heterozygosity ratio}

}

\references{

Skrbinsek, T., Jelencic, M., Waits, L. P., Potocnik, H., Kos, I., & Trontelj, P. (2012). Using a reference population yardstick to calibrate and compare genetic diversity reported in different studies: an example from the brown bear. Heredity, 109(5), 299-305. doi:10.1038/hdy.2012.42

}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}

\examples{
\dontrun{
# To run an example use the following format

test_results <- divBasic(infile = Test_data, outfile = 'out', 
                         gp = 3, pop_stats = NULL, refPos = NULL, 
                         bootstraps = 1000, parallel = TRUE)
}
}
