\name{basicStats}
\alias{basicStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Claculate basic descriptive population parameters from a genepop genotype file
}

\description{
A new and improved version of \code{divBasic}, with a more comprehensive set of control arguments, as well as better memory and CPU efficiency. Additional options for testing homozygosity/heterozygosity excess are also provided.
}

\details{

HWE tests are carried out using methods originally implemented in the package \code{HWxtest}. Test of HWE proportions are either carried out using either exact testing (when the number of genotype tables is less than 10 million), or through Monte Carlo simulations, where the number of random trials can be specified using the argument \code{mc_reps}. HWE tests will also be automatically carried out to assess homozygosity/hetrozygosity excess (i.e. one tailed tests). The test results returned is automatically determined by the proporties of the original data. See \code{?HWx.test} for more details.

}

\usage{
basicStats(infile = NULL, outfile = NULL, fis_ci = FALSE,
           ar_ci = FALSE, fis_boots = NULL, ar_boots = NULL, 
           mc_reps = 9999, rarefaction = TRUE, ar_alpha = 0.05, 
           fis_alpha = 0.05)
}
\arguments{

\item{infile}{Specifying the name of the \emph{`genepop'}(Rousset, 2008) file from which the statistics are to be calculated. This file can be in either the 3 digit of 2 digit format, and must contain only one \emph{whitespace} separator (e.g. \dQuote{space} or \dQuote{tab}) between each column including the individual names column. The number of columns must be equal to the number of loci + 1 (the individual names column). If this file is not in the \code{working directory} the file path must be given. The name must be a character string (i.e. enclosed in \dQuote{} or `').}

\item{outfile}{A character string specifying the name to be given to the output tab delimited file. Default value of \code{NULL} means that no file is written.}

\item{fis_ci}{A logical argument specifying whether confidence limits should be calculated for the inbreeding coefficient.}

\item{ar_ci}{A logical argument specifying whether confidence limits should be calculated for resampleed allelic richness. Only valid when \code{rarefaction = FALSE}.}

\item{fis_boots}{An integer specifying the number of bootstrap replicate used to calculate confidence limits for the inbreeding coefficient if \code{fis_ci = TRUE}}

\item{ar_boots}{An integer specifying the number of bootstrap replicate used to calculate confidence limits for allelic richness if \code{ar_ci = TRUE} and \code{rarefaction = FALSE}}

\item{mc_reps}{An integer specifying the number of Monte Carlo replicates used when carrying out psudo-exact HWE tests}

\item{rarefaction}{A logical indicating whether allelic richness should be calculated using rarefaction. If \code{FALSE} user should specify the number of \code{ar_boots} to use when estimating allelic richness using the resample method}

\item{ar_alpha}{A numeric argument specifying the alpha used to calculate the confidence interval for resample based allelic richness. Valid if \code{rarefaction = FALSE}, \code{ar_boots > 2} and \code{ar_ci = TRUE}. Confidence limits are defined as alpha/2 and 1-(alpha/2).}

\item{fis_alpha}{A numeric argument specifying the alpha used to calculate the confidence interval for resample based allelic richness. Valid if \code{fis_boots > 2} and \code{fis_ci = TRUE}. Confidence limits are defined as alpha/2 and 1-(alpha/2).}

}

\value{
A list of results
}

\references{
Rousset, F., ``genepop'007: a complete re-implementation of the genepop
software for Windows and Linux.,'' Molecular ecology resources, vol. 8,
no. 1, pp. 103-6, (2008).

Bill Engels (2014). HWxtest: Exact Tests for Hardy-Weinberg Proportions. R
package version 1.0.5. http://CRAN.R-project.org/package=HWxtest

}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}

\examples{
\dontrun{
# To run an example use the following format
library("diveRsity")
data(Test_data)
test_results <- basicStats(infile = Test_data, outfile = "test_run", 
                           fis_ci = TRUE, ar_ci = TRUE, fis_boots = 999, 
                           ar_boots = 999, mc_reps = 9999, 
                           rarefaction = FALSE, ar_alpha = 0.05, 
                           fis_alpha = 0.05)
}
}