\name{inCalc}
\alias{inCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate locus informative for the inference of ancestry
}

\description{
\code{inCalc} allows the calculation of locus informativeness for ancestry (\emph{In}), (Rosenberg \emph{et al.,} 2003), both across all population samples under consideration and for all pairwise combinations of population samples. These data can be bootstrapped using the same procedure as above, to obtain 95\% confidence intervals.
}

\details{
All results will be written to a user defined folder ("working_directory/outfile"). The format of outputs will vary depending on the availability of the package \code{'xlsx'} in the users local package library. If \code{'xlsx'} is available, results will be written to an Excel workbook. If \code{'xlsx'} is not available, results will be written to \code{.txt} files.
}

\usage{
inCalc(infile, outfile = NULL, gp = 3, 
       bs_locus = FALSE, bs_pairwise = FALSE, 
       bootstraps = 0, plot = FALSE, parallel = FALSE)
}
\arguments{

\item{infile}{Specifying the name of the \emph{`genepop'}(Rousset, 2008) file from which the statistics are to be calculated This file can be in either the 3 digit of 2 digit format, and must contain only one \emph{whitespace} separator (e.g. \dQuote{space} or \dQuote{tab}) between each column including the individual names column. The number of columns must be equal to the number of loci + 1 (the individual names column). If this file is not in the \code{working directory} the file path must be given. The name must be a character string (i.e. enclosed in \dQuote{} or `').}

\item{outfile}{Allows users to specify a prefix for an output folder. Name must a character string enclosed in either ``'' or `'.}

\item{gp}{Specifies the digit format of the \code{infile}. Either 3 (default) or 2.}

\item{bs_locus}{Gives users the option to bootstrap locus statistics. Results will be written to \code{.xlsx} file by default if the package \code{'xlsx'} is installed, and to a \code{.png} file if \code{plot=TRUE}. If \code{'xlsx'} is not installed, results will be written to \code{.txt} files.}

\item{bs_pairwise}{Gives users the option to bootstrap statistics across all loci for each pairwise population comparison. Results will be written to a \code{.xlsx} file by default if the package \code{'xlsx'} is installed. If \code{'xlsx'} is not installed, results will be written to \code{.txt} files.}

\item{bootstraps}{Determines the number of bootstrap iterations to be carried out. The default value is \code{bootstraps = 0}, this is only valid when 
all bootstrap options are false.  Again, there is no limit on the number of bootstrap iterations, however very large numbers of bootstrap iterations (i.e. > 1000) will take quite a long time to run.}

\item{plot}{Optional \code{.png} image file of the plotted bootstrap results for loci if \code{bs_locus = TRUE}. The default option is \code{plot = FALSE}.}

\item{parallel}{Allows for parallel computation of pairwise locus \emph{In}. The number of available core is automatically detected if \code{parallel = TRUE}.}
}

\value{

\item{Allele_In}{A character matrix of allelic \emph{In} values per locus along with locus sums.}

\item{l_bootstrap}{A character matrix of locus \emph{In} values as well as 95\% confidence intervals, calculated from bootstraps (Manly, 1997). Returned when \code{bs_locus} is set to \code{TRUE}.}

\item{PW_bootstrap}{A list of matrices of bootstrapped pairwise locus \emph{In} results corresponding to each pairwise comparison.}
}

\references{
Dragulescu, A.D., ``xlsx: Read, write, formal Excel 2007 and Excel 97/2000/xp/2003 files'', R package version 0.4.2, url:http://CRAN.R-project.org/package=xlsx, (2012).

Manly, F.J., ``Randomization, bootstrap and Monte Carlo methods in biology'', Chapman and Hall, London, 1997.

Rosenberg, N., Li, L., Ward, R., and Pritchard, J., ``Informativeness of
genetic markers for inference of ancestry.,'' American Journal of Human
Genetics, vol. 73, no. 6, pp. 1402-22, (2003).

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
# test_results<-inCalc(infile = Test_data, outfile = 'out',
#                      gp = 3, bs_locus = TRUE, bs_pairwise = TRUE,
#                      bootstraps = 10, plot = FALSE, parallel = TRUE)
}
