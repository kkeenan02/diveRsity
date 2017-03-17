\name{fastDivPart}
\alias{fastDivPart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genetic differentiation statistics and their estimators
}

\description{
\code{fastDivPart} is identical to the \code{divPart} function in regards to what it calculates. The difference with this function is the speed with which it processes pairwise calculations. By using more efficient programming techniques, \code{fastDivPart} can execute commands up to 20X faster than \code{divPart}.
}

\details{
All results will be written to a user defined folder ("working\_directory/outfile"). The format of outputs will vary depending on the availability of the package \code{'xlsx'} in the users local package library. If \code{'xlsx'} is available, results will be written to an Excel workbook. If \code{'xlsx'} is not available, results will be written to \code{.txt} files. 
Multi-locus estimates of Weir and Cockerham's theta are calculated as per Weir and Cockerham, 1984.
}

\usage{

fastDivPart(infile = NULL, outfile = NULL, gp = 3, pairwise = FALSE, 
            fst = FALSE, bs_locus = FALSE, bs_pairwise = FALSE,
            boots = 0, plot = FALSE, para = FALSE)
}
			   
\arguments{

\item{infile}{Specifying the name of the \emph{`genepop'} (Rousset, 2008) file from which the statistics are to be calculated This file can be in either the 3 digit of 2 digit format, and must contain only one whitespace separator (e.g. \dQuote{space} or \dQuote{tab}) between each column including the individual names column. The number of columns must be equal to the number of loci + 1 (the individual names column). If this file is not in the \code{working directory} the file path must be given. The name must be a character string (i.e. enclosed in \dQuote{} or `').}

\item{outfile}{Allows users to specify a prefix for an output folder. Name must a character string enclosed in either \dQuote{} or `'.}

\item{gp}{Specifies the digit format of the \code{infile}. Either 3 (default) or 2.}

\item{pairwise}{A logical argument indicating whether standard pairwise diversity statistics should be calculated and returned as a diagonal matrix.}

\item{fst}{A Logical argument indicating whether Weir & Cockerham's 1984 F-statistics should be calculated. NOTE - Calculating these statistics adds significant time to analysis when carrying out pairwise comparisons.}

\item{bs_locus}{Gives users the option to bootstrap locus statistics. Results will be written to \code{.xlsx} workbook by default if the package \code{'xlsx'} is installed, and to a \code{.html} file if \code{plot=TRUE}. If \code{'xlsx'} is not installed, results will be written to \code{.txt} files.}

\item{bs_pairwise}{Gives users the option to bootstrap statistics across all loci for each pairwise population comparison. Results will be written to a \code{.xlsx} file by default if the package \code{'xlsx'} is installed, and to a \code{.html} file if \code{plot=TRUE}. If \code{'xlsx'} is not installed, results will be written to \code{.txt} files.}

\item{boots}{Determines the number of bootstrap iterations to be carried out. The default value is \code{boots = 0}, this is only valid when all bootstrap options are false. There is no limit on the number of bootstrap iterations, however very large numbers of bootstrap iterations can take some time to run. If time is a limiting factor, the function \code{diffCalc} is up to 10x faster than \code{fastDivPart}, but does not allow users to write results to xlsx workbooks or plot results from the function. If required, these must be done manually.}

\item{plot}{Optional interactive \code{.html} image file of the plotted bootstrap results for loci if \code{bs_locus = TRUE} and pairwise population comparisons if \code{bs_pairwise = TRUE}. The default option is \code{plot = FALSE}.}

\item{para}{A logical input, indicating whether your analysis should be run in parallel mode or sequentially.}
}

\value{

\item{standard}{A matrix containing identical data to the \code{Standard_stats} worksheet in the \code{.xlsx} workbook.}

\item{estimate}{A matrix containing identical data to the \code{Estimated_stats} worksheet in the \code{.xlsx} workbook.}

\item{pairwise}{A group of six matrices containing population pairwise statistics. This object is identical to that written as `pairwise-stats` in the \code{.xlsx} workbook.}

\item{bs_locus}{A list containing six matrices of locus values for \emph{Gst}, \emph{G'st}, \emph{D(Jost)}, \emph{Gst-(est)}, \emph{G'st-(est)}, and \emph{D(Jost)-(est)} along with their respective 95\% confidence interval.}

\item{bs_pairwise}{A list containing three-four matrices (depending on whether Weir & Cockerham's Fst is calculated) of pairwise values for \emph{Gst}, \emph{G'st}, \emph{D(Jost)}, \emph{Gst-(est)}, \emph{G'st-(est)}, and \emph{D(Jost)-(est)} along with their respective 95\% confidence intervals, including bias corrected 95\% confidence intervals.}
}

\references{
Dragulescu, A.D., ``xlsx: Read, write, formal Excel 2007 and Excel 97/2000/xp/2003 files'', R package version 0.4.2, url:http://CRAN.R-project.org/package=xlsx, (2012).

Guile, D.P., Shepherd, L.A., Sucheston, L., Bruno, A.,  and Manly, K.F., ``sendplot: Tool for sending interactive plots with tool-tip content.'', R package version 3.8.10, url: http://CRAN.R-project.org/package=sendplot, (2012).

Hedrick, P., ``A standardized genetic differentiation measure,'' Evolution,
vol. 59, no. 8, pp. 1633-1638, (2005).

Jost, L., ``G ST and its relatives do not measure differentiation,'' Molec-
ular Ecology, vol. 17, no. 18, pp. 4015-4026, (2008).

Manly, F.J., ``Randomization, bootstrap and Monte Carlo methods in biology'', Chapman and Hall, London, 1997.

Nei, M. and Chesser, R., ``Estimation of fixation indices and gene diver-
sities,'' Ann. Hum. Genet, vol. 47, no. Pt 3, pp. 253-259, (1983).

R Development Core Team (2011). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria. ISBN
3-900051-07-0, URL http://www.R-project.org/.

Revolution Analytics (2012). doParallel: Foreach parallel adaptor for the parallel package. R package version 1.0.1. http://CRAN.R-project.org/package=doParallel

Revolution Analytics (2012). foreach: Foreach looping construct for R. R package version 1.4.0. http://CRAN.R-project.org/package=foreach

Rousset, F., ``genepop'007: a complete re-implementation of the genepop
software for Windows and Linux.,'' Molecular ecology resources, vol. 8,
no. 1, pp. 103-6, (2008).

Weir, B.S. & Cockerham, C.C., Estimating F-Statistics, for the Analysis of Population Structure, Evolution, vol. 38, No. 6, pp. 1358-1370 (1984). 
}
\author{Kevin Keenan <kkeenan02@qub.ac.uk>}

\examples{
\dontrun{
# simply use the following format to run the function

test_result <- fastDivPart(infile = 'mydata', outfile = "myresults', 
                           gp = 3, pairwise = TRUE, bs_locus = TRUE, 
                           bs_pairwise = TRUE, boots = 1000, 
                           plot = TRUE)
}
}
