\name{bigDivPart}
\alias{bigDivPart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genetic differentiation statistics and their estimators for high throughput data (e.g. RAD-seq derived SNPS)
}

\description{
\code{bigDivPart} allows for the calculation of four main diversity partition statistics and their respective estimators from large genepop files.
}

\details{
All results are optionally written to a user defined directory in either .xlsx or .txt format. The function is identical to the \code{divPart} function without the pairwise or bootstrapping functionality.
}

\usage{

bigDivPart(infile = NULL, outfile = NULL,
           WC_Fst = FALSE, format = NULL)
}
			   
\arguments{

\item{infile}{Specifying the name of the \emph{`genepop'} (Rousset, 2008) file from which the statistics are to be calculated This file can be in either the 3 digit of 2 digit format, and must contain only one whitespace separator (e.g. \dQuote{space} or \dQuote{tab}) between each column including the individual names column. The number of columns must be equal to the number of loci + 1 (the individual names column). If this file is not in the \code{working directory} the file path must be given. The name must be a character string (i.e. enclosed in \dQuote{} or `').}

\item{outfile}{Allows users to specify a prefix for an output folder. Name must a character string enclosed in either \dQuote{} or `'.}

\item{WC_Fst}{A Logical argument indicating whether Weir & Cockerham's 1984 F-statistics should be calculated.}

\item{format}{A character string specifying the preferred output format for calculated results. The arguments \code{txt} or \code{xlsx} are valid when \code{outfile} is not \code{NULL}.}
}

\value{

\item{standard}{A matrix containing identical data to the \code{Standard_stats} worksheet in the \code{.xlsx} workbook.}

\item{estimate}{A matrix containing identical data to the \code{Estimated_stats} worksheet in the \code{.xlsx} workbook.}

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
 
test_result <- bigDivPart(infile = 'mydata', outfile = "myresults',
                          WC_Fst = TRUE, format = "txt")
}
}
