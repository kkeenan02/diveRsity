\name{fstOnly}
\alias{fstOnly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A minimal function for the calculate of Weir & Cockerham's (1984) F_{ST} and F_{IT} from codominant molecular data in the genepop format.
}

\description{
This function calculates locus and pairwise confidence intervals for Weir & Cockerham's (1984) F_{ST} and F_{IT}. These statistics can also be calculated using the \code{divPart} function, however, \code{fstOnly} is designed to be more memory efficient for larger datasets (e.g. SNPs).
}

\details{
Because \code{fstOnly} is intended to maximise memory (RAM) efficiency, the function does not provide many of the plotting utilities that \code{divPart} does.
}

\usage{
fstOnly(infile = NULL, outfile = NULL, gp = 3, bs_locus = FALSE,
        bs_pairwise = FALSE, bootstraps = 0, parallel = FALSE)
}

\arguments{

\item{infile}{Specifying the name of the \emph{`genepop'} (Rousset, 2008) file from which the statistics are to be calculated This file can be in either the 3 digit of 2 digit format, and must contain only one whitespace separator (e.g. \dQuote{space} or \dQuote{tab}) between each column including the individual names column. The number of columns must be equal to the number of loci + 1 (the individual names column). If this file is not in the \code{working directory} the file path must be given. The name must be a character string (i.e. enclosed in \dQuote{} or `').}

\item{outfile}{Allows users to specify a prefix for an output folder. Name must a character string enclosed in either \dQuote{} or `'.}

\item{gp}{Specifies the digit format of the \code{infile}. Either 3 (default) or 2.}

\item{bs_locus}{Specifies whether locus bootstrapped confidence intervals should be calculated.}

\item{bs_pairwise}{Specified whether pairwise population bootstrapped confidence intervals should be calculated.}

\item{bootstraps}{Determines the number of bootstrap iterations to be carried out. The default value is \code{bootstraps = 0}, this is only valid when all bootstrap options are false. There is no limit on the number of bootstrap iterations, however very large numbers of bootstrap iterations (< 1000) on even modest data sets (e.g. 265 individuals x 38 loci) will take over 20 minutes to run on a most PCs).}

\item{parallel}{A logical input, indicating whether your analysis should be run in parallel mode or sequentially. \code{parallel} = \code{TRUE} is only valid if the packages, \code{parallel}, \code{doParallel} and \code{foreach} are installed.}

}

\value{

\item{locus}{A list object containing two matrices, F_{ST} and F_{IT}. These matrices contain actual, lower 95\% confidence interval and upper 95\% confidence interval per locus. Global values are also presented with their respective confidence intervals.}

\item{pairwise}{A list object containing two matrices, F_{ST} and F_{IT}. These matrices contain actual, lower 95\% confidence interval and upper 95\% confidence interval per pairwise population comparison.}

}

\note{
This function has become obsolete following improvements to \code{fastDivPart} and \code{diffCalc}. The use of either of these two functions is recommended over \code{fstOnly}. \code{fstOnly} will be deprecated in future \code{diveRsity} releases.
}

\references{
Rousset, F., ``genepop'007: a complete re-implementation of the genepop
software for Windows and Linux.,'' Molecular ecology resources, vol. 8,
no. 1, pp. 103-6, (2008).

Weir, B.S. & Cockerham, C.C., Estimating F-Statistics, for the Analysis of Population Structure, Evolution, vol. 38, No. 6, pp. 1358-1370 (1984). 

}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}
