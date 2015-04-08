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
All results will be written to a user defined folder ("working\_directory/outfile"). The format of outputs will vary depending value of the \code{xlsx} argument. If \code{xlsx = TRUE}, results will be written to a .xlsx workbook using the \code{xlsx} package. If \code{xlsx = FALSE}, results will be written to \code{.txt} files.
}

\usage{
inCalc(infile = NULL, outfile = NULL, pairwise = FALSE, xlsx = FALSE, 
       boots = NULL, para = FALSE)
}
\arguments{

\item{infile}{Specifying the name of the \emph{`genepop'} (Rousset, 2008) file from which the statistics are to be calculated This file can be in either the 3 digit of 2 digit format. See \url{http://genepop.curtin.edu.au/help_input.html} for detail on the genepop file format.}

\item{outfile}{Allows users to specify a prefix for an output folder. Name must a character string enclosed in either ``'' or `'.}

\item{pairwise}{Specified whether pairwise I\_{n} should be calculated.}

\item{xlsx}{A logical argument indicating whether results should be written to an \emph{xlsx} file. If \code{xlsx = FALSE} (default), results will be written to text files.}

\item{boots}{Determines the number of bootstrap iterations to be carried out. The default value is \code{boots = 0}, this is only valid when 
all bootstrap options are false.}

\item{para}{Allows for parallel computation of pairwise locus \emph{In}. The number of available core is automatically detected if \code{para = TRUE}.}
}

\value{

\code{inCalc} return a list object to the R workspace, with elements described below. In addition to this results can be optionally written to file using the \code{outfile} argument. If \code{xlsx = TRUE} results will be written to a multi-sheet xlsx file. If \code{xlsx = FALSE} results are written to multiple text file, the number of which depends on the function arguments used.

\item{global}{A \code{data.frame} containing the \emph{In} values for each locus, calculated across all samples in \code{infile}. If \code{boots} is an integer greater than 0, this data.frame will also contain lower and upper 95\% confidence limits for each locus.}

\item{pairwise}{A \code{data.frame} containing the pairwise locus \emph{In} values for all possible pairwise population comparisons. This object is returned when \code{boots} is an integer greater than 0.}

\item{lower_CI}{If \code{pairwise = TRUE} and \code{boots} is an integer greater than 0, \code{lower_CI} is returned. It is a \code{data.frame} containing the lower 95\% confidence limit for the corresponding pairwise estimate in the \code{pairwise} data.}

\item{upper_CI}{If \code{pairwise = TRUE} and \code{boots} is an integer greater than 0, \code{upper_CI} is returned. It is a \code{data.frame} containing the upper 95\% confidence limit for the corresponding pairwise estimate in the \code{pairwise} data.}

}

\note{
Since version 1.9.0, the speed of this function has been greatly improved. Users can expect up to x10 speed up on previous versions. The output data structure is also slightly different from v1.9.0 onwards.
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
\dontrun{
# To run an example use the following format
library(diveRsity)
data(Test_data)
Test_data[is.na(Test_data)] <- ""

test_results<-inCalc(infile = Test_data, outfile = 'out', pairwise = TRUE,
                     xlsx = FALSE, boots = 1000, para = TRUE)
}
}
