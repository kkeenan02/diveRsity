\name{diffCalc}
\alias{diffCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A faster function for calculating genetic differentiation statistics
}

\description{
This function allows the calculation of pairwise differentiation using a range of population statistics, such as Gst (Nei & Chesser, 1983), G'st (Hedrick, 2005), theta (Weir & Cockerham, 1984) and D (Jost, 2008). These parameters can also be calculated at the global and locus levels. Significance of differentiation can be assessed through the calculation of 95\% confidence limits using a bias corrected bootstrapping method. The functionality \code{diffCalc} is similar to the \code{fastDivPart} function. However \code{diffCalc} is much faster and more memory efficient than \code{fastDivPart}. This function also only allows results to be written to text files rather than xlsx file (as in \code{fastDivPart}. No plotting options are provide in \code{diffCalc}.)
}


\usage{
diffCalc(infile = NULL, outfile = NULL, fst = FALSE, pairwise = FALSE, 
         bs_locus = FALSE, bs_pairwise = FALSE, boots = NULL, 
         ci_type = "individuals", alpha = 0.05, para = FALSE)
}
    	   
\arguments{

\item{infile}{Specifying the name of the \emph{`genepop'} (Rousset, 2008) file from which the statistics are to be calculated This file can be in either the 3 digit of 2 digit format. See \url{http://genepop.curtin.edu.au/help_input.html} for detail on the genepop file format.}

\item{outfile}{A character string specifying the name of the folder to which results should be written.}

\item{fst}{A Logical argument indicating whether Weir & Cockerham's 1984 F-statistics should be calculated. NOTE - Calculating these statistics adds significant time to analysis when carrying out pairwise comparisons.}

\item{pairwise}{A logical argument indicating whether standard pairwise diversity statistics should be calculated and returned as a diagonal matrix.}

\item{bs_locus}{Gives users the option to calculate bias corrected 95\% confidence intervals for locus statistic species.}

\item{bs_pairwise}{Gives users the option to calculate bias corrected 95\% confidence intervals for pairwise statistics.}

\item{boots}{Specified the number of bootstraps for the calculation of 95\% confidence intervals.}

\item{ci_type}{A character string indicating whether bootstrapping should be carried out over individuals within samples (\code{``individuals''}.), or across loci (\code{``loci''}).}

\item{alpha}{A numeric argument, specifying the alpha value used to estimate confidence limits for relevant parameters. Both the alpha/2 and the 1-(alpha/2) quantiles will be returned. Default value results in 95\% CI.}

\item{para}{A logical argument indicating whether computations should be carried out over multiple CPUs, if available.}
}

\value{

If \code{outfile} is given as a character string, all results will be written to text files. The files will be written to a directory under the current working directory. The number of files written depends on the options choose. As well as this a list object is returned to the R workspace, containing the following results:

\item{std_stats}{A \code{data.frame}, containing locus estimates for Gst, G'st, G''st, D (and Weir and Cockerham's F-statistics, if \code{fst = TRUE}). The last row of this dataframe contains the global estimate for each statistic across all samples and loci.}

\item{global_bs}{If \code{bs_locus = TRUE}, this object is returned. It is a \code{dataframe} containing global estimates and lower and upper 95\% CIs for all relevant statistics.}

\item{bs_locus}{If \code{bs_locus = TRUE}, this object is returned. It is a list of either 4 or 7 dataframes the number of which depend on \code{fst}. Each dataframe contains the locus statistics estimate across all samples along with lower and upper 95\% confidence limits.}

\item{pw_locus}{If \code{pairwise = TRUE}, this list of dataframes is returned. The list contains a dataframe for each relevant statistics, where rows correspond to loci and columns correspond to all possible pairwise combinations of samples. If \code{outfile} is provided, these results are also written to file.}

\item{pairwise}{If \code{pairwise = TRUE}, this object is returned. It is a list of either 3 or 4 matrices (depending on the value of \code{fst}) containing pairwise estimate of relevant statistics.}

\item{bs_pairwise}{A list of either 4 or 5 (depending on the value of \code{fst}) containing pairwise estimates and lower and upper 95\% confidence intervals for all relevant statistics.}
}

\references{

Eddelbuettel, D., and Francois, R., (2011). Rcpp: Seamless R and C++ Integration. Journal of Statistical Software, 40(8), 1-18. URL http://www.jstatsoft.org/v40/i08/.

Hedrick, P., ``A standardized genetic differentiation measure,'' Evolution,
vol. 59, no. 8, pp. 1633-1638, (2005).

Jost, L., ``G ST and its relatives do not measure differentiation,'' Molec-
ular Ecology, vol. 17, no. 18, pp. 4015-4026, (2008).

Manly, F.J., ``Randomization, bootstrap and Monte Carlo methods in biology'', Chapman and Hall, London, 1997.

Meirmans, P.G., and Hedrick, P.W., (2011), Assessing population structure: Fst and related measures., Molecular Ecology, Vol. 11, pp5-18. doi: 10.1111/j.755-0998.2010.02927.x

Nei, M. and Chesser, R., ``Estimation of fixation indices and gene diver-
sities,'' Ann. Hum. Genet, vol. 47, no. Pt 3, pp. 253-259, (1983).

R Development Core Team (2011). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria. ISBN
3-900051-07-0, URL http://www.R-project.org/.

Rousset, F., ``genepop'007: a complete re-implementation of the genepop
software for Windows and Linux.,'' Molecular ecology resources, vol. 8,
no. 1, pp. 103-6, (2008).

Weir, B.S. & Cockerham, C.C., Estimating F-Statistics, for the Analysis of Population Structure, Evolution, vol. 38, No. 6, pp. 1358-1370 (1984). 
}
\author{Kevin Keenan <kkeenan02@qub.ac.uk>}

\examples{
\dontrun{
# simply use the following format to run the function
library(diveRsity)
data(Test_data)
Test_data[is.na(Test_data)] <- ""

test_result <- diffCalc(infile = Test_data, outfile = "myresults",
                        fst = TRUE, pairwise = TRUE, bs_locus = TRUE,
                        bs_pairwise = TRUE, boots = 1000, para = TRUE)
}
}
