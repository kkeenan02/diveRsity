\name{divMigrateOnline_test}
\alias{divMigrateOnline_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An experimental function for the detection of directional differentiation from microsatellite data.
}

\description{
\code{divMigrateOnline_test} uses the method described in Sundqvist \emph{et al.,} 2013 to plot the relative migration levels between population samples, from microsatellite allele frequency data. The method is still in the experimental stages, and it is not clear how well it performs under most evolutionary scenarios. Caution should be exercised.
}

\details{
The function will except both Arlequin (.arp) genotype and genepop (.gen/.txt) files, containing co-dominant diploid data. Using the method outlined in Sundqvist \emph{et al.,} 2013, the relative migration levels between all pairs of populations is determined. A weighted network plot is returned, as well as four matrices containing the objects described below.
}

\usage{
divMigrateOnline_test(infile = NULL, nbs = 0, stat = "all", para = FALSE, 
                 true.nm = FALSE, alpha = 0.05)
}
\arguments{

\item{infile}{Input genepop file name or path.}

\item{nbs}{The number of bootstrap iteration to carry out for calculating mean relative migration and confidence intervals.}

\item{stat}{A character string indicating which statistic should be used to estimate relative migration between populations. The argument accepts one of the following: "d" (Jost's D), "gst" (Nei's Gst) or "Nm" (Alcala et al, 2014).}

\item{para}{A logical argument, specifying if multiple CPUs should be used when available.}

\item{true.nm}{A logical indicating whether Nm (effective number of migrants) as estimated by Alcala et al. (2014) should be returned.}

\item{alpha}{The alpha level to use when constructing confidence intervals. The alpha/2 and 1-(alpha/2) quantiles are used when testing the significance of directional gene flow}

}

\value{

\item{Relative migration}{Relative migration matrices}

\item{Significant directional migration}{If \code{nbs > 0} the significance of pairwise relative migration is tested (i.e. non-overlapping 95\% CIs). All non-significant values in the standard relative migration matrices are replaced with 0.}

}

\references{
Lisa Sundqvist, M.Z. & Kleinhans, D., 2013. Directional genetic differentiation and asymmetric migration. arXiv pre-print: arXiv:1304.0118v2

Rousset, F., ``genepop'007: a complete re-implementation of the genepop
software for Windows and Linux.,'' Molecular ecology resources, vol. 8,
no. 1, pp. 103-6, (2008).

Alcala, N., Goudet, J., Vuilleumier, S., 2014, On the transition of genetic differentiation from isolation to panmixia: What we can learn from Gst and D, Theoretical Population Biology, Vol 93, pp75-84.

}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}

\keyword{internal}