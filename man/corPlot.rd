\name{corPlot}
\alias{corPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot the relationship between Gst, G'st, Theta, D (Jost) and the mean number of alleles at a locus. A reimplementation of \code{corPlot}. This function will replace the old function in later version of the package.
}

\description{
\code{corPlot} uses the information calculated by \code{diffCalc}  to plot and calculate the relationship between Gst, G'st, Theta, D (Jost) and the mean number of alleles per locus. This information can then be used to assess the likelihood that data derived from the loci are suitable for the calculation of population demography using Fst or its analogues. This is based on the assumption that where \emph{u} > \emph{m}, demographic process are obscured by mutation.
}

\details{
This function returns four scatter plots showing the relationship between Fst (Weir and Cockerham 1984), Gst (Nei and Chesser 1983), G'st (Hedrick 2005) and D (Jost 2008) and the mean number of alleles per locus. The function allows user to write plots to file, and return a four panelled figure. Plots are generated using the \code{ggplot2} package and arranged in a grid using the \code{multiplot} function by Winston Chang (Chang 2012).
}

\usage{
corPlot(infile = NULL, write = FALSE, plot.format = NULL)
}

\arguments{

\item{infile}{A character string indicating the location and name of a genepop format file to be read. If the file is in the current working directory, only the name must be provided. If the file is in a directory other than the current working directory, either a relative or absolute path to the file must be provided. The genepop file can be in the 2-digit or 3-digit allele format.}

\item{write}{A logical argument indicating whether results should be written to file or not}

\item{plot.format}{A string indicating the format to which plots should be written. Either 'png' or 'eps' are accepted.}

}

\references{

Cheng, W., R Graphics Cookbook, O'Reilly Media inc. 2012.

Hedrick, P., ``A standardized genetic differentiation measure,'' Evolution,
vol. 59, no. 8, pp. 1633-1638, (2005).

Jost, L., ``G ST and its relatives do not measure differentiation,'' Molec-
ular Ecology, vol. 17, no. 18, pp. 4015-4026, (2008).

Nei, M. and Chesser, R., ``Estimation of fixation indices and gene diver-
sities,'' Ann. Hum. Genet, vol. 47, no. Pt 3, pp. 253-259, (1983).

Weir, B.S. & Cockerham, C.C., Estimating F-Statistics, for the Analysis of Population Structure, Evolution, vol. 38, No. 6, pp. 1358-1370 (1984). 

}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}