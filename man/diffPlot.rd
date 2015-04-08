\name{diffPlot}
\alias{diffPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot pairwise statistics calculated by divPart.
}

\description{
This function uses results from \code{divPart} to plot pairwise statistic estimators. 
}

\details{
A number of heatmap style plots of pairwise differentiation are generated. The function takes the output from either \code{fastDivPart} or \code{diffPlot} and writes interactive \code{HTML} plots to file. The number of plots depends on the input structure. For instance, if \code{fst = FALSE} in either \code{fastDivPart} or \code{diffCalc}, then a plot containing pairwise Fst will be produced.
}

\usage{
diffPlot(x, outfile = NULL, interactive = FALSE)
}

\arguments{

\item{x}{Results object returned from the function 'divPart'}
\item{outfile}{A character string indication the folder location to which plot files should be written.}
\item{interactive}{A logical argument indication whether the package 'sendplot' should be used to plot 'divPart' pairwise results.}
}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}
