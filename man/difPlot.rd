\name{difPlot}
\alias{difPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot pairwise statistics calculated by divPart.
}

\description{
This function uses results from \code{divPart} to plot pairwise statistic estimators. 
}

\details{
This function returns four R plots for Fst (theta) (if \code{WC_Fst = TRUE} in \code{divPart}, Gst, G'st and D (Jost) pairwise values. This plot printed to the R graphical device if the argument \code{interactive = FALSE}. N.B. the argument \code{interactive = TRUE} is only valid if the package \code{sendplot} is available.
}

\usage{
difPlot(x, outfile = NULL, interactive = FALSE)
}

\arguments{

\item{x}{Results object returned from the function 'divPart'}
\item{outfile}{A character string indication the folder location to which plot files should be writen.}
\item{interactive}{A logical argument indication whether the package 'sendplot' should be used to plot 'divPart' pairwise results.}
}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}
