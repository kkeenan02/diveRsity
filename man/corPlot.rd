\name{corPlot}
\alias{corPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot the relationship between Gst, G'st, Theta, D (Jost) and the number of alleles at a locus.
}

\description{
\code{corPlot} uses the information calculated by \code{divPart} and \code{readGenepop} to plot and calculate the relationship between Gst, G'st, Theta, D (Jost) and the number of alleles per locus. This information can then be used to assess the likelyhood that data derived from the loci is suitable for the calculation of population demography. This is based on the assumption that where \emph{u} > \emph{m}, demographic process are obscured by mutation.
}

\details{
This function returns an R plot for Fst (theta) vs number of alleles, Gst vs number of alleles, G'st vs number of alleles and D (Jost) vs number of alleles. 
}

\usage{
corPlot(x, y)
}

\arguments{

\item{x}{
Results object returned from the function \code{readGenepop.user}
}

\item{y}{
Results object returned from the function \code{divPart}
}
}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}
