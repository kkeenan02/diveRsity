\name{corPlotter}
\alias{corPlotter}
\title{
An internal function for plotting the relationship between differentiation and locus polymorphism
}
\description{
An internal function for plotting the relationship between differentiation and locus polymorphism
}
\usage{
corPlotter(x, y, dat, write = FALSE, plot.format = NULL, yname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A character string indicating the column name containing mean number of alleles per locus}
  \item{x}{A character string indicating the column name of the differentiation values per locus to be plotted}
  \item{dat}{A dataframe from which variable x and y are sourced}
  \item{write}{A logical argument, indicating whether results should be written to file}
  \item{plot.format}{A string indicating the format to be used when writing results to file}
  \item{yname}{A string or expression to be used as the yaxis label on the scatter plot}
}
\keyword{internal}