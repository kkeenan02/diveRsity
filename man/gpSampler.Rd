\name{gpSampler}
\alias{gpSampler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Randomly sample a genepop file
}
\description{
Randomly re-samples population samples from a genepop file for user defined sizes and returns a new genepop file containing the sub-sample data. The function can easily be integrated into simulation pipelines when exploring sample size effects.
}
\usage{
gpSampler(infile = NULL, samp_size = 10, outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{A character string pointing to a genepop input file. If the file exists in the current working directory, only the file name is required. If the file resides elsewhere, the entire file path is required.}
  \item{samp_size}{This argument specifies the number of individuals to be randomly sampled from each population sample in the input file. Either a single numeric argument can be passed or a numeric vector of length equal to the number of population samples in the genepop file.}
  \item{outfile}{A character string specifying the prefix name for the resulting file. This string will be suffixed with ".gen". The file will be written to the working directory, unless a properly formatted path is passed to \code{outfile}}
}

\author{
Kevin Keenan
}