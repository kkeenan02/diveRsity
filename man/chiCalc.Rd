\name{chiCalc}
\alias{chiCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing sample independence from genotype counts
}

\description{
\code{chiCalc} carries out Fisher's exact tests of sample independence using genotype data.
}

\details{
All results will be written to a user defined folder ("working_directory/outfile"), providing an argument is passed for 'outfile'. Otherwise, results will only be returned to the workspace.

Fisher's exact tests are carried out using the function \code{fisher.test}. Multilocus p values are calculated using Fisher's method for combining p value.
}

\usage{
chiCalc(infile = NULL, outfile = NULL, pairwise = FALSE, mcRep = 2000)
}

\arguments{

\item{infile}{A character string indicating the location and name of a genepop format file to be read. If the file is in the current working directory, only the name must be provided. If the file is in a directory other than the current working directory, either a relative or absolute path to the file must be provided. The genepop file can be in the 2-digit or 3-digit allele format.}

\item{outfile}{A character string indicating the prefix to be added to the results directory created. All results files will be written to this directory.}

\item{pairwise}{A logical argument indicating whether sample independence should be tested between all population pairs.}

\item{mcRep}{An integer specifying the number Monte Carlo test replicates. See \code{?fisher.test} for more information.}

}

\value{

\item{overall}{A data frame containing p values calculated across all population samples, per locus and across all loci.}

\item{multilocus_pw}{Generated if \code{pairwise = TRUE}. The object is a data frame containing multilocus p value calculated for all population pairs.}

\item{locus_pw}{A dataframe containing locus p values calculated for all pairs of populations. Rows represent loci, while columns represent pairwise comparisons.}

}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}

\examples{
\dontrun{
# To run an example use the following format
library(diveRsity)
data(Test_data)
test_results <- chiCalc(infile = Test_data, outfile = NULL,
                        pairwise = TRUE, mcRep = 5000)
}
}
