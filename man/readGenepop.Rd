\name{readGenepop}
\alias{readGenepop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate allele frequencies from genepop files.
}

\description{
\code{readGenepop} allows the calculation of various parameters from 3 digit and 2 digit genepop files. The purpose of the function is mainly as a data manipulation process to allow for easy downstream analysis.
}

\details{
Results from this function allow for the calculation of various population genetics statistics, such as those calculated by \code{div.part} and \code{in.calc}. Users may find it useful for data exploration. For instance by employing the \code{plot {graphics}} function, an \emph{'ad hoc'} assessment of allele size distribution can be carried out using the code in the example section below.
From this example it is clear that the function will be particularly useful for those wishing to develop their own novel analysis methods. 
}

\usage{
readGenepop(infile = NULL, gp = 3, bootstrap = FALSE)
}

\arguments{

\item{infile}{Specifies the name of the \emph{`genepop'}(Rousset, 2008) file from which the statistics are to be calculated. This file can be in either the 3 digit of 2 digit format, and must contain only one \emph{whitespace} separator (e.g. ``space" or ``tab") between each column including the individual names column. The number of columns must be equal to the number of loci + 1 (the individual names column). If this file is not in the \code{working directory} the file path must be given. The name must be a character string (i.e. enclosed in ``" or `').}

\item{gp}{A numeric argument specifying the format of the \code{infile}. Either `3' or `2' are accepted as arguments. Default is \code{ gp = 3}.}

\item{bootstrap}{A logical argument specifying whether the user would like the \code{infile} data bootstrapped. If \code{bootstrap = TRUE} a genepop format object is returned. See \code{bootstrap_file} in the value section below.}

}

\value{

\item{npops}{The number of population samples in \code{infile}.}
\item{nloci}{The number of loci in \code{infile}.}
\item{pop_alleles}{A list of matrices (n = 2 x \code{npops}) containing haploid allele designations. Every two matrices contain the two alleles per individual per population. For example \code{pop_alleles[[1]][1,1]} and \code{pop_alleles[[2]][1,1]} are the two alleles observed in individual `1' in population `1' at locus `1', whereas \code{pop_alleles[[3]][1,1]} and \code{pop_alleles[[4]][1,1]} are the two alleles observed in individual `1' in population `2' at locus `1'.}
\item{pop_list}{A list of matrices (n = \code{npops}) containing the diploid genotypes of individuals per locus. }
\item{loci_lames}{A character vector containing the names of loci from \code{infile}.}
\item{pop_pos}{A numeric vector or the row index locations of the first individual per population in \code{infile}.}
\item{pop_sizes}{A numeric vector of length \code{npops} containing the number of individuals per population sample in \code{infile}.}
\item{allele_names}{A list of \code{npops} lists containing \code{nloci} character vectors of alleles names per locus. Useful for identifying unique alleles.}
\item{all_alleles}{A list of \code{nloci} character vectors of all alleles observed across all population samples in \code{infile}. }
\item{allele_freq}{A list containing \code{nloci} matrices containing allele frequencies per alleles per population sample.}
\item{raw_data}{An unaltered data frame of \code{infile}.}
\item{loci_harm_N}{A numeric vector of length \code{nloci}, containing the harmonic mean number of individuals genotyped per locus.}
\item{n_harmonic}{A numeric value representing the harmonic mean of \code{npops}.}
\item{pop_names}{A character vector containing a four letter population sample name for each population in \code{infile} (the first four letter of the first individual).}
\item{indtyp}{A list of length \code{nloci} containing character vectors of length \code{npops}, indicating the number of individuals per population sample typed at each locus.}
\item{nalleles}{A vector of the total number of alleles observed at each locus.}
\item{bs_file}{A genepop format data frame of bootstrapped \code{infile}. This value is only returned if \code{bootstrap = TRUE}.}
}

\references{
Rousset, F., ``genepop'007: a complete re-implementation of the genepop
software for Windows and Linux.,'' Molecular ecology resources, vol. 8,
no. 1, pp. 103-6, (2008).
}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}

\examples{
# Code to plot ordered allele fragment sizes to assess mutation model.
data(Test_data, package = "diveRsity") # define data
x <- readGenepop(infile = Test_data, gp = 3, bootstrap = FALSE)
locus10_pop1 <- c(x$pop_alleles[[1]][[2]][,10],
                  x$pop_alleles[[1]][[2]][,10])
sort_order <- order(locus10_pop1, decreasing = FALSE) #sort alleles
plot(locus10_pop1[sort_order], col="red", ylab = "Allele size")
}
