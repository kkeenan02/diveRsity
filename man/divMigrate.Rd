\name{divMigrate}
\alias{divMigrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An experimental function for the detection of directional differentiation from microsatellite data.
}

\description{
\code{divMigrate} uses the method described in Sundqvist \emph{et al.,} 2013 to plot the relative migration levels between population samples, from microsatellite allele frequency data. The method is still in the experimental stages, and it is not clear how well it performs under most evolutionary scenarios. Caution should be exercised.
}

\details{
The function will except both Arlequin (.arp) genotype and genepop (.gen/.txt) files, containing co-dominant diploid data. Using the method outlined in Sundqvist \emph{et al.,} 2013, the relative migration levels between all pairs of populations is determined. A weighted network plot is returned, as well as four matrices containing the objects described below.
}

\usage{
divMigrate(infile = NULL, stat = "d_jost")
}
\arguments{

\item{infile}{Should either be a character string of the input file name or pathway to the input file. This file should be in either the Arlequin genotype file format (.arp), or in the genepop file format.}

\item{stat}{A character vector of the value \code{stat = 'gst'}, \code{stat = 'd_jost'} (default), or \code{stat = c('gst', 'd_jost')}. This argument specifies which statistic should be used to calculate relative migration between samples.}

}

\value{

\item{D_mig}{A pairwise matrix containing the normalized relative migration values for all population pairs, as determined by $D_{Jost}$ values}

\item{Gst_mig}{A pairwise matrix containing the normalized relative migration values for all population pairs, as determined by $G_{ST}$ values}

}

\references{
Lisa Sundqvist, M.Z. & Kleinhans, D., 2013. Directional genetic differentiation and asymmetric migration. arXiv pre-print: arXiv:1304.0118v2

Rousset, F., ``genepop'007: a complete re-implementation of the genepop
software for Windows and Linux.,'' Molecular ecology resources, vol. 8,
no. 1, pp. 103-6, (2008).
}

\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}