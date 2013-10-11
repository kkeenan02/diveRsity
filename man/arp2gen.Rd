\name{arp2gen}
\alias{arp2gen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast and simple conversion of arlequin (.arp) genotype files to genepop files, form within the \code{R} environment.}

\description{
\code{arp2gen} allows simple file conversion from the arlequin genotype format to the genepop genotype format. Arlequin files can gave 2-digit or 3-digit allele records. No file size limit is imposed, however, system RAM is a limiting factor.
}

\details{
Following the input of a .arp file, \code{arp2gen} will write a .gen file to the same directory as the original .arp file. The output file will have the same name as the original \code{infile}, with the exception of the file extension, which will be .gen following file conversion. The original .arp file will remain unmodified.
}

\usage{
arp2gen(infile)
}
\arguments{

\item{infile}{Specifying the name of the \emph{`.arp'} arlequin genotype file. The argument must be a character string of the file name if it is located in the current working directory, or the file path (relative or absolute) if not.}

}


\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}