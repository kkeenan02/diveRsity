\name{microPlexer}
\alias{microPlexer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Launches a \code{shiny} app for the arrangement of microsatellite loci into size and fluorophore based multiplex groups.
}

\description{
This function will launch a browser based app designed using the package, \code{shiny}, which allows users to design microsatellite multiplex groups. 
}

\details{
The application provides flexibility in marker organisation through the use of two distinct algorithms. The \emph{high-throughput} algorithm will attempt to group as many loci into as few multiplex groups as allowable based on locus size and the available fluorophore labels, while the \emph{balanced} algorithm attempts to organise loci into multiplex group of roughly equal density, so as to offset possible primer interactions etc. As input, the application accepts a \code{.csv} file with three named columns. The structure of this file is as follows:

\code{nms} - contains the names of loci

\code{lower} - contains the lower size range of loci

\code{upper} - contains the upper size range of loci



}

\usage{
microPlexer()
}


\author{
Kevin Keenan <kkeenan02@qub.ac.uk>
}

