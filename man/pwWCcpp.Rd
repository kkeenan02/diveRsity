\name{pwWCcpp}
\alias{pwWCcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{C++ level function for pairwise Weir & Cockerham F-stats}
\description{C++ level function for pairwise Weir & Cockerham F-stats}

\usage{
pwWCcpp(hsum1, af1, indtyp1, pw)
}

\arguments{
  \item{hsum1}{hsum returned by internal function \code{statCalc}}
  \item{af1}{af returned by internal function \code{statCalc}}
  \item{indtyp1}{indtyp returned by internal function \code{statCalc}}
  \item{pw}{Integer matrix indicating pairwise indexes}
}
\keyword{internal}