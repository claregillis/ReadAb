% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readAntibody.R
\name{.IsValidChain}
\alias{.IsValidChain}
\title{\strong{All following functions are helpers}
Determine if a chain ID is valid}
\usage{
.IsValidChain(chainIdentifier, realChains)
}
\arguments{
\item{chainIdentifier}{a chain identifier object passed by the user}

\item{realChains}{the list of chains present in the PDB file passed by the user}
}
\value{
True iff the chainIdentifier is a single character present in
realChains
}
\description{
Determine whether a chain identifier is valid (a character identifier that
is present in the list of real chain identifiers)
}
\details{
Chain identifiers in PDB files must be a single character
}
\examples{
.isValidChain('B', c('A', 'B', 'C'))

.isValidChain('D', c('A', 'B', 'C'))

}
