% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Status.R
\name{DFO_status}
\alias{DFO_status}
\title{Summarize current stock status according to standard DFO reference regions}
\usage{
DFO_status(Hist, legpos = "left", legpos2 = "bottomleft")
}
\arguments{
\item{Hist}{An object of class Hist created by runMSE(OM, Hist=T)}

\item{legpos}{Position of legend 1 (see ?legend)}

\item{legpos2}{Position of legend 2 (see ?legend)}
}
\description{
Distribution of estimated spawning biomass and fishing mortality rate relative to reference levels
}
\examples{
DFO_status(Hist)
}
\author{
T. Carruthers
}
