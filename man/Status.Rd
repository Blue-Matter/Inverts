% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Status.R
\name{Status}
\alias{Status}
\title{Summarize current stock status}
\usage{
Status(Hist, plot = T, col = "cadetblue", digits = 4)
}
\arguments{
\item{Hist}{An object of class Hist created by runMSE(OM, Hist=T)}

\item{plot}{Logical, should the results be plotted?}

\item{col}{Color of bars}
}
\description{
Histograms of spawning biomass and fishing mortatlity rate relative to reference levels
}
\examples{
Status(Hist)
}
\author{
T. Carruthers
}
