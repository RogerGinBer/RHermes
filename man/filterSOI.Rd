% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SOI_functions.R
\name{filterSOI}
\alias{filterSOI}
\alias{filterSOI,RHermesExp,numeric-method}
\title{filterSOI}
\usage{
filterSOI(struct, id, minint = 10000, isofidelity, minscore = 0.8)

\S4method{filterSOI}{RHermesExp,numeric}(struct, id, minint = 10000, isofidelity, minscore = 0.8)
}
\arguments{
\item{struct}{The RHermesExp object to which the SOI list will be saved.}

\item{id}{ID of the SOI list to be filtered/checked.}

\item{minint}{Minimun SOI intensity. All SOIs below this value will be removed
from the SOI list}

\item{isofidelity}{Boolean. Whether to perform an isotopic fidelity check.}

\item{minscore}{Numeric. Minimum value (between 0 and 1) of isofidelity to
retain an entry. Defaults to 0.8.}
}
\value{
A filtered SOI list.
}
\description{
Performs a series of filters and quality checks to a given SOI
 list, removing unwanted SOIs in the process.
}
\examples{
\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
                             package = "RHermes"))}
struct <- filterSOI(struct, id = 1, minint = 10000, isofidelity = TRUE)
}
\author{
Roger Gine
}
