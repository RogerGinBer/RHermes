% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot_functions.R
\name{plotIL}
\alias{plotIL}
\alias{plotIL,RHermesExp,numeric-method}
\title{plotIL}
\usage{
plotIL(struct, id)

\S4method{plotIL}{RHermesExp,numeric}(struct, id)
}
\arguments{
\item{struct}{An RHermesExp object}

\item{id}{Number of the inclusion list to plot}
}
\value{
An interactive plot_ly object
}
\description{
Draws a mz-rt representation of the inclusion list entries,
coloured by their intensity value.
}
\examples{
\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
                             package = "RHermes"))}
plotIL(struct, 1)
}
\seealso{
Other Plots: 
\code{\link{plotCoverage}()},
\code{\link{plotISF}()},
\code{\link{plotPL}()}
}
\author{
Roger Gine
}
\concept{Plots}
