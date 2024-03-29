% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General_methods.R
\name{PL}
\alias{PL}
\alias{PL,RHermesExp,numeric-method}
\title{PL}
\usage{
PL(struct, id)

\S4method{PL}{RHermesExp,numeric}(struct, id)
}
\arguments{
\item{struct}{An RHermesExp object}

\item{id}{The number corresponding to the PL you want to access}
}
\value{
The RHermes_PL object at position id.
}
\description{
PL
}
\examples{
\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
                             package = "RHermes"))}
PL(struct, 1)
}
\seealso{
Other Getters: 
\code{\link{DB}()},
\code{\link{IL}()},
\code{\link{Ident}()},
\code{\link{MS2Data}()},
\code{\link{SOI}()},
\code{\link{adlist}()}
}
\author{
Roger Gine
}
\concept{Getters}
