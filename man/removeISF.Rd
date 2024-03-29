% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SOI_functions.R
\name{removeISF}
\alias{removeISF}
\title{removeISF}
\usage{
removeISF(
  struct,
  id,
  DBpath = "D:/MS2ID_B2R_20201113_083214.rds",
  justAnnotate = FALSE
)
}
\arguments{
\item{struct}{The RHermesExp object to which the SOI list will be saved.}

\item{id}{ID of the SOI list to be filtered/checked.}

\item{DBpath}{Path to a MS2 spectral database with a specific format}

\item{justAnnotate}{Whether to actually remove the ISF or just keep them
annotated. Defaults to FALSE (removes the ISF). If set to TRUE, you can then
plot the ISF network with plotISF()}
}
\value{
An updated RHermesExp object where the selected SOI has had its ISF
removed. If using justAnnotate = TRUE, the ISF aren't removed.
}
\description{
Detect and remove ISF signals from a SOI list using low
collision energy MS2 data
}
\examples{
if(FALSE){
    #You need a MS2 Database with the proper format
    #(we are working on releasing one to the public soon enough)
    struct <- removeISF(struct, id = 1, DBpath = "MS2Database.rds")
}
}
\seealso{
plotISF
}
