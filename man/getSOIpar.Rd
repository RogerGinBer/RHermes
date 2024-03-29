% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SOI_methods.R
\name{getSOIpar}
\alias{getSOIpar}
\alias{getSOIpar,ANY-method}
\title{getSOIpar}
\usage{
getSOIpar(tag = "double", mode = "regular", cwp = NA)

\S4method{getSOIpar}{ANY}(tag = "double", mode = "regular", cwp = NA)
}
\arguments{
\item{tag}{A character string that tells which premade SOI parameter
object to use. Currently the following tags are available: 'single',
'double', 'triple', and their extended counterparts for longer
chromatography experiments, 'single-x', 'double-x' and 'triple-x'.
These are all stored in /app/www/SOIFilterParams.csv, feel free to
locally change them or add new ones for your use (if you know what
you're doing).}

\item{mode}{Whether SOI detection should use the regular density-based
algorithm or xcms peak detection for defining the SOIs}

\item{cwp}{A CentWaveParam object used for either SOI detection (xcms mode)
or long SOI splitting (regular mode)}
}
\value{
A SoiParam object
}
\description{
Returns a SOIParam object for posterior SOI detection.
This function does NOT require the RHermesExp object as multple SOIParam
can be used at once. See \link[RHermes]{findSOI} for more info.
}
\examples{
if(FALSE){
par <- getSOIpar('double')
par2 <- getSOIpar('triple-x') #Etc. etc.
}
}
