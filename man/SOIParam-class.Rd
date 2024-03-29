% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classes.R, R/SOI_methods.R
\docType{class}
\name{SOIParam-class}
\alias{SOIParam-class}
\alias{SOIParam}
\alias{show,SOIParam-method}
\title{SOIParam}
\usage{
\S4method{show}{SOIParam}(object)
}
\arguments{
\item{object}{A SOIParam object}
}
\value{
A SOIParam object
}
\description{
The SOIParam class contains all details regarding SOI detection,
 such as the number of RT bins and their width, minimum data point intensity,
 usage of blank substraction, etc.
}
\details{
It's a part of the RHermesSOI object, and is an input for findSOI().
 This object is usually generated by getSOIParam(), which features
 already-made templates.
}
\section{Slots}{

\describe{
\item{\code{specs}}{A data frame containing the RT bin information. Each row
represents a binning step that will be applied to the annotated data points
to detect the SOIs.}

\item{\code{maxlen}}{Maximum SOI length in seconds. If the SOI is longer tha maxlen, a
broad peak picking will divide the long SOI into smaller sized SOIs.}

\item{\code{minint}}{Minimum data point intensity. If a point is smaller than minint,
it will not be considered in the SOI detection.}

\item{\code{blanksub}}{Logical. Whether to perform blanksub. It should NOT be entered
by the user in this point, since findSOI will fill in this information.}

\item{\code{blankname}}{Character. The corresponding blank file name. As before, it
should not be entered by the user, findSOI will fill it automatically.}
}}

\examples{
SOIParam()
}
\seealso{
RHermesSOI findSOI getSOIParam

Other Params: 
\code{\link{ExpParam-class}},
\code{\link{ILParam-class}}
}
\author{
Roger Gine
}
\concept{Params}
