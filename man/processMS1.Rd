% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MS1_functions.R
\name{processMS1}
\alias{processMS1}
\alias{processMS1,RHermesExp,character-method}
\title{processMS1}
\usage{
processMS1(struct, files, labelled = FALSE)

\S4method{processMS1}{RHermesExp,character}(struct, files, labelled = FALSE)
}
\arguments{
\item{struct}{RHermesExp S4 object to update with the processed data.
Important: The objects needs to have the experimental parameters already
set before calling this function.}

\item{files}{Character vector of all the paths to the files to be processed.}

\item{labelled}{Logical, whether to check for all 13C isotopic signals.
Defaults to FALSE}
}
\value{
An RHermesExp object with the processed PLs.
}
\description{
Main function to process the input MS1 mzML files. It
accesses the files with mzR and processes them with RHermes internal
functions. For each file a PL object is generated inside the
object_data_PL slot.
}
\examples{
\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
                             package = "RHermes"))
         MS1file <- system.file("extdata", "MS1TestData.mzML",
                             package = "RHermes")
         }
processMS1(struct, MS1file, FALSE) #For regular use
processMS1(struct, MS1file, TRUE) #For labelled data

}
