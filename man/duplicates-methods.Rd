\name{duplicates-methods}
\docType{methods}
\alias{duplicates-methods}
\alias{duplicates}
\alias{duplicates,ChIPQCexperiment-method}
\alias{duplicates,ChIPQCsample-method}
\alias{duplicates,list-method}
\title{Retrieve numbers of duplicate reads.}
\description{
Retrieve the numbers of duplicate reads.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", bFiltered)}}{
Retrieve a \code{vector} of the numbers of duplicate reads for each sample in an experiment. A read is considered duplicated if another read maps to the same genomic location. If \code{bFiltered=TRUE} (or is missing), this will be the number of duplicates that pass the mapping quality filter for each sample. if \code{bFiltered=FALSE}, it will be the total number of duplicates for each sample.
}

\item{\code{signature(object = "list", bFiltered)}}{
Retrieve a \code{vector} of the numbers of duplicate reads for each sample in a list of ChIPQCsample objects. A read is considered duplicated if another read maps to the same genomic location. If \code{bFiltered=TRUE} (or is missing), this will be the number of duplicates that pass the mapping quality filter for each sample. if \code{bFiltered=FALSE}, it will be the total number of duplicates for each sample.
}

\item{\code{signature(object = "ChIPQCsample", bFiltered)}}{
Retrieve the number of duplicates for a sample. A read is considered duplicated if another read maps to the same genomic location. If \code{bFiltered=TRUE} (or is missing), this will be the number of duplicates that pass the mapping quality filter. if \code{bFiltered=FALSE}, it will be the total number of duplicates for the sample.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
duplicates(exampleExp)
duplicates(QCsample(exampleExp,1))
}

