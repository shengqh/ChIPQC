\name{ribl-methods}
\docType{methods}
\alias{ribl-methods}
\alias{ribl}
\alias{ribl,ChIPQCexperiment-method}
\alias{ribl,ChIPQCsample-method}
\alias{ribl,list-method}
\title{Retrieve numbers of reads overlapping blacklisted regions}
\description{
Retrieve the numbers of reads overlapping blacklisted regions.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of the numbers of reads overlapping blacklisted regions for each sample in an experiment.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of the numbers of reads overlapping blacklisted regions for each sample in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve the number of reads overlapping blacklisted regions in a sample.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
ribl(exampleExp)
ribl(QCsample(exampleExp,1))
}

