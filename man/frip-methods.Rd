\name{frip-methods}
\docType{methods}
\alias{frip}
\alias{frip-methods}
\alias{frip,ChIPQCexperiment-method}
\alias{frip,ChIPQCsample-method}
\alias{frip,list-method}
\title{Retrieve fraction of reads in peaks}
\description{
Retrieve the fraction of reads in peaks
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of values representing the proportion of reads that overlap peaks for each sample in an experiment.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of values representing the proportion of reads that overlap peaks for each sample in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve a value representing the proportion of reads that overlap the peaks for a sample.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
frip(exampleExp)
frip(QCsample(exampleExp,1))
}