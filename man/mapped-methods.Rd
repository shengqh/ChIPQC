\name{mapped-methods}
\docType{methods}
\alias{mapped-methods}
\alias{mapped}
\alias{mapped,ChIPQCexperiment-method}
\alias{mapped,ChIPQCsample-method}
\alias{mapped,list-method}
\title{Retrieve numbers of mapped reads}
\description{
Retrieve the numbers of mapped reads.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of the numbers of mapped (aligned) reads for each sample in an experiment.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of the numbers of mapped (aligned) reads for each sample in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve the number of mapped (aligned) reads in a sample.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
mapped(exampleExp)
mapped(QCsample(exampleExp,1))
}

