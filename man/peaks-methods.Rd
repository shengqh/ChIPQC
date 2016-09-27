\name{peaks-methods}
\docType{methods}
\alias{peaks}
\alias{peaks-methods}
\alias{peaks,ChIPQCexperiment-method}
\alias{peaks,ChIPQCsample-method}
\alias{peaks,list-method}
\title{Retrieve peaks}
\description{
Retrieve peaks.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{\link{GRangesList}} of the peaks associated with all the samples in an experiment.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{\link{GRangesList}} of the peaks associated with all the samples in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve a \code{\link{GRanges}} object containing the peaks associated with a sample.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
expPeaks = peaks(exampleExp)
length(expPeaks)
peaks(QCsample(exampleExp,1))
}