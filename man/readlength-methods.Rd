\name{readlength-methods}
\docType{methods}
\alias{readlength-methods}
\alias{readlength}
\alias{readlength,ChIPQCexperiment-method}
\alias{readlength,ChIPQCsample-method}
\alias{readlength,list-method}
\title{Retrieve read length values}
\description{
Retrieve read length values.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of read length values, one for each sample in an experiment.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of read length values, one for each sample in a list of ChIPQCsample objects.
}


\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve the read length value for a sample.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
readlength(exampleExp)
readlength(QCsample(exampleExp,1))
}
