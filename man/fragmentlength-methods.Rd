\name{fragmentlength-methods}
\docType{methods}
\alias{fragmentlength-methods}
\alias{fragmentlength}
\alias{fragmentlength,ChIPQCexperiment-method}
\alias{fragmentlength,ChIPQCsample-method}
\alias{fragmentlength,list-method}
\title{Retrieve the estimated fragment length}
\description{
Retrieve the estimated fragment length.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of estimated fragments sizes, one for each sample in the experiment.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of estimated fragments sizes, one for each sample in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample", width)}}{
Retrieve the estimated fragment length for a sample. If \code{width} is missing, the readlength derived from the bam file is used as the read length.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
fragmentlength(exampleExp)
fragmentlength(QCsample(exampleExp,1))
}