\name{rip-methods}
\docType{methods}
\alias{rip-methods}
\alias{rip}
\alias{rip,ChIPQCexperiment-method}
\alias{rip,ChIPQCsample-method}
\alias{rip,list-method}
\title{Retrieve numbers of reads overlapping peaks}
\description{
Retrieve the numbers of reads overlapping peaks.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of the numbers of reads overlapping peaks for each sample in an experiment.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of the numbers of reads overlapping peaks for each sample in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve the number of reads overlapping peaks in a sample.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
rip(exampleExp)
rip(QCsample(exampleExp,1))
}

