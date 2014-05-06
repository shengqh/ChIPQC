\name{QCsample-methods}
\docType{methods}
\alias{QCsample}
\alias{QCsample-methods}
\alias{QCsample,ChIPQCexperiment-method}
\title{Retrieve sample objects associated with an experiment}
\description{
Retrieves a list of \code{\link{ChIPQCsample}} objects, or one specific \code{\link{ChIPQCsample}} object, from a \code{\link{ChIPQCexperiment}} object. 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", sampleID)}}{

If \code{sampleID} is missing or equal to \code{0}, the full \code{list} of \code{\link{ChIPQCsample}} objects is returned. If \code{sampleID} is an integer \code{n}, the \code{\link{ChIPQCsample}} object corresponding to the \code{n}th sample is returned. If \code{sampleID} is a character string, the \code{\link{ChIPQCsample}} object corresponding to the sample with that ID is returned.
}
}}

\keyword{methods}

\author{Rory Stark and Thomas Carroll}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCsample}, \link{ChIPQCexperiment}
}
\examples{
data(example_QCexperiment)
samples = QCsample(exampleExp)
length(samples)
names(samples)
samples$CTCF_1
}