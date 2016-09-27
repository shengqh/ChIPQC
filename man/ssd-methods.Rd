\name{ssd-methods}
\docType{methods}
\alias{ssd-methods}
\alias{ssd}
\alias{ssd,ChIPQCexperiment-method}
\alias{ssd,ChIPQCsample-method}
\alias{ssd,list-method}
\title{Retrieve SSD (squared sum of deviations) values of peak coverage}
\description{
Retrieve SSD (squared sum of deviations) values of peak coverage density.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of SSD values, one for each sample in an experiment.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of SSD values, one for each sample in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve the SSD (squared sum of deviations) for a sample, computed from the standard deviation of the \code{\link{coveragehistogram}}.
}
}}

\note{
uses the SSD calculation from the \code{chipseq} package.
}
\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
ssd(exampleExp)
ssd(QCsample(exampleExp,1))
}