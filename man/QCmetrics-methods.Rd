\name{QCmetrics-methods}
\docType{methods}
\alias{QCmetrics}
\alias{QCmetrics-methods}
\alias{QCmetrics,ChIPQCexperiment-method}
\alias{QCmetrics,ChIPQCsample-method}
\alias{QCmetrics,list-method}
\title{Retrieve consolidated set of QC metrics}
\description{
Retrieves a consolidated set of QC metrics.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieves a \code{matrix} of QC metrics for all the samples in an experiment, with a column of values for each sample.
}

\item{\code{signature(object = "list")}}{
Retrieves a \code{matrix} of QC metrics for all the samples in a list of ChIPQCsample objects, with a column of values for each sample.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieves a \code{vector} of QC metrics for a sample.
}
}}
\keyword{methods}

\author{Rory Stark and Thomas Carroll}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCsample}, \link{ChIPQCexperiment}
}
\examples{
data(example_QCexperiment)
QCmetrics(exampleExp)

data(tamoxifen_QC)
QCmetrics(tamoxifen)

QCmetrics(QCsample(tamoxifen,1))
}
