\name{ReadLengthReadLengthCrossCoverage-methods}
\docType{methods}
\alias{ReadLengthCrossCoverage}
\alias{ReadLengthCrossCoverage-methods}
\alias{ReadLengthCrossCoverage,ChIPQCexperiment-method}
\alias{ReadLengthCrossCoverage,ChIPQCsample-method}
\alias{ReadLengthCrossCoverage,list-method}

\title{Retrieve the cross coverage values without extending reads}
\description{
Retrieve the cross coverage values without extending reads.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of cross-coverage values for all samples in a ChIP-seq experiment, with no shift.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of cross-coverage values for all samples  in a list of ChIPQCsample objects, with no shift.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve the cross-coverage value for a sample, with no shift.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
ReadLengthCrossCoverage(exampleExp)
ReadLengthCrossCoverage(QCsample(exampleExp,1))
}