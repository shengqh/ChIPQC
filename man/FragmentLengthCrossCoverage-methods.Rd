\name{ReadLengthFragmentLengthCrossCoverage-methods}
\docType{methods}
\alias{FragmentLengthCrossCoverage}
\alias{FragmentLengthCrossCoverage-methods}
\alias{FragmentLengthCrossCoverage,ChIPQCexperiment-method}
\alias{FragmentLengthCrossCoverage,ChIPQCsample-method}
\alias{FragmentLengthCrossCoverage,list-method}
\title{Retrieve the cross coverage values when extending reads to the optimal fragment length.}
\description{
Retrieve the cross coverage values when extending reads to the optimal fragment length.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of cross-coverage values for all samples in a ChIP-seq experiment, when all reads are shifted by the optimal fragment length (the maximum cross-coverage value).
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of cross-coverage values for all samples in a list of ChIPQCsample objects, when all reads are shifted by the optimal fragment length (the maximum cross-coverage value).
}


\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve the cross-coverage value for a sample, when all reads are shifted by the optimal fragment length (the maximum cross-coverage value).
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
FragmentLengthCrossCoverage(exampleExp)
FragmentLengthCrossCoverage(QCsample(exampleExp,1))
}