\name{RelativeCrossCoverage-methods}
\docType{methods}
\alias{RelativeCrossCoverage}
\alias{RelativeCrossCoverage-methods}
\alias{RelativeCrossCoverage,ChIPQCexperiment-method}
\alias{RelativeCrossCoverage,ChIPQCsample-method}
\alias{RelativeCrossCoverage,list-method}
\title{Retrieve the relative cross coverage values for a range of shift sizes}
\description{
Retrieve the relative cross-coverage values for a range of shift sizes 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{vector} of relative cross-coverage values for all samples in a ChIP-seq experiment, computed based on the maximal value (when extending the reads to the  optimal fragment length) versus the cross-coverage values using non-extended reads.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{vector} of relative cross-coverage values for all samples in a list of ChIPQCsample objects, computed based on the maximal value (when extending the reads to the  optimal fragment length) versus the cross-coverage values using non-extended reads.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve the relative cross-coverage value for a sample, computed based on the maximal value (when extending the reads to the  optimal fragment length) versus the cross-coverage values using non-extended reads.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
RelativeCrossCoverage(exampleExp)
RelativeCrossCoverage(QCsample(exampleExp,1))
}