\name{reads-methods}
\docType{methods}
\alias{reads-methods}
\alias{reads}
\alias{reads,ChIPQCexperiment-method}
\alias{reads,ChIPQCsample-method}
\alias{reads,list-method}
\title{Retrieve numbers of reads}
\description{
Retrieve the numbers of reads.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", bFiltered)}}{
Retrieve a \code{vector} of the numbers of reads for each sample in an experiment. If \code{bFiltered=TRUE} (or is missing), this will be the number of reads that pass the mapping quality filter for each sample. if \code{bFiltered=FALSE}, it will be the total number of reads for each sample.
}

\item{\code{signature(object = "list", bFiltered)}}{
Retrieve a \code{vector} of the numbers of reads for each sample in a list of ChIPQCsample objects. If \code{bFiltered=TRUE} (or is missing), this will be the number of reads that pass the mapping quality filter for each sample. if \code{bFiltered=FALSE}, it will be the total number of reads for each sample.
}

\item{\code{signature(object = "ChIPQCsample", bFiltered)}}{
Retrieve the number of reads for a sample. If \code{bFiltered=TRUE} (or is missing), this will be the number of reads that pass the mapping quality filter. if \code{bFiltered=FALSE}, it will be the total number of reads for the sample.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
reads(exampleExp)
reads(QCsample(exampleExp,1))
}

