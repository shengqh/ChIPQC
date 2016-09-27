\name{duplicateRate-methods}
\docType{methods}
\alias{duplicateRate-methods}
\alias{duplicateRate}
\alias{duplicateRate,ChIPQCexperiment-method}
\alias{duplicateRate,ChIPQCsample-method}
\alias{duplicateRate,list-method}
\title{Retrieve duplication rates}
\description{
Retrieve duplication rates.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", bFiltered)}}{
Retrieve a \code{vector} of the duplication rates for each sample in an experiment. A read is considered duplicated if another read maps to the same genomic location; the duplication rate is the number of duplicated reads divided by the total number of reads for a sample. If \code{bFiltered=TRUE} (or is missing), only reads that pass the mapping quality filter for each sample are included. if \code{bFiltered=FALSE}, all reads for each sample will be included.
}

\item{\code{signature(object = "list", bFiltered)}}{
Retrieve a \code{vector} of the duplication rates for each sample in a list of ChIPQCsample objects. A read is considered duplicated if another read maps to the same genomic location; the duplication rate is the number of duplicated reads divided by the total number of reads for a sample. If \code{bFiltered=TRUE} (or is missing), only reads that pass the mapping quality filter for each sample are included. if \code{bFiltered=FALSE}, all reads for each sample will be included.
}

\item{\code{signature(object = "ChIPQCsample", bFiltered)}}{
Retrieve the duplication rate for a sample. A read is considered duplicated if another read maps to the same genomic location; the duplication rate is the number of duplicated reads divided by the total number of reads for the sample. If \code{bFiltered=TRUE} (or is missing), only reads that pass the mapping quality filter for the sample are included. if \code{bFiltered=FALSE}, all reads for the sample will be included.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
duplicateRate(exampleExp)
duplicateRate(QCsample(exampleExp,1))
}

