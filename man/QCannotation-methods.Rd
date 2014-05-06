\name{QCannotation-methods}
\docType{methods}
\alias{QCannotation}
\alias{QCannotation-methods}
\alias{QCannotation,ChIPQCexperiment-method}

\title{Retrieve an annotation description, or a processed annotation}
\description{
Retrieve an annotation description, or a processed annotation, from a \code{\link{ChIPQCexperiment}} object.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", bRetrieve=FALSE)}}{
Retrieve the annotation. If \code{bRetrieve=FALSE} (default), the character string describing the annotation is returned (currently only \dQuote{hg19} is supported). If \code{bRetrieve=TRUE}, a processed annotation is returned (in the form of a \code{list}). This can be used in subsequent calls to \code{\link{ChIPQC}} and/or \code{\link{ChIPQCsample}} for efficiency purposes. 
}
}}

\keyword{methods}
\author{Rory Stark and Thomas Carroll}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCsample}, \link{ChIPQCexperiment}
}
\examples{
data(example_QCexperiment)
QCannotation(exampleExp)
}
