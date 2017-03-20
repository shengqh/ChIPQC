\name{QCmetadata-methods}
\docType{methods}
\alias{QCmetadata}
\alias{QCmetadata-methods}
\alias{QCmetadata,ChIPQCexperiment-method}
\alias{QCmetadata,list-method}
\title{Retrieve metadata associated with an experiment}
\description{
Retrieve metadata for a \code{\link{ChIPQCexperiment}} object.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a data frame containing metadata for all the samples in a ChIP-seq experiment represented by a \code{\link{ChIPQCexperiment}} object.
}

\item{\code{signature(object = "list")}}{
Retrieve a data frame containing metadata for all the samples in a list of ChIPQCsample objects.
}
}}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
\keyword{methods}
\author{Rory Stark and Thomas Carroll}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCsample}, \link{ChIPQCexperiment}
}
\examples{
data(tamoxifen_QC)
meta = QCmetadata(tamoxifen)
meta
}