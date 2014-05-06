\name{plotPrincomp-methods}
\docType{methods}
\alias{plotPrincomp-methods}
\alias{plotPrincomp}
\alias{plotPrincomp,ChIPQCexperiment-method}
\title{Generate Principal Components Analysis plot for ChIP samples}
\description{
Generate principal components analysis plot for ChIP samples.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", attributes, ...)}}{

Generate principal components analysis plot, for all the samples in an experiment that are associated with a peakset.
\tabular{ll}{
attributes \tab character string, or vector of character strings, representing metadata field names, for use grouping samples by color.\cr
... \tab additional parameters passed to \code{\link{dba.plotPCA}}\cr
}
}
}
}

\note{
plotPrincomp uses \code{\link{dba.plotPCA}} for plotting.
}

\keyword{methods}
\author{Rory Stark and Thomas Carroll}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{dba.plotPCA}
}
\examples{
data(tamoxifen_QC)
plotPrincomp(tamoxifen,attributes=c("Condition"))
}

