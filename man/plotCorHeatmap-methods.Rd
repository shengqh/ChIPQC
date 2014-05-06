\name{plotCorHeatmap-methods}
\docType{methods}
\alias{plotCorHeatmap-methods}
\alias{plotCorHeatmap}
\alias{plotCorHeatmap,ChIPQCexperiment-method}
\title{Generate Correlation Heatmap for ChIP samples}
\description{
Generates correlation heatmap for ChIP samples.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", attributes, ...)}}{
Generate correlation heatmap, including clustering dendrogram, for all the samples in an experiment that are associated with a peakset.
\tabular{ll}{
attributes \tab character string, or vector of character strings, representing metadata field names, for use in labeling\cr
... \tab additional parameters passed to \code{\link{dba.plotHeatmap}}\cr
}
}
}
}

\note{
plotCorHeatmap uses \code{\link{dba.plotHeatmap}} for plotting.
}

\keyword{methods}
\author{Rory Stark and Thomas Carroll}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{dba.plotHeatmap}
}
\examples{
data(tamoxifen_QC)
plotCorHeatmap(tamoxifen,attributes=c("Tissue","Condition","Replicate"))
}

