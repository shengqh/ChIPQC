\name{plotCoverageHist-methods}
\docType{methods}
\alias{plotCoverageHist-methods}
\alias{plotCoverageHist}
\alias{plotCoverageHist,ChIPQCexperiment-method}
\alias{plotCoverageHist,ChIPQCsample-method}
\alias{plotCoverageHist,list-method}
\title{Generate coverage histogram plot}
\description{
Generate coverage histogram plot.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Generate coverage histogram plots for all the samples in an experiment. 
}

\item{\code{signature(object = "list")}}{
Generate coverage histogram plots for all the samples in a list of ChIPQCsamples. 
}

\item{\code{signature(object = "ChIPQCsample")}}{
Generate coverage histogram plots for a sample.
}
}}

\note{
Uses \code{ggplot2} for plotting, and returns a \code{ggplot2} plot dataframe.
}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
plotCoverageHist(exampleExp,facetBy=c("Tissue","Factor"))
plotCoverageHist(QCsample(exampleExp,1))
}

