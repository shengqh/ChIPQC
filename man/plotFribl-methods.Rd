\name{plotFribl-methods}
\docType{methods}
\alias{plotFribl-methods}
\alias{plotFribl}
\alias{plotFribl,ChIPQCexperiment-method}
\alias{plotFribl,ChIPQCsample-method}
\alias{plotFribl,list-method}
\title{Generate fraction of reads in blacklist plot}
\description{
Generate fraction of reads in blacklist plot.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", type="barstacked",facet=T,
                     facetBy=c("Tissue","Factor"),AsPercent=TRUE)}}{
Generate fraction of reads in blacklist plots for all the samples in an experiment.
}

\item{\code{signature(object = "list", type="barstacked",facet=T,
                     facetBy=c("Sample"),AsPercent=TRUE)}}{
Generate fraction of reads in blacklist plots for all the samples in a list of ChIPQCsample objects..
}

\item{\code{signature(object = "ChIPQCsample", type="barstacked",AsPercent=TRUE)}}{
Generate fraction of reads in blacklist plots for a sample.
}
}}

\note{
plotFribl uses \code{ggplot2} for plotting, and returns a \code{ggplot2} plot dataframe.
}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
plotFribl(exampleExp)
plotFribl(QCsample(exampleExp,1))
}

