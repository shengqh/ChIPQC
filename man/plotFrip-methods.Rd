\name{plotFrip-methods}
\docType{methods}
\alias{plotFrip-methods}
\alias{plotFrip}
\alias{plotFrip,ChIPQCexperiment-method}
\alias{plotFrip,ChIPQCsample-method}
\alias{plotFrip,list-method}
\title{Generate fraction of reads in peaks plot}
\description{
Generate fraction of reads in peaks plot.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", type="barstacked",facet=T,
                      facetBy=c("Tissue","Factor"),AsPercent=TRUE)}}{
Generate fraction of reads in peaks plots for all the samples in an experiment.
}

\item{\code{signature(object = "list", type="barstacked",facet=T,
                      facetBy=c("Sample"),AsPercent=TRUE)}}{
Generate fraction of reads in peaks plots for all the samples in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample",type="barstacked",facet=T,
                                facetBy=c("Tissue","Factor"),AsPercent=TRUE)}}{
Generate fraction of reads in peaks plots for a sample.
}
}}

\note{
plotFrip uses\code{ggplot2} for plotting, and returns a \code{ggplot2} plot dataframe.
}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
plotFrip(exampleExp)
plotFrip(QCsample(exampleExp,1))
}

