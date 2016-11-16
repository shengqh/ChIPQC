\name{plotRap-methods}
\docType{methods}
\alias{plotRap-methods}
\alias{plotRap}
\alias{plotRap,ChIPQCexperiment-method}
\alias{plotRap,ChIPQCsample-method}
\alias{plotRap,list-method}
\title{Generate reads in peaks plot}
\description{
Generate reads in peaks plot.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", type="barstacked",facet=T,
                      facetBy=c("Tissue","Factor"))}}{
Generate reads in peaks plots for all the samples in an experiment.
}

\item{\code{signature(object = "list", type="barstacked",facet=T,
                      facetBy=c("Sample"))}}{
Generate reads in peaks plots for all the samples in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCexperiment", type="barstacked",facet=T,
                      facetBy=c("Tissue","Factor"))}}{
Generate reads in peaks plots for all the samples in an experiment.
}

\item{\code{signature(object = "ChIPQCsample",type="barstacked",facet=T,
                                facetBy=c("Tissue","Factor"))}}{
Generate reads in peaks plots for a sample.
}
}}

\note{
plotRap uses \code{ggplot2} for plotting, and returns a \code{ggplot2} plot dataframe.
}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
plotRap(exampleExp,facetBy=c("Tissue","Factor"))
plotRap(QCsample(exampleExp,1))
}

