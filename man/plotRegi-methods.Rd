\name{plotRegi-methods}
\docType{methods}
\alias{plotRegi-methods}
\alias{plotRegi}
\alias{plotRegi,ChIPQCexperiment-method}
\alias{plotRegi,ChIPQCsample-method}
\alias{plotRegi,list-method}
\title{Generate relative enrichment of genomic features plot}
\description{
Generate relative enrichment of genomic features plot.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", facet=T,
                                facetBy=c("Tissue","Factor"))}}{
Generate relative enrichment of genomic features plots for all the samples in an experiment.
}

\item{\code{signature(object = "list",facet=T,
                                facetBy=c("Sample"))}}{
Generate relative enrichment of genomic features plots for all the samples in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Generate relative enrichment of genomic features plots for a sample.
}
}}

\note{
plotRegi uses \code{ggplot2} for plotting, and returns a \code{ggplot2} plot dataframe.
}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
plotRegi(exampleExp,facetBy=c("Tissue","Factor"))
plotRegi(QCsample(exampleExp,1))
data(tamoxifen_QC)
plotRegi(tamoxifen,facetBy=c("Tissue","Condition"))
}

