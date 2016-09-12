\name{plotSSD-methods}
\docType{methods}
\alias{plotSSD-methods}
\alias{plotSSD}
\alias{plotSSD,ChIPQCexperiment-method}
\alias{plotSSD,ChIPQCsample-method}
\alias{plotSSD,list-method}
\title{Generate SSD plot}
\description{
Generate SSD metric plot. If blacklists supplied, will generate SSD prior and post blacklisting
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", facet=T,
                                facetBy=c("Tissue","Factor"))}}{
Generate SSD metric plot for all samples in experiment. If blacklists supplied, will generate SSD prior and post blacklisting
}

\item{\code{signature(object = "list", facet=T,
                                facetBy=c("Tissue","Factor"))}}{
Generate SSD metric plot for list of samples in experiment. If blacklists supplied, will generate SSD prior and post blacklisting
}

\item{\code{signature(object = "ChIPQCsample")}}{
Generate SSD metric plot for single sample. If blacklists supplied, will generate SSD prior and post blacklisting
}
}}

\note{
plotSSD uses \code{ggplot2} for plotting, and returns a \code{ggplot2} plot gg object.
}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
plotSSD(exampleExp,facetBy=c("Tissue","Factor"))
plotSSD(QCsample(exampleExp,1))
plotSSD(QCsample(exampleExp)[1:3])
data(tamoxifen_QC)
plotSSD(tamoxifen,facetBy=c("Tissue","Condition"))
}

