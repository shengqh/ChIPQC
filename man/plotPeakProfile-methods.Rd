\name{plotPeakProfile-methods}
\docType{methods}
\alias{plotPeakProfile-methods}
\alias{plotPeakProfile}
\alias{plotPeakProfile,ChIPQCexperiment-method}
\alias{plotPeakProfile,ChIPQCsample-method}
\alias{plotPeakProfile,list-method}
\title{Generate peak profile plots}
\description{
Generate peak profile plots.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", method)}}{
Generate peak profile plots for all the samples in an experiment. 
}

\item{\code{signature(object = "list", method)}}{
Generate peak profile plots for all the samples in a list of ChIPQCsample objects.. 
}
\item{\code{signature(object = "ChIPQCsample", method)}}{
Generate peak profile plots for a sample.
}
}}

\note{
plotPeakProfile uses \code{ggplot2} for plotting, and returns a \code{ggplot2} plot dataframe.
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
data(tamoxifen_QC)
plotCoverageHist(tamoxifen,facetBy=c("Tissue","Condition"))
}

