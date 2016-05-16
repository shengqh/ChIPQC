\name{plotCC-methods}
\docType{methods}
\alias{plotCC-methods}
\alias{plotCC}
\alias{plotCC,ChIPQCexperiment-method}
\alias{plotCC,ChIPQCsample-method}
\title{Generate Cross-Coverage plots}
\description{
Generate Cross-Coverage plots.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", method)}}{
Generate cross-coverage plots for all the samples in an experiment. 
}

\item{\code{signature(object = "ChIPQCsample", methods)}}{
Generate cross-coverage plots for a sample.
Supported methods include:
\tabular{ll}{
"Coverage" [default] \tab Coverage plot\cr
}}
}}

\note{
plotCC uses \code{ggplot2} for plotting, and returns a \code{ggplot2} plot dataframe.
}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
plotCC(exampleExp)
plotCC(QCsample(exampleExp,1))
}

