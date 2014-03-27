\name{QCcontrol-methods}
\docType{methods}
\alias{QCcontrol}
\alias{QCcontrol-methods}
\alias{QCcontrol,ChIPQCexperiment-method}
\title{Retrieve control objects associated with a sample}
\description{
Retrieve the \code{\link{ChIPQCsample}} objects representing controls, or a specific \code{\link{ChIPQCsample}} representing the control for a specific sample, from a \code{\link{ChIPQCexperiment}} object.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment", sampleID)}}{Get the control sample (as a \code{\link{ChIPQCsample}} object) associated with a ChIP sample, or, if \code{sampleID} is missing, a list of all samples used as controls. \code{sampleID} is a character string or an integer.}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
\keyword{methods}
\author{Rory Stark and Thomas Carroll}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCsample}, \link{ChIPQCexperiment}
}
\examples{
data(tamoxifen_QC)
controls = QCcontrol(tamoxifen)
length(controls)
names(controls)
controls[[1]]
bt474control = QCcontrol(tamoxifen,"BT4741")
bt474control
}