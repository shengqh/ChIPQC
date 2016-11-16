\name{crosscoverage-methods}
\docType{methods}
\alias{crosscoverage}
\alias{crosscoverage-methods}
\alias{crosscoverage,ChIPQCexperiment-method}
\alias{crosscoverage,ChIPQCsample-method}
\alias{crosscoverage,list-method}
\title{Retrieve the cross coverage values for a range of shift sizes}
\description{
Retrieves the cross-coverage values for a range of shift sizes.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{matrix} of cross-coverage data for all samples in an ChIP-seq experiment. Each column represents a sample, and each row a shift size, with the value representing the cross-coverage using that size read.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{matrix} of cross-coverage data for all samples in a list of ChIPQCsample objects. Each column represents a sample, and each row a shift size, with the value representing the cross-coverage using that size read.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve a \code{vector} of cross-coverage data for a sample. Each position in the vector corresponds to a shift size, with the value representing the cross-coverage using that size read.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
CTCFcoverage = crosscoverage(QCsample(exampleExp,1))
length(CTCFcoverage)
plot(CTCFcoverage,type='l',
     ylab="Cross-coverage",
     xlab="Fragment length") 

allcoverages = crosscoverage(exampleExp)
dim(allcoverages)
for(i in 1:ncol(allcoverages)) lines(allcoverages[,i],col=i)
}