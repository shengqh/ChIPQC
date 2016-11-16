\name{Normalisedaveragepeaksignal-methods}
\docType{methods}
\alias{Normalisedaveragepeaksignal}
\alias{Normalisedaveragepeaksignal-methods}
\alias{Normalisedaveragepeaksignal,ChIPQCexperiment-method}
\alias{Normalisedaveragepeaksignal,ChIPQCsample-method}
\alias{Normalisedaveragepeaksignal,list-method}
\title{Retrieve normalised average peak profiles}
\description{
Retrieve normalised average peak profiles
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{matrix} of normalised average peak signal data for all samples in a ChIP-seq experiment. Each column represents a sample, and each row a base pair position, centered on peak summits. 
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{matrix} of normalised average peak signal data for all samples in a list of ChIPQCsample objects. Each column represents a sample, and each row a base pair position, centered on peak summits. 
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve a \code{vector} representing the normalised average peak profile for a sample. Each column represents a basepair position, centered on the peak summits.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
CTCFprofile = Normalisedaveragepeaksignal(QCsample(exampleExp,1))
length(CTCFprofile)
plot(CTCFprofile,type='l',ylab="normalised mean pileup")

allprofiles = Normalisedaveragepeaksignal(exampleExp)
dim(allprofiles)
for(i in 1:ncol(allprofiles)) lines(allprofiles[,i],col=i)
}