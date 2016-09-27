\name{averagepeaksignal-methods}
\docType{methods}
\alias{averagepeaksignal}
\alias{averagepeaksignal-methods}
\alias{averagepeaksignal,ChIPQCexperiment-method}
\alias{averagepeaksignal,ChIPQCsample-method}
\alias{averagepeaksignal,list-method}
\title{Retrieve average peak profiles}
\description{
Retrieve the average peak profile for a sample or set of samples.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{matrix} of the average peak profiles for all of the samples in an ChIP-seq experiment. Each column represents a sample, and each row a base pair position, centered on peak summits. 
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{matrix} of the average peak profiles for all of the samples in a list of ChIPQC samples. Each column represents a sample, and each row a base pair position, centered on peak summits. 
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve a \code{vector} representing the average peak profile for a sample. Each column represents a basepair position, centered on the peak summits.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
CTCFprofile = averagepeaksignal(QCsample(exampleExp,1))
length(CTCFprofile)
plot(CTCFprofile,type='l',ylab="mean pileup")

allprofiles = averagepeaksignal(exampleExp)
dim(allprofiles)
for(i in 1:ncol(allprofiles)) lines(allprofiles[,i],col=i)

}