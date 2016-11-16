\name{coveragehistogram-methods}
\docType{methods}
\alias{coveragehistogram}
\alias{coveragehistogram-methods}
\alias{coveragehistogram,ChIPQCexperiment-method}
\alias{coveragehistogram,ChIPQCsample-method}
\alias{coveragehistogram,list-method}
\title{Retrieve histogram data representing densities of coverage pileups}
\description{
Retrieve histogram data representing densities of coverage pileups.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{matrix} of coverage histogram data for all samples in a ChIP-seq experiment. Each column represents a sample, and each row a pileup height, with the value representing the number of basepair positions that report this pileup height.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{matrix} of coverage histogram data for all ChIPQCsamples in a list. Each column represents a sample, and each row a pileup height, with the value representing the number of basepair positions that report this pileup height.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve a \code{vector} representing coverage histogram data for a sample. Values represent the number of base pairs positions that report the pileup value. The value in position 1 of the vector the contains the number of examined basepair positions that are overlapped by exactly zero reads, while position 2 shows the number of basepair positions overlapped by exactly one read, etc.
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
CTCFcoverage = coveragehistogram(QCsample(exampleExp,1))
length(CTCFcoverage)
plot(log10(CTCFcoverage),type='l',ylab="log10 Density",
                          xlab="Pileup")

allcoverages = coveragehistogram(exampleExp)
dim(allcoverages)
for(i in 1:ncol(allcoverages)) lines(log10(allcoverages[,i]),col=i)
}
