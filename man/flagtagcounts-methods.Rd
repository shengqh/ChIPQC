\name{flagtagcounts-methods}
\docType{methods}
\alias{flagtagcounts}
\alias{flagtagcounts-methods}
\alias{flagtagcounts,ChIPQCexperiment-method}
\alias{flagtagcounts,ChIPQCsample-method}
\alias{flagtagcounts,list-method}
\title{Retrieve numbers of reads that pass various filters}
\description{
Retrieve numbers of reads that pass various filters
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{matrix} of counts passing various filters for all the samples in an experiment. Each column represents the counts in a sample.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{matrix} of counts passing various filters for all the samples in a list of ChIPQCsample objects. Each column represents the counts in a sample.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve a \code{vector} of counts passing various filters. The values are:
\tabular{ll}{
UnMapped \tab Number of reads that are not mapped (aligned)\cr
Mapped \tab Number of reads that are mapped (aligned)\cr
Duplicates \tab Number of reads that align to exactly the same place as another read\cr
MapQPass \tab Number of reads with a mapping quality score greater than or equal to the specified threshold\cr
MapQPassandDup \tab Number of reads that are mapped (aligned) and not duplicates\cr
}
}
}}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
flagtagcounts(exampleExp)
flagtagcounts(QCsample(exampleExp,1))
}