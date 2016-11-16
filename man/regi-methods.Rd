\name{regi-methods}
\docType{methods}
\alias{regi}
\alias{regi-methods}
\alias{regi,ChIPQCexperiment-method}
\alias{regi,ChIPQCsample-method}
\alias{regi,list-method}
\title{Retrieve genomic profile information}
\description{
Retrieve genomic profile information in terms of relative enrichment over background genomic distribution.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ChIPQCexperiment")}}{
Retrieve a \code{matrix} of relative enrichment values for a variety of genomic features. Each column represents the enrichment values for one sample in the experiment.
}

\item{\code{signature(object = "list")}}{
Retrieve a \code{matrix} of relative enrichment values for a variety of genomic features. Each column represents the enrichment values for one sample in a list of ChIPQCsample objects.
}

\item{\code{signature(object = "ChIPQCsample")}}{
Retrieve a \code{vector} of relative enrichment values for a variety of genomic features for a sample. Relative enrichment is computed as the proportion of reads overlapping a genomic feature type compared to the overall proportion of base pairs in the genome comprising those features. Genomic features include:
\tabular{ll}{
3UTRs \tab 3' UTRs\cr
5UTRs \tab 5' UTRs\cr
Introns \tab Intronic (non-coding) portions of gene bodies\cr
Transcripts \tab Transcribed regions, including exons\cr
Promoters500 \tab 500bp regions immediately upstream of annotated TSSs\cr
Promoters2000to500 \tab 2500bp regions from 2000bp immediately upstream of annotated TSSs to 500bp downstream\cr
Promoters20000to2000 \tab 22000bp regions from 20000bp immediately upstream of annotated TSSs to 2000bp downstream\cr
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
genomicprofile = regi(exampleExp)
heatmap(genomicprofile)
regi(QCsample(exampleExp,1))
}