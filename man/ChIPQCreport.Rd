\name{ChIPQCreport-methods}
\docType{methods}
\alias{ChIPQCreport-methods}
\alias{ChIPQCreport}
\alias{ChIPQCreport,ChIPQCexperiment-method}
\alias{ChIPQCreport,ChIPQCsample-method}
\title{Generate a summary QC report}
\description{
Generate an HTML summary QC report.
}
\section{Methods}{
\describe{

\item{\code{signature(experiment = "ChIPQCexperiment", facet=TRUE, 
                      reportName="ChIPQC", reportFolder="ChIPQCreport",
                      facetBy=c("Tissue","Factor"),colourBy=c("Replicate"))}}{

Generates an summary QC report for the experiment in HTML format.

\tabular{ll}{
experiment \tab \code{\link{ChIPQCexperiment}} object\cr
facet \tab logical indicating whether or not to facet using experimental metadata.\cr
reportName \tab filename of main report file (.html).\cr
reportFolder \tab directory name where plot graphics will be saved\cr
facetBy \tab metadata fileds to use for faceting\cr
colourBy \tab metadata feild to color by\cr
}
}

\item{\code{signature(sample = "ChIPQCsample", 
                      reportName="ChIPQC", reportFolder="ChIPQCreport",)}}{
 
Generate a summary QC report for a sample in HTML format.
 
\tabular{ll}{
sample \tab \code{\link{ChIPQCsample}} object\cr
reportName \tab filename of main report file (.html).\cr
reportFolder \tab directory name where plot graphics will be saved\cr
}
}
}}

\note{
ChIPQCreport uses \code{Nozzle.R2} for generating HTML.
}

\keyword{methods}
\author{Thomas Carroll and Rory Stark}

\seealso{
\link{ChIPQC-package}, \link{ChIPQCexperiment}, \link{ChIPQCsample}
}
\examples{
data(example_QCexperiment)
ChIPQCreport(exampleExp,facetBy=c("Tissue","Factor"))
#report in ChIPQCreport/Example.html

data(tamoxifen_QC)
ChIPQCreport(tamoxifen,facetBy="Tissue")
#report in ChIPQCreport/ChIPQC.html

ChIPQCreport(tamoxifen,facetBy=c("Tissue","Condition"))
#report in ChIPQCreport/ChIPQC.html
}

