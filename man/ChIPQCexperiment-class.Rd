\name{ChIPQCexperiment-class}
\Rdversion{1.1}
\docType{class}
\alias{ChIPQCexperiment-class}
\alias{ChIPQCexperiment}

\alias{show,ChIPQCexperiment-method}

\title{ChIPQCexperiment instances}

\description{
The \code{ChIPQCexperiment} class is built around a list of \code{\link{ChIPQCsample}} objects, each representing a ChIP or control sample in a ChIP-seq experiment. These objects are created using the \code{\link{ChIPQC}} function.
}

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}: internal}
    \item{\code{Samples}:}{Object of class \code{"list"}: List of \code{\link{ChIPQCsample}} objects.}
    \item{\code{DBA}:}{Object of class \code{"DBA"}: \code{\link{DBA}} object (from package \code{\link{DiffBind}})}
    \item{\code{annotation}:}{Object of class \code{"list"} : annotation data}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}.
}
\section{Methods}{
  \describe{

    \item{QCmetadata}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{QCmetadata}}.}
\item{QCmetrics}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{QCmetrics}}.}
    \item{QCsample}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{QCsample}}.}
    \item{QCcontrol}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{QCcontrol}}.}
    \item{QCannotation}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{QCannotation}}.}
    \item{QCdba}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{QCdba}}.}
  
  
    \item{averagepeaksignal}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{averagepeaksignal}}.}
    \item{coveragehistogram}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{coveragehistogram}}.}
    \item{crosscoverage}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{crosscoverage}}.}
    \item{flagtagcounts}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{flagtagcounts}}.}
    \item{fragmentlength}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{fragmentlength}}.}
    \item{FragmentLengthCrossCoverage}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{FragmentLengthCrossCoverage}}.}
    \item{frip}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{frip}}.}
    \item{mapped}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{mapped}}.}
    \item{reads}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{reads}}.}
    \item{duplicates}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{duplicates}}.}
    \item{duplicateRate}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{duplicateRate}}.}
    \item{Normalisedaveragepeaksignal}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{Normalisedaveragepeaksignal}}.}
    \item{peaks}{\code{signature(object = "ChIPQCexperiment")}:see \code{\link{peaks}}. }
    \item{readlength}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{readlength}}.}
    \item{ReadLengthCrossCoverage}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{ReadLengthCrossCoverage}}.}
    \item{RelativeCrossCoverage}{\code{signature(object = "ChIPQCexperiment")}:see \code{\link{RelativeCrossCoverage}}.}
    \item{ribl}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{ribl}}. }
    \item{rip}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{rip}}.}
    \item{show}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{show}}.}
    \item{ssd}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{ssd}}.}
    \item{regi}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{regi}}.}

    \item{plotCC}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{plotCC}}.}
    \item{plotCoverageHist}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{plotCoverageHist}}.}
    \item{plotFribl}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{plotFribl}}.}    
    \item{plotPeakProfile}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{plotPeakProfile}}.}
    \item{plotRap}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{plotRap}}.}   
    \item{plotRegi}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{plotRegi}}.}   
    \item{plotCorHeatmap}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{plotCorHeatmap}}.}   
    \item{plotPrincomp}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{plotPrincomp}}.}   
        
        \item{ChIPQCreport}{\code{signature(object = "ChIPQCexperiment")}: see \code{\link{ChIPQCreport}}.} 
    }
}

%%\references{
%%Frontiers paper?
%%}

\author{Thomas Carroll and Rory Stark}


\seealso{

\link{ChIPQC-package}, \link{ChIPQCsample}, \link{DiffBind}
}
\examples{
## Not run: exampleExp = ChIPQC(samples)
data(example_QCexperiment)
exampleExp

## Not run: tamoxifen = ChIPQC(samples, consensus=TRUE, bCounts=T)
data(tamoxifen_QC)
tamoxifen
}

\keyword{classes}
