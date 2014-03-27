\name{ChIPQCsample-class}
\Rdversion{1.1}
\docType{class}
\alias{ChIPQCsample-class}
\alias{ChIPQCsample}

\alias{show,ChIPQCsample-method}

\title{Class \code{"ChIPQCsample"}}
\description{
Object containing quality metrics computed for a ChIP-seq (or associated control) sample.
}
\section{Objects from the Class}{
Objects can be created using the \code{\link{ChIPQCsample}} function.

}
\section{Slots}{
  \describe{
    \item{\code{AveragePeakSignal}:}{Object of class \code{"list"} ~~ }
    \item{\code{CrossCoverage}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{SSD}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{CountsInPeaks}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{CountsInBlackList}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{CountsInFeatures}:}{Object of class \code{"list"} ~~ }
    \item{\code{PropInFeatures}:}{Object of class \code{"list"} ~~ }
    \item{\code{CoverageHistogram}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{FlagAndTagCounts}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{readlength}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{seqnames}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{ranges}:}{Object of class \code{"IRanges"} ~~ }
    \item{\code{strand}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{elementMetadata}:}{Object of class \code{"DataFrame"} ~~ }
    \item{\code{seqinfo}:}{Object of class \code{"Seqinfo"} ~~ }
    \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{GRanges}"}
}
\section{Methods}{
  \describe{ 
    \item{averagepeaksignal}{\code{signature(object = "ChIPQCsample")}: see \code{\link{averagepeaksignal}}.}
    \item{coveragehistogram}{\code{signature(object = "ChIPQCsample")}: see \code{\link{coveragehistogram}}.}
    \item{crosscoverage}{\code{signature(object = "ChIPQCsample")}: see \code{\link{crosscoverage}}.}
    \item{flagtagcounts}{\code{signature(object = "ChIPQCsample")}: see \code{\link{flagtagcounts}}.}
    \item{fragmentlength}{\code{signature(object = "ChIPQCsample")}: see \code{\link{fragmentlength}}.}
    \item{FragmentLengthCrossCoverage}{\code{signature(object = "ChIPQCsample")}: see \code{\link{FragmentLengthCrossCoverage}}.}
    \item{frip}{\code{signature(object = "ChIPQCsample")}: see \code{\link{frip}}.}
    \item{mapped}{\code{signature(object = "ChIPQCsample")}: see \code{\link{mapped}}.}
    \item{reads}{\code{signature(object = "ChIPQCsample")}: see \code{\link{reads}}.}
    \item{duplicates}{\code{signature(object = "ChIPQCsample")}: see \code{\link{duplicates}}.}
    \item{duplicateRate}{\code{signature(object = "ChIPQCsample")}: see \code{\link{duplicateRate}}.}
    \item{Normalisedaveragepeaksignal}{\code{signature(object = "ChIPQCsample")}: see \code{\link{Normalisedaveragepeaksignal}}.}
    \item{peaks}{\code{signature(object = "ChIPQCsample")}:see \code{\link{peaks}}. }
    \item{readlength}{\code{signature(object = "ChIPQCsample")}: see \code{\link{readlength}}.}
    \item{ReadLengthCrossCoverage}{\code{signature(object = "ChIPQCsample")}: see \code{\link{ReadLengthCrossCoverage}}.}
    \item{RelativeCrossCoverage}{\code{signature(object = "ChIPQCsample")}:see \code{\link{RelativeCrossCoverage}}.}
    \item{ribl}{\code{signature(object = "ChIPQCsample")}: see \code{\link{ribl}}. }
    \item{rip}{\code{signature(object = "ChIPQCsample")}: see \code{\link{rip}}.}
    \item{show}{\code{signature(object = "ChIPQCsample")}: see \code{\link{show}}.}
    \item{ssd}{\code{signature(object = "ChIPQCsample")}: see \code{\link{ssd}}.}
    \item{regi}{\code{signature(object = "ChIPQCsample")}: see \code{\link{regi}}.}
    
    \item{plotCC}{\code{signature(object = "ChIPQCsample")}: see \code{\link{plotCC}}.}
    \item{plotCoverageHist}{\code{signature(object = "ChIPQCsample")}: see \code{\link{plotCoverageHist}}.}
    \item{plotFribl}{\code{signature(object = "ChIPQCsample")}: see \code{\link{plotFribl}}.}    
    \item{plotPeakProfile}{\code{signature(object = "ChIPQCsample")}: see \code{\link{plotPeakProfile}}.}
    \item{plotRap}{\code{signature(object = "ChIPQCsample")}: see \code{\link{plotRap}}.}   
    \item{plotRegi}{\code{signature(object = "ChIPQCsample")}: see \code{\link{plotRegi}}.}       
    }
}
\references{
Frontiers paper?
}
\author{Thomas Carroll and Rory Stark}

\seealso{

\link{ChIPQC-package}, \link{ChIPQCsample}
}
\examples{
bamFile <- system.file("extdata", "ex1.bam", 
                        package="Rsamtools")
ex1 <- ChIPQCsample(bamFile,annotation=NULL)
readlength(ex1)
fragmentlength(ex1)
}

\keyword{classes}
