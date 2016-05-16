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

\section{Constructor Function}{
\code{
ChIPQCsample(reads, peaks, annotation, chromosomes = NULL, mapQCth = 15, 
             blacklist, profileWin = 400, fragmentLength = 125, shifts = 1:300, 
             runCrossCor = FALSE,verboseT=FALSE)
}
\itemize{
  \item{\code{reads}}{character string filename of .bam file}
  \item{\code{peaks}}{\code{\link{GRanges}} object or character string filename of peaks. If present, peak-based metrics will be computed.
}

\item{\code{annotation}}{
Either a character string indicating the genome and version to use as a genomic annotation, or a previously defined annotation (obtained using \code{\link{QCannotation}} on a previously defined \code{\link{ChIPQCexperiment}} object.)  May be left unspecified, in which case no genomic feature analysis is performed. The following annotation specifiers are supported:

\tabular{ll}{
"hg20" \tab Human, version 20 \cr
"hg19" \tab Human, version 19 \cr
"hg18" \tab Human, version 18 \cr
"mm10" \tab Mouse, version 10 \cr
"mm9" \tab Mouse, version 19 \cr
"rn4" \tab Rat, version 4 \cr
"ce6" \tab C. Elgans, version 6 \cr
"dm3"  \tab D. Melanogaster, version 3 \cr
}

Alternatively, you can construct your own annotation; 
see the package Vignette for more information.
}

\item{\code{chromosomes}}{
Specification of which chromosomes to use for computing QC statistics. If missing, the first chromosome which has a peak is checked. If NULL, all chromosomes will be checked (which may be time-consuming). This can be a character string (e.g. \dQuote{chr18}) or a vector or list of character strings. If it is an integer or vector of integers, the chromosomes will be checked based on the order that they are listed in a peak set.
}
  \item{\code{mapQCth}}{
An integer representing a mapping quality score threshold. Only reads with mapping quality scores above this threshold will be used for some statistics.
}
  \item{\code{blacklist}}{
A \code{\link{GRanges}} object or filename specifying a bed file containing genomic regions that should be excluded from the analysis. If missing and the \code{annotation} is \dQuote{hg19}, a default blacklist, \code{\link{blacklist_hg19}} derived from the UCSC list, will be used. No blacklist is used if this is set to NULL, or is left missing and the annotation is not \dQuote{hg19}.
}
  \item{\code{profileWin}}{
An integer indicating the width, in base pairs, of the window to be used for peak profiles. Peaks will be centered on their summits, and include half the window size upstream and half downstream of this point.
}
  \item{\code{fragmentLength}}{
An integer indicating the expected fragment length of the libraries. Optional, as this value will be computed.
}
  \item{\code{shifts}}{
A vector of values to try when computing optimal shift sizes. 
}
  \item{\code{runCrossCor}}{
Compute cross-correlation in addition to cross-coverage. This will take more compute time, and is currently not used in the final report.
}
  \item{\code{verboseT}}{
  TRUE or FALSE, specifying whether to report progress. Default is TRUE. When set to FALSE ChIPQC does not report any progress until complete.
}
}
}

\section{Slots}{
  \describe{
    \item{\code{AveragePeakSignal}:}{Object of class \code{"list"} }
    \item{\code{CrossCoverage}:}{Object of class \code{"numeric"}  }
    \item{\code{CrossCorrelation}:}{Object of class \code{"numeric"}  }
    \item{\code{SSD}:}{Object of class \code{"numeric"}  }
    \item{\code{SSDBL}:}{Object of class \code{"numeric"}  }
    \item{\code{CountsInPeaks}:}{Object of class \code{"numeric"} }
    \item{\code{CountsInBlackList}:}{Object of class \code{"numeric"} }
    \item{\code{CountsInFeatures}:}{Object of class \code{"list"}  }
    \item{\code{PropInFeatures}:}{Object of class \code{"list"}  }
    \item{\code{CoverageHistogram}:}{Object of class \code{"numeric"} }
    \item{\code{FlagAndTagCounts}:}{Object of class \code{"numeric"}  }
    \item{\code{readlength}:}{Object of class \code{"numeric"}  }
    \item{\code{seqnames}:}{Object of class \code{"Rle"}  }
    \item{\code{ranges}:}{Object of class \code{"IRanges"}  }
    \item{\code{strand}:}{Object of class \code{"Rle"}  }
    \item{\code{elementMetadata}:}{Object of class \code{"DataFrame"}  }
    \item{\code{seqinfo}:}{Object of class \code{"Seqinfo"} }
    \item{\code{metadata}:}{Object of class \code{"list"} }
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
Carroll TS, Liang Z, Salama R, Stark R and Santiago Id (in press). Impact of artefact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data. Frontiers in Genetics.
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