\name{ChIPQC-package}
\alias{ChIPQC-package}

\docType{package}
\title{
ChIPQC - Quality metrics for ChIPseq data
}
\description{
ChIPQC analyzes aligned reads (in .bam format) for ChIP-seq samples and their associated controls, computing a variety of quality control metrics and statistics, and providing reporting and plotting functions to enable assessment of experimental data for further analysis. 
}
\details{
\tabular{ll}{
Package: \tab ChIPQC\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-03-01\cr
License: \tab GPL3\cr
}

\code{ChIPQC} primarily uses two object classes: \code{\link{ChIPQCsample}}, which encapsulates the information about individual samples, and \code{\link{ChIPQCexperiment}}, which encapsulates information about larger ChIP-seq experiments (consisting of a number of samples). The primary entry point is the constructor function \code{\link{ChIPQC}}, which takes a description of an entire experiment, constructs objects for all the samples, and computes the quality metrics.

}
\author{
Tom Carroll and Rory Stark

Maintainers: Tom Carroll <tc.infomatics@gmail.com> and Rory Stark <rory.stark@cruk.cam.ac.uk>

}
\references{
Frontiers?
}

\seealso{
\code{ChIPQC} is designed to work closely with the \code{\link{DiffBind}} package, which provides functionality for analyzing ChIP-seq experiments, including performing differential binding analysis to identify significantly differentially bound peaks.
}
