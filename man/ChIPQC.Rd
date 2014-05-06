\name{ChIPQC}
\alias{ChIPQC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct a \code{\link{ChIPQCexperiment}} object
}
\description{
Constructs a new \code{\link{ChIPQCexperiment}} object.
}
\usage{
ChIPQC(experiment, annotation, chromosomes, samples,
       consensus=FALSE, bCount=FALSE, mapQCth=15, blacklist=NULL, 
       profileWin=400, fragmentLength=125, shifts=1:300, ...)
}

\arguments{
  \item{experiment}{
A specification of the ChIP-seq experiment to evaluate. 
This can either be a dataframe, a filename for a .csv file, 
or a \code{\link{DBA}} object as defined in the \code{\link{DiffBind}} package.
Columns names in sample sheet may include:
\itemize{
\item {SampleID:}       {Identifier string for sample}
\item {Tissue:} {Identifier string for tissue type}
\item {Factor:} {Identifier string for factor}
\item {Condition:}      {Identifier string for condition}
\item {Treatment:}      {Identifier string for treatment}
\item {Replicate:}      {Replicate number of sample}
\item {bamReads:}       {file path for bam file containing aligned reads for ChIP sample}
\item {bamControl:}     {file path for bam file containing aligned reads for control sample}
\item {ControlID:}      {Identifier string for control sample}
\item {Peaks:}  {path for file containing peaks for sample. Format determined by PeakCaller field or caller parameter}
\item {PeakCaller:}     {Identifier string for peak caller used. If Peaks is not a bed file, this will determine how the Peaks file is parsed. If missing, will use default peak caller specified in caller parameter. Possible values:}
\itemize{
\item {\dQuote{raw}:}   {text file file; peak score is in fourth column}
\item {\dQuote{bed}:}   {.bed file; peak score is in fifth column}
\item {\dQuote{narrow}:}  {default peak.format: narrowPeaks file}
\item {\dQuote{macs}:}  {MACS .xls file}
\item {\dQuote{swembl}:}        {SWEMBL .peaks file}
\item {\dQuote{bayes}:} {bayesPeak file}
\item {\dQuote{fp4}:}   {FindPeaks v4}
}
\item {PeakFormat:}     {string indicating format for peak files; see PeakCaller and  \code{\link{dba.peakset}}}
\item {ScoreCol:}       {column in peak files that contains peak scores}
\item {LowerBetter:}    {logical indicating that lower scores signify better peaks}
}

See the documentation for the \code{sampleSheet} parameter of \code{\link{dba}} for details.
}

\item{annotation}{
Either a character string indicating the genome and version to use as a genomic annotation, or a previously defined annotation (obtained using \code{\link{QCannotation}} on a previously defined \code{\link{ChIPQCexperiment}} object.)  May be left unspecified, in which case no genomic feature analysis is performed. The following annotation specifiers are supported:
\itemize{
\item {\dQuote{hg19}:} {Human, version 19}
\item {\dQuote{hg18}:} {Human, version 18}
\item {\dQuote{mm10}:} {Mouse, version 10}
\item {\dQuote{mm9} :} {Mouse, version 19}
\item {\dQuote{rn4} :} {Rat, version 4}
\item {\dQuote{ce6} :} {C. Elgans, version 6}
\item {\dQuote{dm3} :} {D. Melanogaster, version 3}
}
Alternatively, you can construct your own annotation; see the package vignette for more information.
}

\item{chromosomes}{
Specification of which chromosomes to use for computing QC statistics. If missing, the first chromosome  which has a peak is checked. If NULL, all chromosomes will be checked (which may be time-consuming). This can be a character string (e.g. \dQuote{chr18}) or a vector or list of character strings. If it is an integer or vector of integers, the chromosomes will be checked based on the order that they are listed in a peak set.
}
\item{samples}{
\code{list} of \code{ChIPCsample} objects. If present, the sample objects will be taken directly from this list instead of being computed using the \code{\link{ChIPQCsample}} constructor.
}
\item{consensus}{
If \code{consensus} is a \code{\link{GRanges}} object, all samples will use this peakset when computing peak-based metrics. If \code{consensus=TRUE}, a consensus peakset will be generated and used for all samples, derived by merging overlapping peaks in all provided peaksets, keeping any peaks that overlap in at least two samples To avoid this behavior, set \code{consensus=FALSE}; this will result in only supplied peaksets being used for calculation of peak-based metrics (and no peak-based metric being computed for samples with no peakset specified, such as controls).
}

\item{bCount}{
if \code{TRUE}, the peak scores for all samples will be based on read counts using \code{\link{dba.count}} using a consensus peakset. If \code{consenus} is missing, any samples (such as controls) that are not already associated with a peakset will be associated with the consensus peakset (if \code{consensus} is not missing, all samples will be associated with the consensus peakset). Note that the re-counting process may be time-consuming.
}

\item{mapQCth}{
An integer representing a mapping quality score threshold. Only reads with mapping quality scores above this threshold will be used for some statistics.
}
\item{blacklist}{
A \code{\link{GRanges}} object or filename specifying a bed file containing genomic regions that should be excluded from the analysis. If missing and the \code{annotation} is \dQuote{hg19}, a default blacklist, \code{\link{blacklist_hg19}} derived from the UCSC list, will be used. No blacklist is used if this is set to NULL, or is left missing and the annotation is not \dQuote{hg19}.
}
\item{profileWin}{
An integer indicating the width, in base pairs, of the window to be used for peak profiles. Peaks will be centered on their summits, and include half the window size upstream and half downstream of this point.
}
\item{fragmentLength}{
An integer indicating the expected fragment length of the libraries. Optional, as this value will be computed for each library.
}
\item{shifts}{
A vector of values to try when computing optimal shift sizes. 
}

\item{...}{
additional parameters passed to \code{\link{dba.count}} if \code{bCount=TRUE}.
}
}
\details{
ChIPQC first constructs a new \code{\link{DBA}} object if one is not provided. Next it computes the annotation if one is not provided. The main loop constructs new \code{\link{ChIPQCsample}} objects for each sample (and unique control sample).
}
\value{
A \code{\link{ChIPQCexperiment}} object.
}

\author{Thomas Carroll and Rory Stark}
\note{
ChIPQC uses \code{\link{bplapply}} from the \code{\link{BiocParallel}} package to build the \code{\link{ChIPQCsample}} object in parallel, if supported. Control of the parallelization can be effected using \code{\link{BiocParallel}} functions, such as \code{\link{register}}.
}

\seealso{
\link{ChIPQC-package},  \link{ChIPQCexperiment}, \link{ChIPQCsample}, \link{DiffBind}
}


\examples{
## Not run: exampleExp = ChIPQC(samples,annotation="hg19")
data(example_QCexperiment)
exampleExp

## Not run: tamoxifen = ChIPQC(samples, ,annotation="hg19", consensus=TRUE, bCounts=T)
data(tamoxifen_QC)
tamoxifen
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

