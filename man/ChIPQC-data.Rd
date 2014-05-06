\name{ChIPQC-data}
\alias{tamoxifen_QC}
\alias{example_QCexperiment}
\alias{blacklist_hg19}
\docType{data}
\title{
Example data sets for \code{\link{ChIPQC}} package, each containing a \code{\link{ChIPQCexperiment}} object, as well as a pre-compiled blacklist for hg19.
}
\description{
The \code{tamoxifen_QC} example data set contains a \code{\link{ChIPQCexperiment}} object named \code{tamoxifen}. This data set, based on Ross-Innes et al (2012), includes 11 ER ChIP libraries, and their input controls, divided into tamoxifen responsive and tamoxifen resistant sample groups. Only data for chromosome 18 (chr18) are included.

The \code{example_QCexperiment}  data set contains a \code{\link{ChIPQCexperiment}} object named \code{exampleExp}. This data set, derived from ENCODE data, includes 6 ChIP libraries. Only data for chromosome 22 (chr22) are included.

\code{blacklist_hg19} data set includes a \code{\link{GRanges}} object named \code{blacklist_hg19} containing blacklisted regions for the human genome, derived from the UCSC blacklist.
}
\usage{
data(tamoxifen_QC)
data(example_QCexperiment)
data(blacklist_hg19)
}
\format{
\code{tamoxifen_QC}: A single \code{\link{ChIPQCexperiment}} object named \code{tamoxifen} is loaded. This object is used for the \code{\link{ChIPQC-package}} examples and vignette. This object can also be used with the \code{\link{DiffBind}} package (see related data objects \code{\link{tamoxifen}}).

\code{example_QCexperiment}: A single \code{\link{ChIPQCexperiment}} object named \code{exampleExp} is loaded. This object is used for the \code{\link{ChIPQC-package}} examples and vignette.

\code{blacklist_hg19}: A single \code{\link{GRanges}} object named \code{blacklist.hg19} that is used by default when processing hg19 data sets.
}

\source{
  Ross-Innes, C. S., Stark, R., Teschendorff, A. E., Holmes, K. A.,
  Ali, H. R., Dunning, M. J., Brown, G. D., Gojis, O., Ellis, I. O.,
  Green, A. R., Ali, S., Chin, S.-F., Palmieri, C., Caldas, C., and
  Carroll, J. S. (2012). Differential oestrogen receptor binding is
  associated with clinical outcome in breast cancer. Nature 481,
  389-393.
}

\examples{
data(tamoxifen_QC)
tamoxifen
plotRegi(tamoxifen)

data(example_QCexperiment)
exampleExp

data(blacklist_hg19)
blacklist.hg19
}
\keyword{datasets}
