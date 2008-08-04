\name{nudSegmentation}
\alias{nudSegmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Assigns each segment of each sample to either 'normal', 'up', 'down', 'amplified', or 'deleted'.}
\description{
  Classifies each segment in each sample in a SegList according to the
  absolute magnitude of the segment and the noisiness of the data for
  that sample. In doing so, converts the SegList to a ClassifySegList.
}
\usage{
nudSegmentation(seg, factor.change = 0.75, amplified.max.width = 10,
deleted.max.width = 10, amplified.magnitude = 1, deleted.magnitude = 1,
cellularity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{Object of class \code{\link[snapCGH:SegList]{SegList}}}
  \item{factor.change}{Factor of estimated noise by which a segment should deviate from zero to be classified as gained or lost.}
  \item{amplified.max.width}{Maximum width permitted for a segment to be classified as an amplicon.}
  \item{deleted.max.width}{Maximum width permitted for a segment to be classified as a deletion.}
  \item{amplified.magnitude}{Minimum (absolute) segmentation value permitted for a segment to be classified as an amplicon.}
  \item{deleted.magnitude}{Minimum (absolute) segmentation value permitted for a segment to be classified as a deletion.}  	
  \item{cellularity}{If present, a vector of the same length as the number of samples in seg, giving cellularity values as a fraction between 0-1. NAs should be used where these are not known.}
}
\value{
  Object of class \code{\link[bayesCGH]{ClassifySegList}}. This class
    extends the class \code{\link[snapCGH]{SegList}}.
}

\details{
  The function classifies each segment as either 'normal', 'up', 'down', 'amplified' or 'deleted'. A region is classified as 'up' or 'down' if the segmentation value is more (in absolute value) than 'factor.change' multiplied by the interquartile range of the difference between observed and predicted values for each observation on the genome. A region is classified as amplified or deleted if it falls within the constraints of the arguments given. It is felt that regions of amplification or deletion should be sufficiently striking not to be affected by the noise of the data.
  The 'cellularity' argument, if given, adjusts the thresholds for identifying 'up' and 'down' regions in terms of the level of cellularity given. Thus, for a sample of low observed cellularity, smaller changes should be identified as gains or losses. In cases of sufficiently low cellularity, this will inevitably lead to greater misclassification and so the classification on these samples should be checked and, if necessary, the sample discarded.
}
\seealso{
  \code{\link[snapCGH:SegList]{SegList}}
  \code{\link[bayesCGH:ClassifySegList]{ClassifySegList}}
}
\author{Thomas Hardcastle}

\keyword{manip}% at least one, from doc/KEYWORDS
