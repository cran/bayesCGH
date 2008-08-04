\name{selectRegions}
\alias{selectRegions}
\alias{selectAmpsDels}
\title{Selects regions with high proportions of chromosomal gains and losses.}
\description{
The 'selectRegions' function takes a 'ClassifySegList' object (or subset) and
identifies those regions on the genome that show high likelihoods of
high proportions of chromosomal gains or losses.

The 'selectAmpsDels' function simply selects those regions on the genome
that have previously been classified as amplified or deleted.
}
\usage{
selectRegions(seg, arrays, plot = FALSE, summarise = TRUE,
                  prob.prop = 0.9, prop.greater.than = 0.2, group.name =
"", ...)

selectAmpsDels(seg, arrays, group.name = "")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{Object of class \code{\link[bayesCGH:ClassifySegList]{ClassifySegList}}}
  \item{arrays}{Vector describing which arrays of 'seg' should be
    considered in region selection.}
  \item{plot}{Should the selected regions be plotted (\link[bayesCGH:plotProportions]{plotProportions})?}
  \item{summarise}{Should the selected regions be summarised (see Details).}
  \item{prob.prop}{Threshold probability; see Details.}
  \item{prop.greater.than}{Threshold proportion; see Details.}
  \item{group.name}{Name for group of arrays selected in 'arrays'.}
  \item{...}{Additional arguments to be passed to plotProportions.}
}
\details{
This function takes a 'ClassifySegList' object (or subset) and
identifies those regions on the genome that show high likelihoods of
high proportions of chromosomal gains or losses. Specifically, given a threshold
proportion 'prop.greater.than', empirical Bayes methods are applied to calculate the
posterior probability that the proportion of chromosomal gains or losses
exceeds this threshold for each probe on the arrays.

Probe selection is then made on the basis of those regions for which
this posterior probability exceeds some probability threshold
'prob.prop'.

If 'summarise = TRUE' the selected probes are summarised into continuous
regions (\link[bayesCGH:summariseRegions]{summariseRegions}).
}
\value{
  List containing:
  \item{data}{Description of data used.}
  \item{gains}{Probability for each probe that the proportion of samples
    with a chromosomal gain at that probe exceeds that given as
    'prop.greater.than'.}
  \item{losses}{Probability for each probe that the proportion of samples
    with a chromosomal loss at that probe exceeds that given as
    'prop.greater.than'.}
  \item{summary}{If 'summarise = TRUE', a summary of regions (and genes) that have a probability
    greater than 'prob.prop' of having a proportion of samples with a
    chromosomal gain or loss exceeding that of 'prop.greater.than'.}
}
\references{Hardcastle 2008}
\seealso{\code{\link[bayesCGH:nudSegmentation]{nudSegmentation}}}
\author{Thomas Hardcastle}

\keyword{manip}% at least one, from doc/KEYWORDS

