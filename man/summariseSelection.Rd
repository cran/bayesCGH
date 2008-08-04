\name{summarySelection}
\alias{summarySelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produce summary of selected probes.}
\description{
When given a 'ClassifySegList' object, a vector of probes selected from
this object, and a set of sets of arrays, this function summarises the
probes into a set of regions and gives data on proportions of
classification at each region.}
\usage{
summarySelection(seg, regions.clone.no, stats, group.arrays, group.names,
probe.classes = list(c("Up", "Amplified"), c("Down", "Deleted"), c("Normal")),
probe.class.name = c("Normal", "Gain", "Loss"), local.maximum = FALSE,
cytobands = cytobands)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{Object of class \code{\link[bayesCGH:ClassifySegList]{ClassifySegList}}.}
  \item{regions.clone.no}{Vector giving the probes to be summarised.}
  \item{stats}{Vector of some summary statistic (for all probes).}
  \item{group.arrays}{A list object of which each item should be a numeric
  vector giving the indices of the arrays in 'seg' belonging to a
  particular clinical set.}
  \item{group.names}{Names of the clinical sets defined by 'sel.groups'.}
  \item{probe.classes}{A list object of which each item is a character
    vector that describes the different groups of classifications of the
    aCGH data. Unlikely to change from default.}
    \item{probe.class.name}{A vector of names for the sets of probe
      classifications defined in 'probe.classes'. Again, unlikely to
      change from default.}
    \item{local.maximum}{If TRUE, then only those regions with locally
      maximal summary statistic will be reported.}
    \item{cytobands}{Location of cytobands on genome. Given in data
      attached to bayesCGH package.}
}
\details{
Probes from an aCGH array that are adjacent on the genome will often
have exactly the same patterns of gain and loss across multiple arrays,
and  in the case of dense arrays, will come from the same gene. This
function is intended to summarise the selected probes into contigious
regions with identical characteristics.
}
\value{
  List containing
  \item{summary}{Dataframe summarising the probes selected.}
  \item{genes}{Details of all genes associated with the selected probes.}
}
\author{Thomas Hardcastle}
}
\seealso{\code{\link[bayesCGH:writeSummaries]{writeSummaries}}.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

