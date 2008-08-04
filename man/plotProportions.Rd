\name{plotProportions}
\alias{plotProportions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots proportions of gains and losses observed in a set of array CGH
  data.}
\description{
  Plots proportions of gains and losses observed in a set of array CGH
  data held in a \code{\link[snapCGH:ClassifySegList]{ClassifySegList}}
  object.  Proportions of gain are plotted on the positive axis,
  proportions of loss are plotted on the negative axis. Regions which
  are frequently gained or lost in the sample set can thus be visually identified.
}
\usage{
plotProportions(seg, arrays, chrom.info = chrominfo.Mb, ylimit=c(-1,1),
                    main="Frequency of genetic changes across tumours", num.chr = 24,
                    amps.dels = TRUE, errors = TRUE, select = FALSE,
                    colors = c("green", "red", "dark green", "dark red"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{Object of class
    \code{\link[snapCGH:ClassifySegList]{ClassifySegList}}.}
  \item{arrays}{Which arrays to plot}
  \item{chrom.info}{}
  \item{ylimit}{Limits on Y-axis. Typically these will be c(-1,1) as
    this is the maximum range of proportion of loss and gain at any
    region.}
  \item{main}{Title for the plot}
  \item{num.chr}{Number of chromosomes to plot; the X and Y chromosomes
    may be omitted by setting 'num.chr = 22'}
  \item{amps.dels}{Plot proportions of amplifications and deletions in
    addition to proportions of gain and loss.}
  \item{errors}{Calculate and plot credibility intervals around
    proportions. If TRUE, 'amps.dels' is ignored.}
  \item{select}{If given, should be a vector identifying the probes to plot}
  \item{colors}{The four colors to be used in the plot.}
}
\details{
  This function plots, for each region along the genome, the proportions
  of samples in 'seg' showing chromosomal gain (positive) or loss (negative) at
  that point.

  If 'errors = TRUE' then 50\% and 90\% credibility intervals for the proportions of
  gain and loss are calculated using an empirical Bayes method with a
  Dirichlet prior (Hardcastle 2008).
}

\references{Hardcastle 2008}
\author{Thomas Hardcastle}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}