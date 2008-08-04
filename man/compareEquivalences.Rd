\name{compareEquivalences}
\alias{compareEquivalences}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian comparison of proportions of chromosomal gain and loss
  between clinical sets.}
\description{
This function compares a set of equivalence relations which define
hypotheses of similarity within and differences between equivalence
classes. The function is intended for use on array CGH data in the form
of a \code{\link[bayesCGH:ClassifySegList]{ClassifySegList}} object.}
\arguments{
  \item{seg}{Object of class \code{\link[bayesCGH:ClassifySegList]{ClassifySegList}}.}
  \item{group.arrays}{A list object of which each item should be a numeric
  vector giving the indices of the arrays in 'seg' belonging to a
  particular clinical set.}
  \item{group.names}{Names of the clinical sets defined by 'sel.groups'.}
  \item{probe.classes}{A list object of which each item is a character
    vector that describes the different groups of classifications of the
    aCGH data. Unlikely to change from default.}
  \item{classes}{A list object, of which each item is another list
    object that defines an equivalence relation on the set of sets
    defined by 'sel.groups' (see Details).}
  \item{class.names}{Names for each equivalence relation defined in 'classes'}
  \item{class.priors}{Priors for each equivalence relation defined in
    'classes' (see Details).}
  \item{X.distinct}{Should the X chromosome be omitted from estimation
    of the prior?}
  \item{Y.distinct}{Should the Y chromosome be omitted from estimation of
  the prior?}
  \item{same.priors}{Do all equivalence relations defined in 'classes'
    have the same prior Dirichlet distribution on the different
    classifications defined in 'probe.classes'?}
}
\details{
This function takes a number of subsets of the 'seg' ClassifySegList
structure and simultaneously tests a number of hypotheses defined in
terms of equivalence classes upon them.

The subsets of 'seg' to be used are defined by 'sel.groups'. This is a
list object in which the nth item is a numerical vector giving the
arrays belonging to the nth subset. Names of the subsets can be defined
in 'group.namess', a string vector.

The 'classes' argument is a list object of which each item is another
list object defining an equivalence relation on the sets defined in
'sel.groups'. Names for each equivalence relation can be defined in
'class.names', a string vector. Prior likelihoods of each equivalence
relation are defined in 'class.priors'.


}
\value{
  List object containing
  \item{data}{Description of data submitted to compareEquivalences.}
  \item{equivalence.classes}{Description of equivalence classes
    submitted to compareEquivalences.}
  \item{probabilities}{Matrix of probabilities for each probe belonging to each
    equivalence class.}
}
\references{Hardcastle 2008}
\author{Thomas Hardcastle}

\seealso{\code{\link[bayesCGH:ClassifySegList]{ClassifySegList}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

