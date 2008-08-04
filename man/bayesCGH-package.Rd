\name{bayesCGH-package}
\alias{bayesCGH-package}
\alias{bayesCGH}
\docType{package}
\title{bayesCGH - Empirical bayesian analysis of aCGH data}
\description{
This package is intended to take data analysed by the 'snapCGH' package
and, by means of the 'compareEquivalences' function, analyse the data in
terms of a set of biological hypotheses established on the data. The
output of this analysis are the posterior probabilities of each
hypothesis for each region of the genome.
}
\details{
\tabular{ll}{
Package: \tab bayesCGH\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-04-30\cr
License: \tab GPL-3\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Thomas Hardcastle

Maintainer: Thomas Hardcastle (tjh48@cam.ac.uk)
}
\references{
Empirical Bayesian methods for model testing on array CGH data. Thomas
Hardcastle, Maria J. Garcia, James D. Brenton, Simon Tavare. In preparation.
}
\keyword{ package }

\examples{
#########################
## not run
# library(bayesCGH)
#
## load segInfo data from 'snapCGH' package.
# load("segInfo.DNACopy.merged.Rdata")
# seg <- nudSegmentation(segInfo.DNACopy.merged)
#
## Define sets of non-intersecting clinical sets on samples
# sel.groups <- list()
# sel.groups[[1]] <- which(seg$samples$Treatment == "T0" & seg$samples$CA125.response == "Responder" & seg$samples$Rejected == FALSE)
# sel.groups[[2]] <- which(seg$samples$Treatment == "T0" & seg$samples$CA125.response == "Non-responder" & seg$samples$Rejected == FALSE)
# sel.groups[[3]] <- which(seg$samples$Treatment == "T1" & seg$samples$CA125.response == "Responder" & seg$samples$Rejected == FALSE)
# sel.groups[[4]] <- which(seg$samples$Treatment == "T1" & seg$samples$CA125.response == "Non-responder" & seg$samples$Rejected == FALSE)
#
## Form hypotheses on data
# classes <- list()
# classes[[1]] <- list(c(1,2,3,4))
# classes[[2]] <- list(c(1,2), c(3,4))
# classes[[3]] <- list(c(1,3), c(2,4))
# classes[[4]] <- list(c(1), c(2, 3,4))
# class.names <- c("conserved effects", "treatment effects", "response effects", "selected response effects")
# class.prior <- rep(1/4, length(classes))
#
## Find posterior probabilities
# compeq <- compareEquivalences(seg, sel.groups, c("pre-sensitive",  "pre-resistant",  "post-sensitive", "post-resistant"), probe.classes = list(c("Normal"), c("Down", "Deleted"), c("Up", "Amplified")), class.priors = class.prior, classes = classes, class.names = class.names)
#
## Summarising regions associated with one hypothesis
# conserved <- summarySelection(seg, which(compeq$probabilities[1,] > 0.5), carbo.compeq$probabilities[1,], carbo.sel.groups,
#                 c("pre-sensitive",  "pre-resistant", "post-sensitive", "post-resistant"))
## Plot posterior probabilities for each hypothesis
# for(eq.class in 1:4)
#   plot.up.down(cbind(compeq$probabilities[eq.class,], 0), seg$genes, main = class.names[eq.class], ylim = c(0,1), num.chr = 24)
}
