\name{sic.bayes}
\alias{sic.bayes}
\title{Hierarchical posterior survivor interaction contrast (SIC)}
\description{
  Extracts the posterior survivor interaction contrast \eqn{\mathrm{SIC}(t) =
  S_{LL}(t) - S_{LH}(t) - S_{HL}(t) + S_{HH}(t)} from an OR-centred hierarchical
  Bayesian capacity model fitted with a salience split. The curve is returned at
  the subject (partially pooled), population (typical-subject parameter),
  \code{population_finite_mean}, and \code{new_subject} (posterior-predictive)
  levels, each with pointwise posterior sign probabilities. It is only
  identified from a salience-split fit.
}
\usage{
sic.bayes(object, ...)
}
\arguments{
  \item{object}{A fitted \code{sft_bayes} object, or a canonical SFT
    data frame to fit with \code{\link{semiparametricSFT.bayes}}.}
  \item{...}{Passed to \code{\link{semiparametricSFT.bayes}} when \code{object} is a data
    frame. A \code{salience_split} is required to identify the four cells and
    defaults to \code{TRUE}.}
}
\value{
  A list with a tidy \code{summary} data frame (columns include \code{Level},
  \code{Subject}, \code{Time}, posterior \code{Mean} and highest-density
  interval, and \code{Prob_super} \eqn{= P(\mathrm{SIC}(t) > 0)} and
  \code{Prob_limited} \eqn{= P(\mathrm{SIC}(t) < 0)}), per-level views
  (\code{subject}, \code{population}, \code{new_subject}), and metadata.
}
\details{
  Unlike a fit-each-subject-then-average approach, the subject curves are
  partially pooled toward the population, the \code{population} level is a model
  parameter (the typical subject, with random effects at zero), and
  \code{new_subject} is the posterior-predictive SIC for an unobserved
  participant, which carries the between-subject variance. The classic SIC
  signatures still apply: an all-positive curve indicates parallel
  self-terminating or coactive processing, an all-negative curve indicates
  parallel exhaustive processing, and a negative-then-positive S shape with zero
  net area indicates serial processing.
}
\seealso{\code{\link{semiparametricSFT.bayes}}, \code{\link{mic.bayes}},
  \code{\link{sic.test}}}
