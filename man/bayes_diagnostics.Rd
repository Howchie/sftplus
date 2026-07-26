\name{bayes_diagnostics}
\alias{prior_sensitivity}
\alias{spike_slab}
\alias{sft_waic}
\title{Prior Sensitivity, Per-Participant Nulls, and Model Comparison}
\description{
  Follow-up analyses for fitted Bayesian SFT models. \code{prior_sensitivity} traces how the
  null Bayes factors depend on the alternative-prior width. \code{spike_slab} asks which
  individual participants depart from the null. \code{sft_waic} returns WAIC and PSIS-LOO
  for a semiparametric hazard fit.
}
\usage{
prior_sensitivity(fit, prior_sd = NULL, ndraws = 4000L, burnin = 500L,
                  chains = 2L, seed = NULL, rope = NULL)

spike_slab(fit, w_shape1 = 1, w_shape2 = 1, ndraws = 10000L, burnin = 1000L,
           thin = 1L, chains = 4L, seed = NULL, hdi = NULL)

sft_waic(object)
}
\arguments{
  \item{fit}{A fitted object from \code{\link{capacityGroup.bayes}},
    \code{\link{micGroup.bayes}}, \code{\link{resilienceGroup.bayes}},
    \code{\link{ucip.bayes}}, or \code{\link{resilience.bayes}}. Both functions refit from
    the stored per-participant estimates and precisions, so neither rescores response times.}
  \item{prior_sd}{Grid of alternative-prior standard deviations; \code{NULL} spans an
    eightfold range either side of the fitted value.}
  \item{rope}{Practical-equivalence threshold for the interval null; defaults to the fit's.}
  \item{w_shape1, w_shape2}{Beta prior on the inclusion rate \eqn{w}. The default
    \code{Beta(1, 1)} gives a prior inclusion probability of one half.}
  \item{ndraws, burnin, thin, chains, seed, hdi}{Sampler and summary controls.}
  \item{object}{A fitted \code{sft_bayes} object from
    \code{\link{semiparametricSFT.bayes}}.}
}
\details{
  \strong{Prior sensitivity.} A point-null Bayes factor is not prior-free: a more diffuse
  alternative spreads its mass and mechanically favours the null (Jeffreys-Lindley), so the
  honest summary is a curve rather than a number. The interval null is prior-sensitive too,
  and not necessarily less so, since widening the alternative also lowers the prior mass
  inside the region. The two depend on the prior differently -- the point null through the
  prior density at zero, the interval null through the prior mass in the ROPE -- so
  agreement across the grid is meaningful in a way either number alone is not.

  \strong{Spike-and-slab.} The hierarchy's point null asks whether the population mean is
  zero with \eqn{\tau} free; it cannot say which participants individually depart from the
  benchmark, and a large \eqn{\tau} can absorb the very signal one wants to detect. Here
  each participant's effect is either exactly zero (the spike) or drawn from the population
  slab, and the posterior inclusion probability \eqn{P(z_i = 1)} is that participant's
  evidence. Because the inclusion rate \eqn{w} is estimated, participants shrink toward one
  another, so an isolated departure in an otherwise null sample is discounted -- the
  multiplicity control that independent per-participant tests lack.

  \strong{WAIC and LOO.} Computed from the pointwise Poisson log-likelihood of the fitted
  hazard model, for comparing fits that differ in structure (pooled versus salience-split,
  bin count, basis dimension, smoothness priors). They are \emph{cell-wise}: they estimate
  predictive accuracy for a new time bin of an observed participant, not for a new
  participant. Compare only fits built on the same pooled grid and observation set. PSIS-LOO
  requires the optional \pkg{loo} package and is preferable where available, because its
  Pareto-\eqn{k} diagnostic warns when the pointwise approximation fails and WAIC has no
  equivalent.
}
\value{
  \code{prior_sensitivity} returns a data frame of class \code{sft_prior_sensitivity} with
  one row per grid point. \code{spike_slab} returns a list with a per-participant
  \code{summary} (inclusion probability and Bayes factor), a \code{population} summary of
  \eqn{w}, \eqn{\mu}, and \eqn{\tau}, the \code{draws}, and \code{diagnostics}.
  \code{sft_waic} returns a list with \code{waic} and optional \code{loo}.
}
\seealso{\code{\link{capacityGroup.bayes}}, \code{\link{semiparametricSFT.bayes}}}
