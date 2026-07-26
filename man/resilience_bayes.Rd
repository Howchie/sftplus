\name{resilience_bayes}
\alias{resilience.bayes}
\alias{resilienceGroup.bayes}
\title{Bayesian Resilience Analysis}
\description{
  Bayesian companions to \code{\link{resilience.test}}. \code{resilience.bayes} gives one
  participant an analytic normal-normal posterior; \code{resilienceGroup.bayes} pools
  participants with the same hierarchy as \code{\link{capacityGroup.bayes}} and, when the
  data carry more than one condition, estimates one population mean per condition together
  with their pairwise contrasts.
}
\usage{
resilience.bayes(inData, CR = NULL, at = NULL, at_quantile = 0.5, ndraws = 10000L,
                 seed = NULL, hdi = .94, prior_mean = 0, prior_sd = NULL, chains = 4L,
                 rope = NULL, Condition = NULL, Subject = NULL,
                 score_method = c("multiplicative", "difference"),
                 var_method = c("analytic", "bootstrap"), n_boot = 2000L)

resilienceGroup.bayes(inData, CR = NULL, at = NULL, at_quantile = 0.5, ndraws = 10000L,
                      prior_shape = NULL, prior_rate = NULL, seed = NULL, hdi = .94,
                      prior_mean = 0, prior_sd = NULL, burnin = 1000L, thin = 1L,
                      chains = 4L, rope = NULL,
                      method = c("InvGamma", "HalfNormal", "Centered"),
                      prior_tau_sd = NULL, adapt_delta = .95, max_treedepth = 12L,
                      stan_control = list(), Condition = NULL, Subject = NULL,
                      score_method = c("multiplicative", "difference"),
                      var_method = c("analytic", "bootstrap"), n_boot = 2000L)
}
\arguments{
  \item{inData}{Canonical SFT data frame, or a list of AB, A, and B response-time vectors
    (nested by participant for the group function).}
  \item{CR}{Optional correctness indicators; omit when \code{inData} is a data frame.}
  \item{at}{Evaluation horizon. \code{NULL} derives it from \code{at_quantile}.}
  \item{at_quantile}{Pooled response-time quantile defining the horizon when \code{at} is
    \code{NULL} (default the median). See \sQuote{Details} for why the terminal time is not
    the default.}
  \item{score_method}{\code{"multiplicative"} (default), the log resilience ratio
    \eqn{\psi = \log(H_{AB} / (H_A + H_B))}; or \code{"difference"}, the cumulative-hazard
    contrast \eqn{\delta = H_{AB} - H_A - H_B}. Both are exactly zero under the UCIP-OR
    benchmark.}
  \item{var_method, n_boot}{\code{"analytic"} Nelson-Aalen variance (delta-method for
    \eqn{\psi}) or a within-condition \code{"bootstrap"} variance.}
  \item{Condition}{Condition selector. \code{NULL} fits every condition present in a data
    frame; supplying one value reproduces the single-condition fit.}
  \item{Subject}{Optional subject selector.}
  \item{ndraws, burnin, thin, chains, seed, hdi}{Sampler and summary controls.}
  \item{prior_mean, prior_sd, prior_shape, prior_rate, prior_tau_sd}{Hierarchical prior
    hyperparameters; \code{NULL} resolves to a scale-aware default.}
  \item{rope}{Optional region of practical equivalence.}
  \item{method}{Hierarchical sampler: \code{"InvGamma"} (conjugate Gibbs),
    \code{"HalfNormal"}, or \code{"Centered"} (both Stan).}
  \item{adapt_delta, max_treedepth, stan_control}{NUTS adaptation controls.}
}
\details{
  The resilience contrast compares the dual-target cumulative hazard against the UCIP-OR
  prediction \eqn{H_{AB} = H_A + H_B}. Both effect scales are zero exactly when that holds,
  so the point null is the UCIP-OR benchmark on either.

  \strong{The evaluation horizon matters, and the terminal time is a bad choice.} For
  complete (uncensored) data the Nelson-Aalen cumulative hazard at a sample's own maximum is
  exactly the harmonic number \eqn{\sum_{k=1}^{n} 1/k} \emph{whatever the distribution}.
  All three series therefore coincide there, the contrast collapses to \eqn{-H_n} for every
  dataset, and \eqn{\psi} is \eqn{-\log 2} identically. Censoring softens this but the tail
  stays uninformative. These functions accordingly default to a central quantile of the
  pooled response times, clipped to stay inside every participant's observed range.
  \code{\link{resilience.test}} evaluates at the terminal time and is affected by this.

  \eqn{\psi} is the default scale for pooling because it is horizon-robust: it compares two
  cumulative hazards to each other rather than reporting an absolute level.  \eqn{\delta} is
  an absolute quantity and is comparable across participants only at a shared horizon, which
  is why the group function always evaluates every participant at one time.
}
\value{
  A list mirroring \code{\link{ucip.bayes}} or \code{\link{capacityGroup.bayes}}, with the
  evaluation horizon in \code{at} and, for multi-condition fits,
  \code{condition_contrasts}.
}
\seealso{\code{\link{resilience.test}}, \code{\link{capacityGroup.bayes}},
  \code{\link{spike_slab}}, \code{\link{prior_sensitivity}}}
