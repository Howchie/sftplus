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
resilience.bayes(inData, CR = NULL, ndraws = 10000L,
                 seed = NULL, hdi = .94, prior_mean = 0, prior_sd = NULL, chains = 4L,
                 rope = NULL, Condition = NULL, Subject = NULL,
                 score_method = c("score", "standardized", "multiplicative"),
                 var_method = c("analytic", "bootstrap"), n_boot = 2000L)

resilienceGroup.bayes(inData, CR = NULL, ndraws = 10000L,
                      prior_shape = NULL, prior_rate = NULL, seed = NULL, hdi = .94,
                      prior_mean = 0, prior_sd = NULL, burnin = 1000L, thin = 1L,
                      chains = 4L, rope = NULL,
                      method = c("InvGamma", "HalfNormal", "Centered"),
                      prior_tau_sd = NULL, adapt_delta = .95, max_treedepth = 12L,
                      stan_control = list(), Condition = NULL, Subject = NULL,
                      score_method = c("score", "standardized", "multiplicative"),
                      var_method = c("analytic", "bootstrap"), n_boot = 2000L)
}
\arguments{
  \item{inData}{Canonical SFT data frame, or a list of AB, AY, and XB response-time vectors
    -- the double-target display and each single target accompanied by a distractor --
    nested by participant for the group function.}
  \item{CR}{Optional correctness indicators; omit when \code{inData} is a data frame.}
  \item{score_method}{\code{"score"} (default), the log-hazard-ratio effect
    \eqn{\theta = U/V} with precision \eqn{V}; \code{"standardized"}, the same effect on the
    reference-Cz scale \eqn{\phi = \theta\sqrt{V_{ref}}}; or \code{"multiplicative"},
    \eqn{\eta = \log(A/B)}. All three are exactly zero under the UCIP-OR benchmark.}
  \item{var_method, n_boot}{\code{"analytic"} martingale variance or a within-condition
    \code{"bootstrap"} variance.}
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
  Resilience (Little, Eidels, Fific & Wang, 2015) shares the functional form of the OR
  capacity coefficient,
  \deqn{R(t) = H_{AB}(t) / (H_{AY}(t) + H_{XB}(t)),}
  differing only in which trials enter it: the denominator uses single-target displays that
  still carry distracting information rather than lone targets. Houpt and Little (2016)
  accordingly build its statistical test on the Houpt-Townsend (2012) weighted-logrank score
  statistic, and these functions reuse that machinery verbatim -- so
  \code{resilience.bayes} is \code{\link{ucip.bayes}} applied to the resilience conditions,
  and \code{theta_hat * sqrt(V)} reproduces \code{\link{resilience.test}} exactly.

  There is consequently \strong{no evaluation horizon}. The score statistic accumulates
  risk-set-weighted increments across the whole observed time course. That is what makes it
  well behaved where a cumulative hazard read off at a single late time is not: for complete
  data the Nelson-Aalen cumulative hazard at a sample's own maximum is exactly the harmonic
  number \eqn{\sum_{k=1}^{n} 1/k} \emph{whatever the distribution}, so all three series
  coincide there and the contrast carries no information about the UCIP-OR benchmark.
}
\value{
  A list mirroring \code{\link{ucip.bayes}} or \code{\link{capacityGroup.bayes}}, with the
  per-participant scores in \code{score}, the reference information in \code{V_ref}, and,
  for multi-condition fits, \code{condition_contrasts}.
}
\references{
  Little, D. R., Eidels, A., Fific, M., & Wang, T. S. L. (2015). Understanding the influence
  of distractors on workload capacity. \emph{Journal of Mathematical Psychology},
  \bold{68--69}, 25--36.

  Houpt, J. W., & Little, D. R. (2016). Statistical analysis of the resilience function.
  \emph{Behavior Research Methods}.

  Houpt, J. W., & Townsend, J. T. (2012). Statistical measures for workload capacity
  analysis. \emph{Journal of Mathematical Psychology}, \bold{56}, 341--355.
}
\seealso{\code{\link{resilience.test}}, \code{\link{ucip.bayes}},
  \code{\link{capacityGroup.bayes}}, \code{\link{spike_slab}},
  \code{\link{prior_sensitivity}}}
