\name{bayes_group}
\alias{capacityGroup.bayes}
\alias{ucip.bayes}
\alias{micGroup.bayes}
\title{Bayesian Hierarchical Analysis for SFT}
\description{
  Bayesian companions to standard SFT non-parametric measures. \code{capacityGroup.bayes} 
  estimates the population effect of the UCIP/Cz Z-test. \code{ucip.bayes} handles single-subject 
  capacity evaluation using exact analytic posteriors. \code{micGroup.bayes} estimates the 
  normalized population Mean Interaction Contrast (MIC) for cognitive architectures.
}
\usage{
capacityGroup.bayes(inData, CR = NULL, stopping.rule = c("OR", "AND", "STST"), ndraws = 10000L, prior_shape = NULL, prior_rate = NULL, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = NULL, burnin = 1000L, thin = 1L, chains = 4L, rope = NULL, method = c("InvGamma", "HalfNormal", "Centered"), prior_tau_sd = NULL, adapt_delta = .95, max_treedepth = 12L, stan_control = list(), tau2_prior_shape = NULL, tau2_prior_rate = NULL, Condition = NULL, Subject = NULL, score_method = c("standardized", "multiplicative"), var_method = c("analytic", "bootstrap"), n_boot = 2000L)
ucip.bayes(inData, CR = NULL, stopping.rule = c("OR", "AND", "STST"), ndraws = 10000L, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = NULL, chains = 4L, rope = NULL, Condition = NULL, Subject = NULL, score_method = c("standardized", "multiplicative"), var_method = c("analytic", "bootstrap"), n_boot = 2000L)
micGroup.bayes(inData, Condition = NULL, Subject = NULL, correct = TRUE, var_method = c("analytic", "bootstrap"), n_boot = 2000L, ndraws = 10000L, prior_shape = 2, prior_rate = 0.01, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = 0.10, burnin = 1000L, thin = 1L, chains = 4L, rope = NULL, method = c("InvGamma", "HalfNormal", "Centered"), prior_tau_sd = 0.10, adapt_delta = .95, max_treedepth = 12L, stan_control = list(), tau2_prior_shape = NULL, tau2_prior_rate = NULL)
}
\arguments{
  \item{inData}{Standard SFT data frame. Within data frames, \code{subjects},
    \code{rt}, and \code{LogicalRule} are accepted as aliases for
    \code{Subject}, \code{RT}, and \code{Condition}.}
  \item{CR}{Optional correctness indicators.}
  \item{stopping.rule}{Logical stopping rule: \code{"OR"}, \code{"AND"}, \code{"STST"}.}
  \item{ndraws, burnin, thin, chains, seed}{MCMC sampling size and structure controls.}
  \item{prior_shape, prior_rate, prior_mean, prior_sd, prior_tau_sd, tau2_prior_shape, tau2_prior_rate}{Hyperparameters for the normal and inverse-gamma/half-normal hierarchical priors. For \code{capacityGroup.bayes}/\code{ucip.bayes}, \code{prior_sd}, \code{prior_tau_sd}, \code{prior_shape}, and \code{prior_rate} default to \code{NULL} and are resolved to scale-aware values that depend on \code{score_method}: the standardized reference-information (phi) scale is the anchor (\eqn{\mu \sim N(0, 0.5^2)}, \eqn{\tau \sim} HalfNormal\eqn{(0, 1)}, \eqn{\tau^2 \sim} InvGamma\eqn{(2, 0.5)}), and the \code{"multiplicative"} (\eqn{\eta = \log(A/B)}) scale uses \eqn{0.35} for both SDs and \eqn{0.35^2} for the InvGamma rate. Supplying any of these explicitly overrides its default. \code{prior_mean} defaults to \eqn{0} (the UCIP benchmark).}
  \item{hdi}{Highest density interval width.}
  \item{rope}{Region of practical equivalence.}
  \item{method}{Hierarchical sampler: \code{"InvGamma"} (conjugate Gibbs), \code{"HalfNormal"} (non-centred Stan), or \code{"Centered"} (centred Stan).}
  \item{Condition, Subject}{Subsetting controls for the data frame.}
  \item{score_method}{Effect scale: \code{"standardized"} (the default reference-information phi scale) or \code{"multiplicative"} (\eqn{\eta = \log(A/B)}).}
  \item{var_method, n_boot}{Method for estimating participant variance (\code{"analytic"} or \code{"bootstrap"}).}
  \item{correct}{Logical; if TRUE, only correct trials are included in the MIC.}
  \item{adapt_delta, max_treedepth, stan_control}{NUTS adaptation controls.}
}
\value{
  A list including a summary of posterior statistics, the fitted score table (with standard errors, 
  z-scores, and credible intervals), and Bayes Factors for point null and interval null hypotheses.
}
\details{
  These Bayesian equivalents provide partial pooling and rigorous hypothesis tests for capacity and 
  architecture via the Savage-Dickey density ratio. They avoid the problems of strict threshold 
  reliance common in standard frequentist K-S or Z-tests.
}
