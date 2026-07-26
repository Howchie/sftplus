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
capacityGroup.bayes(inData, CR = NULL, stopping.rule = c("OR", "AND", "STST"), ndraws = 10000L, prior_shape = NULL, prior_rate = NULL, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = NULL, burnin = 1000L, thin = 1L, chains = 4L, rope = NULL, method = c("InvGamma", "HalfNormal", "Centered"), prior_tau_sd = NULL, adapt_delta = .95, max_treedepth = 12L, stan_control = list(), tau2_prior_shape = NULL, tau2_prior_rate = NULL, Condition = NULL, Subject = NULL, score_method = c("score", "standardized", "multiplicative"), var_method = c("analytic", "bootstrap"), n_boot = 2000L)
ucip.bayes(inData, CR = NULL, stopping.rule = c("OR", "AND", "STST"), ndraws = 10000L, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = NULL, chains = 4L, rope = NULL, Condition = NULL, Subject = NULL, score_method = c("score", "standardized", "multiplicative"), var_method = c("analytic", "bootstrap"), n_boot = 2000L)
micGroup.bayes(inData, Condition = NULL, Subject = NULL, correct = TRUE, var_method = c("analytic", "bootstrap"), n_boot = 2000L, ndraws = 10000L, prior_shape = 2, prior_rate = 0.01, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = 0.10, burnin = 1000L, thin = 1L, chains = 4L, rope = NULL, method = c("InvGamma", "HalfNormal", "Centered"), prior_tau_sd = 0.10, adapt_delta = .95, max_treedepth = 12L, stan_control = list(), tau2_prior_shape = NULL, tau2_prior_rate = NULL)
}
\arguments{
  \item{inData}{Standard SFT data frame. Within data frames, \code{subjects},
    \code{rt}, and \code{LogicalRule} are accepted as aliases for
    \code{Subject}, \code{RT}, and \code{Condition}.}
  \item{CR}{Optional correctness indicators.}
  \item{stopping.rule}{Logical stopping rule: \code{"OR"}, \code{"AND"}, \code{"STST"}.}
  \item{ndraws, burnin, thin, chains, seed}{MCMC sampling size and structure controls.}
  \item{prior_shape, prior_rate, prior_mean, prior_sd, prior_tau_sd, tau2_prior_shape, tau2_prior_rate}{Hyperparameters for the normal and inverse-gamma/half-normal hierarchical priors. For \code{capacityGroup.bayes}/\code{ucip.bayes} these default to \code{NULL} and are resolved once the data are scored. The anchor is the \code{"score"} (\eqn{\theta}) log-hazard-ratio scale: \eqn{\mu \sim N(0, 0.5^2)}, \eqn{\tau \sim} HalfNormal\eqn{(0, 0.5)}, \eqn{\tau^2 \sim} InvGamma\eqn{(2, 0.25)}. A prior SD of \eqn{0.5} on a log hazard ratio spans roughly \eqn{[0.37, 2.7]} in hazard-ratio terms at 95\%, which is weakly informative on the conventional scale. For \code{"standardized"} every width is pushed through the \eqn{\sqrt{V_{ref}}} map (variances through \eqn{V_{ref}}), so the two scales describe one model. The \code{"multiplicative"} (\eqn{\eta = \log(A/B)}) scale uses \eqn{0.35} for both SDs and \eqn{0.35^2} for the InvGamma rate. Supplying any of these explicitly overrides its default. \code{prior_mean} defaults to \eqn{0} (the UCIP benchmark).

    Note that a \emph{fixed} width on the \eqn{\phi} scale is not a sample-size-free prior: it implies a prior on \eqn{\theta} of width \eqn{prior\_sd / \sqrt{V_{ref}}}, which tightens at the same rate the likelihood sharpens and so is never outweighed by data. Supply \code{prior_sd} on the \eqn{\phi} scale only if you intend that.}
  \item{hdi}{Highest density interval width.}
  \item{rope}{Region of practical equivalence.}
  \item{method}{Hierarchical sampler: \code{"InvGamma"} (conjugate Gibbs), \code{"HalfNormal"} (non-centred Stan), or \code{"Centered"} (centred Stan).}
  \item{Condition, Subject}{Subsetting controls for the data frame. \code{Condition = NULL} now fits \emph{every} condition present, giving one population mean \eqn{\mu_c} per condition over a shared between-participant scale \eqn{\tau}, plus every pairwise contrast in \code{condition_contrasts}; supplying a single value reproduces the one-condition fit exactly. Participants are treated as independent given \eqn{(\mu, \tau)}, so a participant appearing in several conditions is not linked across them: the contrasts are between-condition population comparisons, not within-subject differences.}
  \item{score_method}{Effect scale for the pooled UCIP effect.

    \code{"score"} (default) pools \eqn{\theta_i = U_i / V_i}, the Peto one-step estimate of the log hazard ratio for the UCIP contrast, with precision \eqn{V_i}. This is the effect the \eqn{Cz} statistic tests: \eqn{Cz_i = \theta_i \sqrt{V_i}} reproduces \code{\link{ucip.test}} exactly. Its location is stable as trials accumulate while \eqn{V_i} grows linearly with them, which is what makes a fixed prior width a genuine effect-size prior.

    \code{"standardized"} reports the same fit on the reference-\eqn{Cz} scale \eqn{\phi_i = \theta_i \sqrt{V_{ref}}}, where \eqn{V_{ref}} is the median \eqn{V_i} over the fitted participants. Because that is a deterministic linear map given the data, and the default priors are mapped with it, this returns the identical posterior in \eqn{Cz} units. Note that \eqn{\phi} is a z-statistic: it grows as \eqn{\sqrt{n}} with its precision pinned at 1, so a \code{rope} on this scale is design-specific.

    \code{"multiplicative"} instead uses \eqn{\eta = \log(A/B)}, the log weighted-capacity ratio, with its delta-method precision. This is a different estimand, not a reparameterisation of \eqn{Cz}. \code{"capacity"} remains an accepted alias.}
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
