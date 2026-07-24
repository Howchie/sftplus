\name{sft_plus-extensions}
\alias{assessment_ta_or}
\alias{assessment_ta_and}
\alias{assessment_donkin_discrimination}
\alias{build_at_df}
\alias{build_cdf_df}
\alias{build_ct_df}
\alias{build_sic_df}
\alias{capacity.altieri}
\alias{correctfun}
\alias{find_equal_vc_yes}
\alias{get_time_range}
\alias{LogicalRules_OU_rfun}
\alias{LogicalRules_OU_rfun_R}
\alias{LogicalRules_OU_rfun_Rcpp}
\alias{LogicalRules_rfun_lba}
\alias{plot_altieri}
\alias{resilience}
\alias{resilience.test}
\alias{sicGroupBF}
\alias{sicDPtest}
\alias{sictestBayes}
\alias{simulate_sft}
\alias{simulate_sft_group}
\alias{smooth_one_cdf}
\alias{capacityGroup.bayes}
\alias{ucip.bayes}
\alias{micGroup.bayes}
\title{CapacityAND extensions for Systems Factorial Technology}
\description{
  Additional estimators, tests, plotting helpers, and simulation tools bundled
  by \code{sft.plus}. The Bayesian functions return posterior probabilities or
  posterior/prior event-probability ratios; they are not K-S p-values.
}
\usage{
assessment_ta_or(RT, CR = NULL, times = NULL)
assessment_ta_and(RT, CR = NULL, times = NULL)
assessment_donkin_discrimination(RT, CR = NULL, times = NULL)
build_at_df(times, series_specs, probs = c(.05, .95))
build_cdf_df(times, series_specs, probs = c(.05, .95), smooth_cdf = c("none", "mono", "spline"), smooth_spar = .65)
build_ct_df(times, series_specs, probs = c(.05, .95))
build_sic_df(times, series_specs, trim_cdf_tails = TRUE, cdf_tail_cut = 5e-4, smooth_cdf = c("none", "mono", "spline"), smooth_spar = .65)
capacity.altieri(RT, CR = NULL, ratio = TRUE, Condition = NULL, Subject = NULL)
correctfun(d)
find_equal_vc_yes(vno, t, A, b, t0, sv, interval = c(1, 4.5))
get_time_range(rt_list, probs = c(.05, .95), by = .001)
LogicalRules_rfun_lba(n, pars, LogicalRule, shared_capacity = NULL, kappa = NULL, tau = NULL)
LogicalRules_OU_rfun(n, pars, LogicalRule, cross_talk = c(pos = 0, neg = 0), dt = .001, max_time = 3, bridge = TRUE, bridge_exp_cutoff = -20, backend = c("auto", "Rcpp", "R"), adaptive = FALSE, adaptive_mode = c("distance", "ptol"), adapt_factor = 32, adapt_eps = 1.5, p_tol = 1e-9, curv_tol = .1)
LogicalRules_OU_rfun_R(n, pars, LogicalRule, cross_talk = c(pos = 0, neg = 0), dt = .001, max_time = 3, bridge = TRUE, bridge_exp_cutoff = -20)
LogicalRules_OU_rfun_Rcpp(n, pars, LogicalRule, cross_talk = c(pos = 0, neg = 0), dt = .001, max_time = 3, bridge = TRUE, bridge_exp_cutoff = -20, adaptive = FALSE, adaptive_mode = c("distance", "ptol"), adapt_factor = 32, adapt_eps = 1.5, p_tol = 1e-9, curv_tol = .1)
plot_altieri(n, p_vec, smooth_cdf = c("none", "mono", "spline"), smooth_spar = .65)
resilience(RT, CR = NULL, ratio = TRUE, Condition = NULL, Subject = NULL)
resilience.test(RT, CR = NULL, Condition = NULL, Subject = NULL)
sicGroupBF(inData, domtest = "ks", plotSIC = TRUE, ...)
sicDPtest(dat, nbin = NULL, nsamp = 10000L, maxn = 500000L, tolSIC = 5e-2, tolMIC = NULL, ci = .95, seed = NULL)
sictestBayes(HH, HL, LH, LL, method = c("DP", "IG"), model = NULL, nbin = NULL, nsamp = 10000L, maxn = 500000L, seed = NULL, ...)
simulate_sft(model = c("lba", "ou"), n, p_vec, design = NULL, salience = FALSE, logical_rules = NULL, sics = FALSE, cdfs = FALSE, a_t = FALSE, subject = 1, backend = c("auto", "Rcpp", "R"), dt = .001, max_time = 2, bridge = TRUE, bridge_exp_cutoff = -20, adaptive = TRUE, adaptive_mode = c("ptol", "distance"), adapt_factor = 32, adapt_eps = 1.5, p_tol = 1e-8, curv_tol = .05, kappa = NULL, tau = NULL)
simulate_sft_group(model = c("lba", "ou"), n, params, nSubjects = NULL, subject_ids = NULL, combine = TRUE, ...)
smooth_one_cdf(t, y, smooth_cdf = "none", smooth_spar = .65)
capacityGroup.bayes(inData, CR = NULL, stopping.rule = c("OR", "AND", "STST"), ndraws = 10000L, prior_shape = 2, prior_rate = 0.5, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = 1, burnin = 1000L, thin = 1L, chains = 4L, rope = NULL, method = c("InvGamma", "HalfNormal", "Centered"), prior_tau_sd = 1, adapt_delta = .95, max_treedepth = 12L, stan_control = list(), tau2_prior_shape = NULL, tau2_prior_rate = NULL, Condition = NULL, Subject = NULL, score_method = c("score", "multiplicative", "standardized"), var_method = c("analytic", "bootstrap"), n_boot = 2000L)
ucip.bayes(inData, CR = NULL, stopping.rule = c("OR", "AND", "STST"), ndraws = 10000L, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = 1, chains = 4L, rope = NULL, Condition = NULL, Subject = NULL, score_method = c("score", "multiplicative", "standardized"), var_method = c("analytic", "bootstrap"), n_boot = 2000L)
micGroup.bayes(inData, Condition = NULL, Subject = NULL, correct = TRUE, var_method = c("analytic", "bootstrap"), n_boot = 2000L, ndraws = 10000L, prior_shape = 2, prior_rate = 0.01, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = 0.10, burnin = 1000L, thin = 1L, chains = 4L, rope = NULL, method = c("InvGamma", "HalfNormal", "Centered"), prior_tau_sd = 0.10, adapt_delta = .95, max_treedepth = 12L, stan_control = list(), tau2_prior_shape = NULL, tau2_prior_rate = NULL)
}
\arguments{
  \item{RT, CR}{Response times and optional correctness indicators in the standard SFT list format.}
  \item{times}{Evaluation times for curves.}
  \item{t, y}{Time and CDF vectors for smoothing.}
  \item{series_specs}{A list of named curve/data specifications.}
  \item{probs}{Quantile probabilities used to trim plotted curves.}
  \item{smooth_cdf, smooth_spar}{CDF smoothing method and spline smoothness.}
  \item{trim_cdf_tails, cdf_tail_cut}{Whether and how aggressively to trim common CDF tails.}
  \item{ratio}{Whether to return ratio capacity or additive capacity.}
  \item{Condition, Subject}{When \code{RT} is a row-wise SFT data frame,
    optionally select one condition and one subject for conversion to list
    input.  \code{Condition} must select one condition; the hierarchical
    \code{capacityGroup.bayes} route can retain multiple subjects.}
  \item{d}{A simulated data frame for correctness coding.}
  \item{vno, A, b, t0, sv, interval}{Parameters and search interval for the LBA matching helper.}
  \item{rt_list, by}{RT vectors and grid spacing for time-range construction.}
  \item{n}{Number of simulated trials.}
  \item{pars}{Accumulator parameter matrix.}
  \item{LogicalRule}{Logical stopping rule.}
  \item{shared_capacity}{Optional list or named vector with \code{kappa}, \code{tau}, and \code{active}. The lower-level LBA function applies it to the two target racers only when explicitly activated.}
  \item{kappa, tau}{Optional LBA shared-capacity parameters. \code{kappa} is the mean multiplier and \code{tau} is the standard deviation of one shared trial-level normal multiplier. In \code{simulate_sft()}, explicit values override same-named entries in \code{p_vec}.}
  \item{cross_talk, dt, max_time, bridge, bridge_exp_cutoff}{OU interaction and numerical integration controls.}
  \item{backend, adaptive, adaptive_mode, adapt_factor, adapt_eps, p_tol, curv_tol}{Backend and adaptive integration controls.}
  \item{p_vec}{Named simulation parameter vector.}
  \item{inData}{Standard SFT data frame.}
  \item{domtest, plotSIC}{Dominance test method and plotting switch.}
  \item{HH, HL, LH, LL}{The four SIC response-time samples.}
  \item{dat}{A four-element list of SIC samples in HH, HL, LH, LL order.}
  \item{method}{For \code{capacityGroup.bayes}, the hierarchical route: \code{InvGamma} is the conjugate Gibbs sampler, \code{HalfNormal} is non-centred Stan, and \code{Centered} is centred Stan. For \code{sictestBayes}, \code{DP} is implemented for Bayesian SIC testing, while the legacy \code{IG} route requires external Stan model files.}
  \item{nbin}{Number of bins for the Bayesian SIC/dominance calculations.}
  \item{model}{Simulation model or optional legacy external Bayesian model object.}
  \item{design, salience, logical_rules, sics, cdfs, a_t, subject}{Simulation design, rule, output, and subject controls.}
  \item{params}{Per-subject simulation parameters for \code{simulate_sft_group}. Either a named numeric vector (one \code{p_vec}, recycled to every subject) or a matrix/data frame with one row per subject and named columns giving each subject's process parameters (drifts, thresholds, non-decision time, etc.).}
  \item{nSubjects}{Number of subjects for \code{simulate_sft_group}. Required (and inferred as 1 if omitted) when \code{params} is a vector; when \code{params} is a matrix it defaults to the number of rows and, if supplied, must match it.}
  \item{subject_ids}{Optional character vector of subject labels for \code{simulate_sft_group}. Defaults to the row names of \code{params} when present, otherwise zero-padded \code{"Subject001"}, \code{"Subject002"}, ... labels.}
  \item{combine}{For \code{simulate_sft_group}, whether to return one row-bound data frame (the default, directly usable by the group-level functions) or the list of per-subject \code{simulate_sft} results.}
  \item{stopping.rule}{The logical stopping rule for the UCIP Bayesian companion.}
  \item{score_method}{Per-participant UCIP effect scale. \code{"score"} (default)
    uses the raw score effect, \code{"multiplicative"} uses the log weighted
    capacity ratio, and \code{"standardized"} uses the reference-information
    \eqn{phi} score effect. The historical \code{"capacity"} aliases remain
    accepted for \code{"multiplicative"}.}
  \item{ndraws, nsamp, maxn}{Monte Carlo sizes. Reduce these for exploratory work and increase them for final analyses.}
  \item{tolSIC, tolMIC, ci}{Model-classification tolerances and Monte Carlo interval mass for \code{sicDPtest}.}
  \item{prior_shape, prior_rate}{Shape and rate of the inverse-Gamma prior on between-subject variance \eqn{tau^2} in \code{capacityGroup.bayes}.}
  \item{tau2_prior_shape, tau2_prior_rate}{Clearer aliases for \code{prior_shape} and \code{prior_rate}; when supplied, these override the older argument names.}
  \item{prior_mean, prior_sd}{Normal prior mean and standard deviation for the population UCIP effect \eqn{mu}. For \code{micGroup.bayes} the effect is the dimensionless relative MIC, so the default \code{prior_sd = 0.10} encodes population MICs within roughly \eqn{\pm 20\%} of mean RT rather than a millisecond scale.}
  \item{var_method, n_boot}{Sampling-variance method for the per-participant effect precision in \code{ucip.bayes}, \code{capacityGroup.bayes}, and \code{micGroup.bayes}. \code{"analytic"} (default) uses the closed-form variance (martingale/delta method for the UCIP score and log-capacity, delta method for the relative MIC); \code{"bootstrap"} replaces it with a within-participant resampling variance over \code{n_boot} replicates, a distribution-robust alternative reported alongside the analytic value.}
  \item{correct}{For \code{micGroup.bayes}, keep only correct trials so the estimand is the conditional-correct interaction contrast.}
  \item{burnin, thin, chains}{Gibbs sampler burn-in, thinning interval, and number of chains for \code{capacityGroup.bayes}.}
  \item{prior_tau_sd}{Positive scale of the half-normal prior on \eqn{tau} for the Stan methods.}
  \item{adapt_delta, max_treedepth, stan_control}{NUTS adaptation controls passed to the Stan methods.}
  \item{rope}{Optional positive practical-equivalence threshold on the active effect scale. For \code{score_method = "standardized"}, it is on the reference-information \eqn{phi} scale; if NULL, no ROPE probability is reported.}
  \item{hdi}{Credible interval mass for the UCIP posterior summary. Intervals are highest-density intervals.}
  \item{seed}{Optional reproducibility seed; the previous random-number state is restored on return.}
  \item{...}{Additional plotting or Bayesian-control arguments.}
}
\value{
  The estimator functions return functions and metadata. Bayesian functions
  return posterior summaries and draws, with no frequentist p-value implied.
}
\details{
  The LBA simulator supports an optional shared-capacity effect on AB cells.
  For each AB trial it draws one \eqn{Z \sim N(0,1)} and uses
  \eqn{V_i = (\kappa + \tau Z)v_i + \epsilon_i} for the active target
  accumulators, where the independent \eqn{\epsilon_i} terms use their usual
  LBA drift standard deviations. Thus the target drift rates are bivariate
  normal before the LBA positive-drift handling, with positive trial-level
  co-fluctuation induced by \eqn{\tau}. Nontarget accumulators and AN, NB, and
  NN cells are unaffected. The defaults are \eqn{\kappa = 1} and
  \eqn{\tau = 0}. The parameter names \code{capacity_kappa},
  \code{capacity_target_kappa}, \code{shared_capacity_kappa},
  \code{shared_kappa}, \code{capacity_tau}, \code{capacity_target_tau},
  \code{shared_capacity_tau}, and \code{shared_tau} are also accepted in
  \code{p_vec}.

  \code{simulate_sft_group} generates data for several participants by calling
  \code{simulate_sft} once per subject and row-binding the results. Supply the
  process parameters as a subject-by-parameter matrix (one row per participant,
  named columns), or as a single named \code{params} vector together with
  \code{nSubjects} to simulate that many participants from a common parameter
  set. Each participant is labelled through the \code{Subject} column and, by
  default, the function returns a single data frame with \code{RT},
  \code{Channel1}, \code{Channel2}, \code{Correct}, \code{Subject}, and
  \code{Condition} columns (plus the per-trial bookkeeping columns produced by
  \code{simulate_sft}), ready to pass to \code{capacityGroup},
  \code{assessmentGroup}, or \code{sicGroup}. Passing \code{combine = FALSE}
  returns the list of full per-subject \code{simulate_sft} results instead.

  \code{capacityGroup.bayes} is a Bayesian companion to the UCIP/Cz z test. It
  extracts the UCIP score numerator \eqn{U_i} and information \eqn{V_i} for
  each participant, then fits the score/information hierarchy
  \eqn{U_i | theta_i ~ N(V_i theta_i, V_i)},
  \eqn{theta_i | mu,tau^2 ~ N(mu,tau^2)}. The population effect \eqn{mu} is
  given a Normal prior. With \code{method = "InvGamma"}, \eqn{tau^2} has an
  inverse-Gamma prior and the model is fit by conjugate Gibbs sampling. With
  \code{method = "HalfNormal"} or \code{method = "Centered"}, \eqn{tau} has
  a half-normal prior and the model is fit in Stan using, respectively, a
  non-centred or centred parameterization. All three methods use partial
  pooling on the score effect,
  not a model on pooled response-time distributions and not a terminal
  cumulative-hazard contrast relabeled as \code{Cz}.

  For one participant, \code{RT} has the ordinary package format: a list of
  condition-wise response-time vectors. For group inference, supply a nested
  list, with one such condition list per participant; \code{CR} must have the
  same nesting when supplied. The returned \code{score} table contains the
  original \eqn{U_i}, \eqn{V_i}, \eqn{U_i/V_i}, standard error, and frequentist
  \code{Cz} score. Posterior probabilities and HDIs are reported separately
  for the population effect and each participant effect.

  For one participant, \code{capacityGroup.bayes} delegates to
  \code{ucip.bayes}, which uses the analytic
  Normal-Normal posterior and returns no population hierarchy. The function
  includes prior- and posterior-predictive score draws, split-chain
  \code{Rhat}, bulk and tail effective sample sizes, and posterior shrinkage
  weights. A practical-equivalence
  probability is reported only when a scientifically chosen \code{rope} is
  supplied.

  The \code{bayes_factor} element reports Bayes factors that, unlike a
  directional posterior probability, quantify evidence \emph{for} as well as
  against the null. Which test appears is matched to the effect scale.
  \code{point_null} is the Savage-Dickey density ratio for the exact UCIP
  prediction (population \eqn{mu = 0}, or the single participant effect \eqn{= 0}
  with \eqn{tau} left free in the hierarchy); for the hierarchy it is
  Rao-Blackwellised from the exact Gibbs full conditional of \eqn{mu}, so it is
  valid for all three samplers. \code{interval_null} is the interval Bayes
  factor for the practical-equivalence region \eqn{|effect| < \code{rope}},
  crediting the prior-to-posterior odds shift into that region; it is reported
  only when a \code{rope} is supplied. The point null suits the raw \code{score}
  effect and the dimensionless \code{multiplicative} (\eqn{eta}) effect; the
  interval null suits \code{multiplicative} and \code{standardized}, but not
  the information-referenced raw \code{score}. \code{BF01} favours the null and
  \code{BF10} the alternative.
  Because the point-null Bayes factor depends on \code{prior_sd} (the
  Jeffreys-Lindley sensitivity), a diffuse alternative prior inflates the
  evidence for the null; interpret it together with the posterior summary.

  \code{score_method} selects the per-participant effect that is pooled.
  \code{"score"} (default) models the raw score effect \eqn{U_i / V_i};
  \code{"multiplicative"} models \eqn{eta_i = \log(A_i / B_i)}, the log weighted
  capacity ratio, with its delta-method precision. \code{"multiplicative"} was
  formerly named \code{"capacity"}, which (with \code{"logcapacity"} and
  \code{"eta"}) is still accepted as an alias. \code{"standardized"} models
  the reference-information score effect \eqn{phi_i = (U_i/V_i)\sqrt{V_ref}},
  where \eqn{V_ref = median(V_i)} is computed over the participants retained
  for the fit. Its observation precision is \eqn{V_i/V_ref}, so the frequentist
  \code{Cz} statistic is unchanged and the effect is invariant to a common
  rescaling of score weights. One \eqn{phi}-unit is a one-\code{Cz}-unit shift
  for a participant with median information, making the alternative prior and
  interval null coherent. This remains an information-referenced statistical
  scale, not a physical response-time scale. With a \code{rope},
  \code{"standardized"} reports both \code{point_null} and
  \code{interval_null} Bayes factors.

  \code{micGroup.bayes} is the mean-interaction-contrast companion. It reads the
  salience \eqn{2\times2} cells from \code{Channel1} and \code{Channel2}
  (\code{2 = High}, \code{1 = Low}), forms each participant's cell means, and
  pools the \emph{relative} MIC
  \eqn{rho_i = MIC_i / M_i} with the same normal-normal hierarchy, where
  \eqn{MIC_i = bar T_{HH} - bar T_{HL} - bar T_{LH} + bar T_{LL}} and
  \eqn{M_i} is the equal-weight grand-mean RT of the four cells. Normalising by
  \eqn{M_i} makes the estimand dimensionless and invariant to multiplicative
  differences in participant speed, while preserving the sign, the exact-null
  boundary (\eqn{rho = 0} iff \eqn{MIC = 0}), and the interpretation of the MIC
  as the (now normalised) signed area under the survivor interaction contrast.
  A positive \eqn{rho} is over-additive (parallel self-terminating or coactive),
  a negative \eqn{rho} is under-additive (parallel exhaustive), and \eqn{rho}
  near zero is additive (serial). The \code{rope} is a proportional
  practical-equivalence threshold: \code{rope = 0.025} treats interactions below
  \eqn{2.5\%} of mean RT as negligible, and it is also the region for the
  \code{interval_null} Bayes factor. Because \eqn{rho = 0} is a knife-edge
  boundary between processing architectures rather than a scientific point
  prediction, \code{micGroup.bayes} reports only the interval Bayes factor and
  no point null. A single participant is fitted with the
  analytic posterior, and the sampler, diagnostics, predictive checks, and
  shrinkage machinery are shared with \code{capacityGroup.bayes}.
}
\keyword{ sft }
