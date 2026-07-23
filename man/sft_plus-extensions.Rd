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
\alias{smooth_one_cdf}
\alias{ucip.bayes.test}
\alias{ucip.bayes.single_subject.test}
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
smooth_one_cdf(t, y, smooth_cdf = "none", smooth_spar = .65)
ucip.bayes.test(RT, CR = NULL, stopping.rule = c("OR", "AND", "STST"), ndraws = 10000L, prior_shape = 2, prior_rate = 1, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = 1, burnin = 1000L, thin = 1L, chains = 4L, rope = NULL, method = c("InvGamma", "HalfNormal", "Centered"), prior_tau_sd = 1, adapt_delta = .95, max_treedepth = 12L, stan_control = list(), tau2_prior_shape = NULL, tau2_prior_rate = NULL, Condition = NULL, Subject = NULL)
ucip.bayes.single_subject.test(RT, CR = NULL, stopping.rule = c("OR", "AND", "STST"), ndraws = 10000L, seed = NULL, hdi = .94, prior_mean = 0, prior_sd = 1, chains = 4L, rope = NULL, Condition = NULL, Subject = NULL)
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
    \code{ucip.bayes.test} route can retain multiple subjects.}
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
  \item{method}{For \code{ucip.bayes.test}, the hierarchical route: \code{InvGamma} is the conjugate Gibbs sampler, \code{HalfNormal} is non-centred Stan, and \code{Centered} is centred Stan. For \code{sictestBayes}, \code{DP} is implemented for Bayesian SIC testing, while the legacy \code{IG} route requires external Stan model files.}
  \item{nbin}{Number of bins for the Bayesian SIC/dominance calculations.}
  \item{model}{Simulation model or optional legacy external Bayesian model object.}
  \item{design, salience, logical_rules, sics, cdfs, a_t, subject}{Simulation design, rule, output, and subject controls.}
  \item{stopping.rule}{The logical stopping rule for the UCIP Bayesian companion.}
  \item{ndraws, nsamp, maxn}{Monte Carlo sizes. Reduce these for exploratory work and increase them for final analyses.}
  \item{tolSIC, tolMIC, ci}{Model-classification tolerances and Monte Carlo interval mass for \code{sicDPtest}.}
  \item{prior_shape, prior_rate}{Shape and rate of the inverse-Gamma prior on between-subject variance \eqn{tau^2} in \code{ucip.bayes.test}.}
  \item{tau2_prior_shape, tau2_prior_rate}{Clearer aliases for \code{prior_shape} and \code{prior_rate}; when supplied, these override the older argument names.}
  \item{prior_mean, prior_sd}{Normal prior mean and standard deviation for the population UCIP effect \eqn{mu}.}
  \item{burnin, thin, chains}{Gibbs sampler burn-in, thinning interval, and number of chains for \code{ucip.bayes.test}.}
  \item{prior_tau_sd}{Positive scale of the half-normal prior on \eqn{tau} for the Stan methods.}
  \item{adapt_delta, max_treedepth, stan_control}{NUTS adaptation controls passed to the Stan methods.}
  \item{rope}{Optional positive practical-equivalence threshold on the \eqn{theta} scale. If NULL, no ROPE probability is reported.}
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

  \code{ucip.bayes.test} is a Bayesian companion to the UCIP/Cz z test. It
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

  For one participant, \code{ucip.bayes.test} delegates to
  \code{ucip.bayes.single_subject.test}, which uses the analytic
  Normal-Normal posterior and returns no population hierarchy. The function
  includes prior- and posterior-predictive score draws, split-chain
  \code{Rhat}, bulk and tail effective sample sizes, and posterior shrinkage
  weights. A practical-equivalence
  probability is reported only when a scientifically chosen \code{rope} is
  supplied. A formal Bayes factor for the UCIP null would require an
  explicitly constrained generative model and is not returned here.
}
\keyword{ sft }
