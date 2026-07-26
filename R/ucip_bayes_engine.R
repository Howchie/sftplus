# Hierarchical Bayesian UCIP estimation engine: samplers, diagnostics,
# predictive checks, and shrinkage helpers. No public entry points.


.sft_subject_ucip_input <- function(RT, CR = NULL) {
  if (!is.list(RT) || !length(RT)) stop("RT must be a non-empty list.")
  nested <- all(vapply(RT, is.list, logical(1)))
  if (!nested) {
    return(list(RT = list(RT), CR = list(CR), subject = "1"))
  }
  if (is.null(CR)) CR <- vector("list", length(RT))
  if (!is.list(CR) || length(CR) != length(RT) ||
      any(vapply(CR, function(x) !is.null(x) && !is.list(x), logical(1)))) {
    stop("For nested RT input, CR must be a matching list of participant CR lists.")
  }
  subject <- names(RT)
  if (is.null(subject) || any(!nzchar(subject))) subject <- as.character(seq_along(RT))
  list(RT = RT, CR = CR, subject = subject)
}


# Generic normal-normal hierarchical Gibbs sampler. It is agnostic to the
# effect being pooled (UCIP score/Cz, log-capacity, resilience, or the mean
# interaction contrast): every caller supplies a per-unit effect estimate and its
# sampling precision, so the updates below are unchanged across models.
#
# `group` optionally assigns each unit to a condition, giving one population mean
# mu_c per condition over a shared between-subject scale tau:
#
#   theta_hat_i ~ Normal(theta_i, 1 / precision_i)
#   theta_i     ~ Normal(mu_{group_i}, tau)
#   mu_c        ~ Normal(prior_mean, prior_sd)
#
# Every full conditional stays conjugate; with one condition this reduces
# exactly to the single-mean sampler. Note that the model treats units as
# independent given mu and tau, so a participant contributing to several
# conditions is not linked across them -- the contrasts below are
# between-condition population comparisons, not within-subject differences.
.sft_hierarchical_normal_chain <- function(theta_hat, precision, ndraws, burnin,
                                           thin, prior_mean, prior_sd,
                                           tau2_prior_shape, tau2_prior_rate,
                                           group = NULL) {
  n <- length(theta_hat)
  if (is.null(group)) group <- rep(1L, n)
  group <- as.integer(group)
  n_group <- max(group)
  counts <- tabulate(group, n_group)
  niter <- burnin + ndraws * thin
  mu_draws <- matrix(NA_real_, nrow = ndraws, ncol = n_group)
  tau2_draws <- numeric(ndraws)
  theta_draws <- matrix(NA_real_, nrow = ndraws, ncol = n)
  # Over-dispersed per-chain starting values so split-R-hat reflects genuine
  # mixing rather than a shared deterministic start. Each chain advances the
  # global RNG, so tau2 (drawn from its prior) and the jittered mu differ across
  # chains; theta is resampled at the top of the loop, so its start is immaterial.
  tau2 <- 1 / stats::rgamma(1L, shape = tau2_prior_shape, rate = tau2_prior_rate)
  mu <- stats::rnorm(n_group, stats::weighted.mean(theta_hat, precision), sqrt(tau2))
  save_i <- 0L

  for (iter in seq_len(niter)) {
    theta_var <- 1 / (precision + 1 / tau2)
    theta_mean <- theta_var * (precision * theta_hat + mu[group] / tau2)
    theta <- stats::rnorm(n, theta_mean, sqrt(theta_var))

    group_sum <- as.numeric(rowsum(theta, group, reorder = TRUE))
    mu_var <- 1 / (counts / tau2 + 1 / prior_sd^2)
    mu_mean <- mu_var * (group_sum / tau2 + prior_mean / prior_sd^2)
    mu <- stats::rnorm(n_group, mu_mean, sqrt(mu_var))

    # tau^2 has an inverse-Gamma full conditional. Sampling its reciprocal
    # as Gamma keeps the update stable and avoids a parameter at zero.
    tau2 <- 1 / stats::rgamma(
      1L, shape = tau2_prior_shape + n / 2,
      rate = tau2_prior_rate + sum((theta - mu[group])^2) / 2
    )

    if (iter > burnin && (iter - burnin) %% thin == 0L) {
      save_i <- save_i + 1L
      mu_draws[save_i, ] <- mu
      tau2_draws[save_i] <- tau2
      theta_draws[save_i, ] <- theta

      if (save_i == ndraws) break
    }
  }
  list(mu = mu_draws, tau = sqrt(tau2_draws), theta = theta_draws)
}


# Resolve the var_method argument shared by ucip.bayes() and micGroup.bayes().
# "analytic" keeps the closed-form sampling variance (the martingale variance
# for the UCIP score, the delta-method variance for log-capacity, and the
# sum of cell mean-variances for the MIC). "bootstrap" replaces it with a
# within-subject resampling variance, a distribution-robust alternative that is
# useful when RTs are skewed, cells are small, or trials are clustered.
.sft_var_method <- function(method) {
  key <- tolower(gsub("[^a-z]", "", as.character(method[[1L]])))
  switch(key,
         analytic = "analytic", analytical = "analytic",
         martingale = "analytic", delta = "analytic", closedform = "analytic",
         bootstrap = "bootstrap", boot = "bootstrap", resample = "bootstrap",
         stop("var_method must be one of 'analytic' or 'bootstrap'."))
}


# Stan data name for the per-subject precision, kept distinct per effect scale so
# the generated model code and its cache key stay unambiguous.
.sft_precision_name <- function(score_method) {
  switch(score_method, multiplicative = "P", standardized = "W", "V")
}


# One-line reading of the effect each score_method pools.
.sft_effect_interpretation <- function(score_method) {
  switch(
    score_method,
    multiplicative = "exp(eta) is the participant weighted capacity ratio",
    standardized = paste("phi is the score effect on the reference-Cz scale",
                         "(Cz for a participant with median information);",
                         "phi = theta * sqrt(V_ref)"),
    paste("theta is the UCIP score effect on the log-hazard-ratio scale",
          "(the Peto one-step estimate); Cz_i = theta_i * sqrt(V_i)"))
}


.sft_validate_boot <- function(var_method, n_boot) {
  if (identical(var_method, "bootstrap") &&
      (length(n_boot) != 1L || !is.finite(n_boot) || n_boot < 100L ||
       n_boot != as.integer(n_boot))) {
    stop("n_boot must be an integer >= 100 when var_method = 'bootstrap'.",
         call. = FALSE)
  }
  invisible(TRUE)
}


# Building blocks for the single-subject analytic normal-normal posterior,
# shared by ucip.bayes() and the one-subject branch of micGroup.bayes(). The
# caller supplies the scalar effect estimate and its precision and adds the
# model-specific labels to the returned pieces.
.sft_normal_analytic_fit <- function(effect_hat, precision, subject, effect_name,
                                     ndraws, chains, prior_mean, prior_sd,
                                     hdi, rope) {
  post_var <- 1 / (precision + 1 / prior_sd^2)
  post_mean <- post_var * (precision * effect_hat + prior_mean / prior_sd^2)
  effect_chain <- matrix(stats::rnorm(ndraws * chains, mean = post_mean,
                                      sd = sqrt(post_var)),
                         nrow = ndraws, ncol = chains)
  effect_draws <- matrix(as.vector(effect_chain), ncol = 1L)
  subject_parameter <- paste0(effect_name, "[", subject, "]")
  summary <- .sft_parameter_summary(effect_draws, subject_parameter, hdi, rope)
  rownames(summary) <- NULL
  subject_summary <- summary
  rownames(subject_summary) <- subject
  prior_effect <- stats::rnorm(ndraws, prior_mean, prior_sd)
  prior_score <- prior_effect + stats::rnorm(ndraws, 0, 1 / sqrt(precision))
  posterior_predictive <- .sft_posterior_predictive(effect_draws, effect_hat,
                                                    precision)
  observation_weight <- precision / (precision + 1 / prior_sd^2)
  shrinkage_summary <- data.frame(
    subject = subject,
    mean_weight_on_individual = observation_weight,
    lower = observation_weight,
    upper = observation_weight,
    mean_weight_on_population = NA_real_,
    stringsAsFactors = FALSE
  )
  draws <- data.frame(effect_draws, check.names = FALSE)
  names(draws) <- subject_parameter
  list(effect_draws = effect_draws, subject_parameter = subject_parameter,
       post_mean = post_mean, post_var = post_var,
       summary = summary, subject_summary = subject_summary,
       posterior_predictive = posterior_predictive,
       prior_effect = prior_effect, prior_score = prior_score,
       shrinkage_summary = shrinkage_summary,
       shrinkage_draws = matrix(observation_weight, nrow = ndraws * chains, ncol = 1L),
       draws = draws, observation_weight = observation_weight,
       diagnostics = data.frame(parameter = subject_parameter, rhat = 1,
                                ess_bulk = ndraws * chains,
                                ess_tail = ndraws * chains,
                                stringsAsFactors = FALSE))
}


# Building blocks for the hierarchical normal-normal fit, shared by
# capacityGroup.bayes() and the multi-subject branch of micGroup.bayes(). It
# runs the conjugate Gibbs sampler (InvGamma) or a Stan variant, then produces
# the parameter summaries, predictive checks, shrinkage weights, and MCMC
# diagnostics. The caller sets the RNG seed and adds the model-specific labels.
.sft_normal_hierarchy_fit <- function(effect_hat, precision, subject, effect_name,
                                      method, ndraws, burnin, thin, chains,
                                      prior_mean, prior_sd, prior_shape, prior_rate,
                                      prior_tau_sd, seed, adapt_delta, max_treedepth,
                                      stan_control, hdi, rope,
                                      group = NULL, group_levels = NULL) {
  prior_predictive <- .sft_prior_predictive(
    effect_hat, precision, ndraws, method, prior_shape, prior_rate,
    prior_mean, prior_sd, prior_tau_sd, effect_name = effect_name
  )
  n <- length(effect_hat)
  if (is.null(group)) group <- rep(1L, n)
  group <- as.integer(group)
  n_group <- max(group)
  if (is.null(group_levels)) group_levels <- as.character(seq_len(n_group))
  mu_names <- if (n_group == 1L) "mu" else paste0("mu[", group_levels, "]")
  stan_fit <- NULL
  sampler_diagnostics <- NULL
  if (method == "InvGamma") {
    chain_draws <- lapply(seq_len(chains), function(chain) {
      .sft_hierarchical_normal_chain(effect_hat, precision, ndraws, burnin, thin,
                                     prior_mean, prior_sd, prior_shape, prior_rate,
                                     group = group)
    })
    chain_array <- array(NA_real_, dim = c(ndraws, chains, n_group + 1L + n))
    for (j in seq_len(chains)) {
      chain_array[, j, seq_len(n_group)] <- chain_draws[[j]]$mu
      chain_array[, j, n_group + 1L] <- chain_draws[[j]]$tau
      chain_array[, j, seq_len(n) + n_group + 1L] <- chain_draws[[j]]$theta
    }
  } else {
    stan <- .sft_run_stan_hierarchy(
      method, effect_hat, precision, ndraws, burnin, thin, chains,
      prior_mean, prior_sd, prior_tau_sd, seed, adapt_delta,
      max_treedepth, stan_control, effect_name = effect_name, group = group
    )
    chain_array <- stan$chain_array
    stan_fit <- stan$fit
    sampler_diagnostics <- stan$sampler_diagnostics
  }
  parameter_names <- c(mu_names, "tau", paste0(effect_name, "[", subject, "]"))
  dimnames(chain_array) <- list(NULL, paste0("chain", seq_len(chains)), parameter_names)
  total <- ndraws * chains
  mu_matrix <- matrix(chain_array[, , seq_len(n_group), drop = FALSE],
                      nrow = total, ncol = n_group,
                      dimnames = list(NULL, group_levels))
  tau_draws <- as.vector(chain_array[, , n_group + 1L, drop = TRUE])
  effect_draws <- matrix(chain_array[, , seq_len(n) + n_group + 1L, drop = FALSE],
                         nrow = total, ncol = n)
  mu_draws <- if (n_group == 1L) mu_matrix[, 1L] else mu_matrix
  parameter_draws <- cbind(mu_matrix, tau = tau_draws, effect_draws)
  summary <- .sft_parameter_summary(parameter_draws, parameter_names, hdi, rope)
  population <- summary[summary$parameter %in% c(mu_names, "tau"), , drop = FALSE]
  subject_summary <- summary[grepl(paste0("^", effect_name, "\\["), summary$parameter), , drop = FALSE]
  rownames(population) <- NULL
  rownames(subject_summary) <- subject
  posterior_predictive <- .sft_posterior_predictive(effect_draws, effect_hat, precision)
  shrinkage_draws <- .sft_shrinkage_draws(tau_draws, precision)
  shrinkage_summary <- .sft_shrinkage_summary(shrinkage_draws, subject, hdi)
  draws <- data.frame(mu_matrix, tau = tau_draws, effect_draws,
                      check.names = FALSE)
  names(draws) <- c(mu_names, "tau", paste0(effect_name, "[", subject, "]"))
  diagnostics <- .sft_mcmc_diagnostics(chain_array, parameter_names)
  list(mu_draws = mu_draws, mu_matrix = mu_matrix, tau_draws = tau_draws,
       effect_draws = effect_draws, group = group, group_levels = group_levels,
       mu_names = mu_names,
       parameter_names = parameter_names, summary = summary, population = population,
       subject_summary = subject_summary, posterior_predictive = posterior_predictive,
       shrinkage_draws = shrinkage_draws, shrinkage_summary = shrinkage_summary,
       draws = draws, diagnostics = diagnostics, prior_predictive = prior_predictive,
       stan_fit = stan_fit, sampler_diagnostics = sampler_diagnostics)
}


# ---- Bayes factors for the population / single-effect null -----------------
#
# Two complementary null tests, both computed from the draws the fit builders
# above already return. Which one is scientifically appropriate depends on the
# effect scale, so the public orchestrators choose at their call sites (an exact
# point null suits every score scale; an interval null suits the dimensionless
# relative MIC, eta, and reference-information phi effects). The generalised
# helpers themselves can compute either from any normal-normal fit.


# Point-null Savage-Dickey density ratio for an exact zero (the population mean
# mu = 0, or a single participant effect = 0). For the hierarchy it is
# Rao-Blackwellised: the exact Gibbs full conditional mu | theta, tau is
# Normal(m_mu, v_mu) with the same v_mu / m_mu used by the sampler, so averaging
# its density at the null over the posterior draws of (theta, tau) is a
# low-variance estimate of the marginal posterior density at zero. It is exact
# for the conjugate sampler and equally valid for the Stan variants, whose tau
# prior only reshapes the draws being averaged. BF01 > 1 favours the null. Its
# magnitude depends on the alternative prior_sd (the Jeffreys-Lindley
# sensitivity): a diffuse prior inflates the evidence for the null.
.sft_point_null_bf <- function(prior_mean, prior_sd, null = 0,
                               post_mean = NULL, post_var = NULL,
                               effect_draws = NULL, tau_draws = NULL,
                               group = NULL) {
  prior_density <- stats::dnorm(null, prior_mean, prior_sd)
  posterior_density <- if (!is.null(post_mean)) {
    stats::dnorm(null, post_mean, sqrt(post_var))
  } else {
    # Restrict to the units in this condition: the Gibbs full conditional for
    # mu_c involves only its own group's effects and count.
    if (!is.null(group)) effect_draws <- effect_draws[, group, drop = FALSE]
    n <- ncol(effect_draws)
    tau2 <- tau_draws^2
    v_mu <- 1 / (n / tau2 + 1 / prior_sd^2)
    m_mu <- v_mu * (rowSums(effect_draws) / tau2 + prior_mean / prior_sd^2)
    mean(stats::dnorm(null, m_mu, sqrt(v_mu)))
  }
  bf01 <- posterior_density / prior_density
  list(null = null,
       posterior_density_at_null = posterior_density,
       prior_density_at_null = prior_density,
       BF01 = bf01, BF10 = 1 / bf01,
       posterior_probability_null = bf01 / (1 + bf01),
       estimator = if (is.null(post_mean)) {
         "Savage-Dickey (Rao-Blackwellised)"
       } else "Savage-Dickey (closed form)")
}


# Interval-null Bayes factor: how much the data move probability into the
# practical-equivalence region |effect| < delta relative to the prior. Unlike
# the raw ROPE probability it credits the prior-to-posterior odds shift. The
# marginal prior on mu (and on a single analytic effect) is Normal(prior_mean,
# prior_sd), so the prior mass in the region is available in closed form.
#
# The posterior mass is a Monte Carlo proportion, so it saturates at 0 and 1
# once the region falls outside the draws. Reporting the raw proportion there
# would hand back BF01 = 0 or Inf as though it were an estimate. Instead the
# odds use a Jeffreys-smoothed proportion, (k + 1/2) / (n + 1), which keeps the
# Bayes factor finite and turns a saturated cell into an explicit resolution
# bound: the reported value is then the strongest evidence the draw count can
# support, flagged by mc_resolution_limited.
.sft_interval_null_bf <- function(effect_draws, delta, prior_mean, prior_sd) {
  x <- as.numeric(effect_draws)
  x <- x[is.finite(x)]
  n <- length(x)
  k <- sum(abs(x) < delta)
  posterior_null <- if (n) k / n else NA_real_
  smoothed <- (k + 0.5) / (n + 1)
  saturated <- n > 0L && (k == 0L || k == n)
  prior_null <- stats::pnorm(delta, prior_mean, prior_sd) -
    stats::pnorm(-delta, prior_mean, prior_sd)
  bf01 <- (smoothed / (1 - smoothed)) / (prior_null / (1 - prior_null))
  if (isTRUE(saturated)) {
    warning("The interval null captured ", if (k == 0L) "none" else "all",
            " of the ", n, " posterior draws, so the Bayes factor is limited by ",
            "Monte Carlo resolution; treat it as a bound and raise ndraws to ",
            "sharpen it.", call. = FALSE)
  }
  list(delta = delta,
       posterior_null_probability = posterior_null,
       prior_null_probability = prior_null,
       BF01 = bf01, BF10 = 1 / bf01,
       n_draws = n,
       mc_resolution_limited = saturated,
       estimator = "interval null, Jeffreys-smoothed posterior mass")
}


# Assemble the bayes_factor element for a single-subject analytic fit. The point
# null always applies; the interval null needs a practical-equivalence
# threshold (rope) and is skipped without one. point / interval let the caller
# match the tests to the effect scale.
.sft_analytic_bayes_factor <- function(fit, prior_mean, prior_sd, rope,
                                       point = TRUE, interval = TRUE) {
  out <- list()
  if (point) {
    out$point_null <- .sft_point_null_bf(prior_mean, prior_sd,
                                         post_mean = fit$post_mean,
                                         post_var = fit$post_var)
  }
  if (interval && !is.null(rope)) {
    out$interval_null <- .sft_interval_null_bf(fit$effect_draws[, 1L], rope,
                                               prior_mean, prior_sd)
  }
  if (length(out)) out else NULL
}


# Assemble the bayes_factor element for a hierarchical fit's population mean mu.
.sft_hierarchy_bayes_factor <- function(fit, prior_mean, prior_sd, rope,
                                        point = TRUE, interval = TRUE) {
  out <- list()
  if (point) {
    out$point_null <- .sft_point_null_bf(prior_mean, prior_sd,
                                         effect_draws = fit$effect_draws,
                                         tau_draws = fit$tau_draws)
  }
  if (interval && !is.null(rope)) {
    out$interval_null <- .sft_interval_null_bf(fit$mu_draws, rope,
                                               prior_mean, prior_sd)
  }
  if (length(out)) out else NULL
}


.sft_bayes_method <- function(method) {
  if (length(method) == 0L || is.null(method)) return("InvGamma")
  key <- tolower(gsub("[^A-Za-z]", "", as.character(method[[1L]])))
  switch(key,
         invgamma = "InvGamma", invgammagibbs = "InvGamma",
         ig = "InvGamma", gibbs = "InvGamma",
         halfnormal = "HalfNormal", hnorm = "HalfNormal",
         halfnormalstan = "HalfNormal", stanhalfnormal = "HalfNormal",
         centered = "Centered", centred = "Centered",
         centeredstan = "Centered", centredstan = "Centered",
         stancentered = "Centered", stancentred = "Centered",
         stop("method must be one of 'InvGamma', 'HalfNormal', or 'Centered'."))
}


# ---- Scale-aware default priors for the UCIP capacity models ---------------
#
# The priors are anchored on the *effect* scale theta = U / V, the Peto one-step
# estimator of the log hazard ratio for the UCIP contrast. theta is stable in
# location as trials accumulate and its precision V grows linearly with them, so
# a fixed prior width on theta is a genuine effect-size prior that the data
# eventually swamp.
#
# score anchor (theta, log hazard ratio):
#   mu_theta  ~ N(0, 0.5^2)         +/- 1 log-HR at 95%: hazard ratios in
#                                   roughly [0.37, 2.7], weakly informative on
#                                   the conventional log-HR scale
#   tau_theta ~ HalfNormal(0, 0.5)  between-subject heterogeneity of the same
#                                   order as the population effect
#   tau_theta^2 ~ InvGamma(2, 0.25) conjugate counterpart matching that second
#                                   moment (E[tau^2] = rate / (shape - 1))
#
# standardized (phi = sqrt(V_ref) * theta) is the *same model* reported on the
# reference-Cz scale, so its priors are the theta priors pushed through that
# linear map: every width is multiplied by sqrt(V_ref) (variances by V_ref).
# This is why V_ref is an argument here. Using a fixed width on phi instead --
# as an earlier version did -- implies a prior on theta of width
# prior_sd / sqrt(V_ref), which tightens at exactly the rate the likelihood
# sharpens and therefore never washes out, no matter how many trials are run.
#
# multiplicative (eta = log(A/B)) is already dimensionless: prior_sd =
# prior_tau_sd = 0.35 keeps the central one-SD prior capacity ratio near
# [0.70, 1.42]; the InvGamma second moment matches with prior_rate = 0.35^2.
.sft_default_priors <- function(score_method, V_ref) {
  if (score_method == "multiplicative") {
    return(list(prior_sd = 0.35, prior_tau_sd = 0.35,
                prior_shape = 2, prior_rate = 0.35^2))
  }
  theta_sd <- 0.5
  theta_tau_sd <- 0.5
  # On the phi scale every location/scale is multiplied by sqrt(V_ref); the
  # InvGamma rate is a variance and so scales by V_ref itself.
  scale <- if (score_method == "standardized") {
    if (is.null(V_ref) || !is.finite(V_ref) || V_ref <= 0) 1 else sqrt(V_ref)
  } else 1
  list(prior_sd = theta_sd * scale,
       prior_tau_sd = theta_tau_sd * scale,
       prior_shape = 2,
       prior_rate = theta_tau_sd^2 * scale^2)
}


# Fill any prior hyperparameter the caller left at NULL with the scale-aware
# default, leaving user-supplied values untouched. V_ref comes from the scored
# data and is retained in the signature for callers that may key defaults on it.
.sft_resolve_priors <- function(score_method, V_ref, prior_sd = NULL,
                                prior_tau_sd = NULL, prior_shape = NULL,
                                prior_rate = NULL) {
  defaults <- .sft_default_priors(score_method, V_ref)
  list(
    prior_sd = if (is.null(prior_sd)) defaults$prior_sd else prior_sd,
    prior_tau_sd = if (is.null(prior_tau_sd)) defaults$prior_tau_sd else prior_tau_sd,
    prior_shape = if (is.null(prior_shape)) defaults$prior_shape else prior_shape,
    prior_rate = if (is.null(prior_rate)) defaults$prior_rate else prior_rate
  )
}


.sft_validate_bayes_args <- function(ndraws, burnin, thin, chains,
                                     prior_shape, prior_rate, prior_mean,
                                     prior_sd, prior_tau_sd, hdi, rope,
                                     adapt_delta, max_treedepth) {
  if (ndraws < 100L || ndraws != as.integer(ndraws)) {
    stop("ndraws must be an integer >= 100.")
  }
  if (burnin < 0L || burnin != as.integer(burnin)) {
    stop("burnin must be a non-negative integer.")
  }
  if (thin < 1L || thin != as.integer(thin)) {
    stop("thin must be a positive integer.")
  }
  if (chains < 2L || chains != as.integer(chains)) {
    stop("chains must be an integer >= 2.")
  }
  # prior_shape / prior_rate / prior_sd / prior_tau_sd may be NULL, meaning
  # "resolve to a scale-aware default once the data are scored"; only
  # user-supplied (non-NULL) values are validated here.
  if ((!is.null(prior_shape) && (!is.finite(prior_shape) || prior_shape <= 0)) ||
      (!is.null(prior_rate) && (!is.finite(prior_rate) || prior_rate <= 0))) {
    stop("prior_shape and prior_rate must be positive finite values.")
  }
  if (!is.finite(prior_mean)) stop("prior_mean must be finite.")
  if (!is.null(prior_sd) && (!is.finite(prior_sd) || prior_sd <= 0)) {
    stop("prior_sd must be positive and finite.")
  }
  if (!is.null(prior_tau_sd) && (!is.finite(prior_tau_sd) || prior_tau_sd <= 0)) {
    stop("prior_tau_sd must be positive and finite.")
  }
  if (!is.finite(hdi) || hdi <= 0 || hdi >= 1) stop("hdi must lie in (0, 1).")
  if (!is.null(rope) && (!is.finite(rope) || rope <= 0)) {
    stop("rope must be NULL or a positive finite practical-equivalence threshold.")
  }
  if (!is.finite(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1) {
    stop("adapt_delta must lie in (0, 1).")
  }
  if (max_treedepth < 1L || max_treedepth != as.integer(max_treedepth)) {
    stop("max_treedepth must be a positive integer.")
  }
  invisible(TRUE)
}


# The effect a single UCIP bootstrap replicate contributes: the log hazard ratio
# theta for method = "score", the reference-Cz rescaling phi for method =
# "standardized", or the log-capacity ratio for the "multiplicative"
# (eta = log(A/B)) method. For the standardized score, v_ref is supplied by the
# original fitted participants and is held fixed across bootstrap replicates.
.sft_ucip_boot_effect <- function(rt, CR, stopping.rule, method,
                                  v_ref = NULL) {
  s <- .sft_ucip_score(rt, CR, stopping.rule = stopping.rule)
  switch(method,
         multiplicative = s$log_capacity,
         standardized = (s$numer / s$variance) * sqrt(v_ref),
         s$numer / s$variance)
}


# Within-condition, with-replacement bootstrap variance of the UCIP effect for
# one participant, returned as a precision (1 / variance). Trials are resampled
# within each condition, keeping each RT paired with its correctness, and the
# effect is recomputed on every replicate. Non-finite replicates (e.g. a
# resample with no correct trials in a channel for log-capacity) are dropped.
.sft_ucip_boot_precision <- function(rt, CR, stopping.rule, method, n_boot,
                                     v_ref = NULL) {
  reps <- vapply(seq_len(n_boot), function(b) {
    rb <- vector("list", length(rt))
    cb <- vector("list", length(rt))
    for (i in seq_along(rt)) {
      ni <- length(rt[[i]])
      idx <- if (ni > 0L) sample.int(ni, replace = TRUE) else integer(0)
      rb[[i]] <- rt[[i]][idx]
      cb[[i]] <- CR[[i]][idx]
    }
    tryCatch(.sft_ucip_boot_effect(rb, cb, stopping.rule, method, v_ref),
             error = function(e) NA_real_)
  }, numeric(1))
  reps <- reps[is.finite(reps)]
  if (length(reps) < 2L) return(NA_real_)
  v <- stats::var(reps)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  1 / v
}


.sft_ucip_score_data <- function(input, stopping.rule,
                                 method = c("score", "standardized",
                                            "multiplicative"),
                                 var_method = "analytic", n_boot = 2000L) {
  method <- .sft_score_method(method)
  var_method <- .sft_var_method(var_method)
  scores <- lapply(seq_along(input$RT), function(i) {
    .sft_ucip_score(input$RT[[i]], input$CR[[i]], stopping.rule = stopping.rule)
  })

  # A participant with a non-finite score or non-positive variance carries no
  # information for either model; the capacity model additionally needs an
  # estimable positive weighted numerator and denominator. Rather than aborting
  # the whole hierarchy on one such participant, drop them with a warning so the
  # remaining subjects can still be pooled. Erroring only when none survive.
  score_ok <- vapply(scores, function(x) {
    is.finite(x$numer) && is.finite(x$variance) && x$variance > 0
  }, logical(1))
  capacity_ok <- vapply(scores, function(x) {
    is.finite(x$log_capacity) && is.finite(x$log_capacity_precision) &&
      x$log_capacity_precision > 0
  }, logical(1))
  usable <- if (method == "multiplicative") score_ok & capacity_ok else score_ok

  if (!any(usable)) {
    stop(if (method == "multiplicative") {
      "No participant has an estimable positive weighted capacity effect with positive variance."
    } else {
      "No participant has a finite UCIP score with positive variance."
    }, call. = FALSE)
  }
  if (!all(usable)) {
    warning("Dropping ", sum(!usable), " participant(s) without an ",
            if (method == "multiplicative") "estimable capacity effect" else "estimable UCIP score",
            ": ", paste(input$subject[!usable], collapse = ", "), ".", call. = FALSE)
    scores <- scores[usable]
    input$RT <- input$RT[usable]
    input$CR <- input$CR[usable]
    input$subject <- input$subject[usable]
  }

  theta_hat <- vapply(scores, function(x) x$numer / x$variance, numeric(1))
  theta_precision <- vapply(scores, function(x) x$variance, numeric(1))
  eta_hat <- vapply(scores, function(x) x$log_capacity, numeric(1))
  eta_precision <- vapply(scores, function(x) x$log_capacity_precision, numeric(1))
  U <- vapply(scores, function(x) x$numer, numeric(1))
  V <- vapply(scores, function(x) x$variance, numeric(1))

  # The reference information V_ref (the median sampling variance) anchors both
  # the standardized phi scale and the information-scaled default priors. It is
  # computed after degenerate-subject dropping so it is based on exactly the
  # participants entering the fit. For one participant this makes phi_hat = Cz
  # and its analytic precision equal to one by construction. It is computed for
  # every method even though only "standardized" rescales the effect by it.
  V_ref <- stats::median(V)
  phi_hat <- theta_hat * sqrt(V_ref)

  effect_hat <- switch(method, multiplicative = eta_hat,
                       standardized = phi_hat, theta_hat)
  precision <- switch(method, multiplicative = eta_precision,
                      standardized = V / V_ref, V)
  effect_name <- switch(method, multiplicative = "eta",
                        standardized = "phi", "theta")

  # An optional within-subject bootstrap variance for the active effect scale.
  # It replaces the closed-form precision used by the hierarchy; the analytic
  # value is retained for any subject whose bootstrap variance is undefined.
  analytic_se <- 1 / sqrt(precision)
  used_se <- analytic_se
  if (var_method == "bootstrap") {
    boot_precision <- vapply(seq_along(input$RT), function(i) {
      cr_i <- .sft_cr_list(input$RT[[i]], input$CR[[i]])
      .sft_ucip_boot_precision(lapply(input$RT[[i]], as.numeric), cr_i,
                               stopping.rule, method, n_boot, V_ref)
    }, numeric(1))
    undefined <- !is.finite(boot_precision)
    boot_precision[undefined] <- precision[undefined]
    precision <- boot_precision
    used_se <- 1 / sqrt(precision)
  }

  score_df <- if (method == "score") {
    data.frame(
      subject = input$subject,
      U = U,
      V = V,
      theta_hat = theta_hat,
      hazard_ratio_hat = exp(theta_hat),
      se = used_se,
      se_analytic = analytic_se,
      Cz = vapply(scores, function(x) unname(x$statistic), numeric(1)),
      stringsAsFactors = FALSE
    )
  } else if (method == "standardized") {
    data.frame(
      subject = input$subject,
      U = U,
      V = V,
      V_ref = rep(V_ref, length(V)),
      theta_hat = theta_hat,
      phi_hat = phi_hat,
      se = used_se,
      se_analytic = analytic_se,
      Cz = vapply(scores, function(x) unname(x$statistic), numeric(1)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      subject = input$subject,
      U = U,
      V = V,
      Cz = vapply(scores, function(x) unname(x$statistic), numeric(1)),
      capacity_numerator = vapply(
        scores, function(x) x$capacity_numerator, numeric(1)
      ),
      capacity_denominator = vapply(
        scores, function(x) x$capacity_denominator, numeric(1)
      ),
      capacity_ratio_hat = exp(eta_hat),
      eta_hat = eta_hat,
      eta_se = used_se,
      eta_se_analytic = analytic_se,
      score_effect_hat = U / V,
      score_effect_se = 1 / sqrt(V),
      stringsAsFactors = FALSE
    )
  }
  if (var_method == "analytic") {
    score_df$se_analytic <- NULL
    score_df$eta_se_analytic <- NULL
  }

  list(scores = scores, subject = input$subject,
       theta_hat = theta_hat, phi_hat = phi_hat, eta_hat = eta_hat,
       V_ref = V_ref,
       precision = precision, theta_precision = theta_precision,
       eta_precision = eta_precision, effect_hat = effect_hat,
       effect_name = effect_name, method = method, var_method = var_method,
       score = score_df)
}


.sft_bayes_seed <- function(seed) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  if (!is.null(seed)) set.seed(seed)
  old_seed
}


.sft_restore_bayes_seed <- function(old_seed) {
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }
  invisible(NULL)
}


.sft_parameter_summary <- function(parameter_draws, parameter_names, hdi, rope) {
  summaries <- lapply(seq_len(ncol(parameter_draws)), function(j) {
    x <- parameter_draws[, j]
    interval <- .sft_hdi(x, hdi)
    # tau is a non-negative between-subject scale, not a signed effect: its sign
    # probabilities are 1 and 0 by construction and a ROPE around zero is not a
    # practical-equivalence statement, so those columns are left NA rather than
    # reported as though they were inferential.
    signed <- !identical(parameter_names[j], "tau")
    rope_probability <- if (is.null(rope) || !signed) NA_real_ else
      mean(abs(x) <= rope)
    data.frame(parameter = parameter_names[j], mean = mean(x),
               lower = interval[["lower"]], upper = interval[["upper"]],
               posterior_positive = if (signed) mean(x > 0) else NA_real_,
               posterior_negative = if (signed) mean(x < 0) else NA_real_,
               posterior_rope = rope_probability,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, summaries)
}


.sft_split_chains <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (n < 4L) return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  half <- floor(n / 2L)
  out <- matrix(NA_real_, nrow = half, ncol = ncol(x) * 2L)
  for (j in seq_len(ncol(x))) {
    out[, 2L * j - 1L] <- x[seq_len(half), j]
    out[, 2L * j] <- x[n - half + seq_len(half), j]
  }
  out
}


.sft_fallback_ess <- function(x) {
  x <- as.numeric(x[is.finite(x)])
  n <- length(x)
  if (n < 3L || stats::sd(x) == 0) return(as.numeric(n))
  lag_max <- min(n - 1L, 1000L)
  ac <- stats::acf(x, plot = FALSE, lag.max = lag_max, demean = TRUE)$acf[, 1L, 1L]
  rho <- ac[-1L]
  pairs <- floor(length(rho) / 2L)
  if (!pairs) return(as.numeric(n))
  gamma <- rho[2L * seq_len(pairs) - 1L] + rho[2L * seq_len(pairs)]
  first_nonpositive <- which(gamma <= 0)
  last <- if (length(first_nonpositive)) min(first_nonpositive) - 1L else length(gamma)
  if (last < 1L) return(as.numeric(n))
  tau_int <- max(1, 1 + 2 * sum(gamma[seq_len(last)]))
  min(as.numeric(n), as.numeric(n) / tau_int)
}


.sft_manual_diagnostics <- function(chain_array, parameter_names) {
  out <- lapply(seq_along(parameter_names), function(j) {
    chains <- chain_array[, , j, drop = TRUE]
    split <- .sft_split_chains(chains)
    if (!length(split)) {
      return(data.frame(parameter = parameter_names[j], rhat = NA_real_,
                        ess_bulk = NA_real_, ess_tail = NA_real_,
                        stringsAsFactors = FALSE))
    }
    within <- mean(apply(split, 2L, stats::var))
    between <- nrow(split) * stats::var(colMeans(split))
    rhat <- if (is.finite(within) && within > 0) {
      sqrt(((nrow(split) - 1) * within + between) /
             (nrow(split) * within))
    } else if (is.finite(between) && between == 0) 1 else NA_real_
    flat <- as.vector(split)
    bulk <- stats::qnorm((rank(flat, ties.method = "average") - .5) / length(flat))
    q <- stats::quantile(flat, c(.05, .95), names = FALSE, type = 8)
    tail <- as.numeric(flat <= q[1L] | flat >= q[2L])
    data.frame(parameter = parameter_names[j], rhat = rhat,
               ess_bulk = .sft_fallback_ess(bulk),
               ess_tail = .sft_fallback_ess(tail),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}


.sft_mcmc_diagnostics <- function(chain_array, parameter_names) {
  if (requireNamespace("posterior", quietly = TRUE)) {
    dimnames(chain_array) <- list(NULL, paste0("chain", seq_len(dim(chain_array)[2L])),
                                  parameter_names)
    ans <- tryCatch({
      draws <- posterior::as_draws_array(chain_array)
      stats <- posterior::summarise_draws(
        draws, posterior::rhat, posterior::ess_bulk, posterior::ess_tail
      )
      data.frame(parameter = as.character(stats$variable),
                 rhat = as.numeric(stats[[2L]]),
                 ess_bulk = as.numeric(stats[[3L]]),
                 ess_tail = as.numeric(stats[[4L]]),
                 stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    if (!is.null(ans)) return(ans)
  }
  .sft_manual_diagnostics(chain_array, parameter_names)
}


.sft_prior_predictive <- function(effect_hat, precision, ndraws, method,
                                  prior_shape, prior_rate, prior_mean,
                                  prior_sd, prior_tau_sd,
                                  effect_name = "phi") {
  prior_mu <- stats::rnorm(ndraws, prior_mean, prior_sd)
  if (method == "InvGamma") {
    prior_tau2 <- 1 / stats::rgamma(ndraws, shape = prior_shape, rate = prior_rate)
    prior_tau <- sqrt(prior_tau2)
  } else {
    prior_tau <- abs(stats::rnorm(ndraws, 0, prior_tau_sd))
    prior_tau2 <- prior_tau^2
  }
  n <- length(effect_hat)
  prior_effect <- matrix(
    stats::rnorm(ndraws * n,
                 mean = rep(prior_mu, each = n),
                 sd = rep(prior_tau, each = n)),
    nrow = ndraws, ncol = n, byrow = TRUE
  )
  prior_score <- prior_effect + sweep(
    matrix(stats::rnorm(length(prior_effect)), nrow = ndraws),
    2L, 1 / sqrt(precision), "*"
  )
  out <- list(mu = prior_mu, tau = prior_tau, score = prior_score)
  out[[effect_name]] <- prior_effect
  out[c("mu", "tau", effect_name, "score")]
}


.sft_posterior_predictive <- function(effect_draws, effect_hat, precision) {
  posterior_score <- effect_draws + sweep(
    matrix(stats::rnorm(length(effect_draws)), nrow = nrow(effect_draws)),
    2L, 1 / sqrt(precision), "*"
  )
  observed_mean <- mean(effect_hat)
  observed_sd <- if (length(effect_hat) > 1L) stats::sd(effect_hat) else 0
  predictive_mean <- rowMeans(posterior_score)
  predictive_sd <- if (ncol(posterior_score) > 1L) {
    apply(posterior_score, 1L, stats::sd)
  } else {
    rep(0, nrow(posterior_score))
  }
  list(
    draws = posterior_score,
    observed = effect_hat,
    discrepancy = data.frame(
      measure = c("mean", "sd"),
      observed = c(observed_mean, observed_sd),
      predictive_mean = c(mean(predictive_mean), mean(predictive_sd)),
      posterior_tail_probability = c(
        mean(predictive_mean >= observed_mean),
        mean(predictive_sd >= observed_sd)
      ),
      stringsAsFactors = FALSE
    )
  )
}


.sft_shrinkage_draws <- function(tau_draws, precision) {
  sweep(
    matrix(tau_draws^2, nrow = length(tau_draws), ncol = length(precision)),
    2L, precision,
    function(tau2, information) information * tau2 / (1 + information * tau2)
  )
}


.sft_shrinkage_summary <- function(shrinkage_draws, subjects, hdi) {
  probs <- c((1 - hdi) / 2, 1 - (1 - hdi) / 2)
  data.frame(
    subject = subjects,
    mean_weight_on_individual = colMeans(shrinkage_draws),
    lower = apply(shrinkage_draws, 2L, stats::quantile, probs = probs[1L], names = FALSE),
    upper = apply(shrinkage_draws, 2L, stats::quantile, probs = probs[2L], names = FALSE),
    mean_weight_on_population = 1 - colMeans(shrinkage_draws),
    stringsAsFactors = FALSE
  )
}


.sft_stan_model_cache <- new.env(parent = emptyenv())


# The generated model carries a condition index so one population mean mu[c] is
# estimated per condition over a shared tau. With C = 1 and every cond[i] = 1
# this is exactly the single-mean model.
.sft_stan_code <- function(method, effect_name = "phi",
                           precision_name = "W") {
  effect_hat_name <- paste0(effect_name, "_hat")
  data_block <- paste0("data {
  int<lower=1> N;
  int<lower=1> C;
  int<lower=1, upper=C> cond[N];
  vector[N] ", effect_hat_name, ";
  vector<lower=0>[N] ", precision_name, ";
  real prior_mean;
  real<lower=0> prior_sd;
  real<lower=0> prior_tau_sd;
}
")
  likelihood <- paste0("  for (i in 1:N)
    ", effect_hat_name, "[i] ~ normal(", effect_name, "[i], 1 / sqrt(",
                       precision_name, "[i]));\n}")
  if (method == "HalfNormal") {
    return(paste0(data_block, "parameters {
  vector[C] mu;
  real<lower=0> tau;
  vector[N] z;
}
transformed parameters {
  vector[N] ", effect_name, ";
  for (i in 1:N)
    ", effect_name, "[i] = mu[cond[i]] + tau * z[i];
}
model {
  mu ~ normal(prior_mean, prior_sd);
  tau ~ normal(0, prior_tau_sd);
  z ~ normal(0, 1);
", likelihood))
  }
  paste0(data_block, "parameters {
  vector[C] mu;
  real<lower=0> tau;
  vector[N] ", effect_name, ";
}
model {
  mu ~ normal(prior_mean, prior_sd);
  tau ~ normal(0, prior_tau_sd);
  for (i in 1:N)
    ", effect_name, "[i] ~ normal(mu[cond[i]], tau);
", likelihood)
}


.sft_stan_sampler_diagnostics <- function(fit, max_treedepth) {
  params <- tryCatch(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                     error = function(e) NULL)
  if (is.null(params)) return(NULL)
  per_chain <- lapply(seq_along(params), function(i) {
    x <- params[[i]]
    divergent <- if ("divergent__" %in% colnames(x)) sum(x[, "divergent__"] > 0) else NA_real_
    treedepth <- if ("treedepth__" %in% colnames(x)) sum(x[, "treedepth__"] >= max_treedepth) else NA_real_
    data.frame(chain = i, divergences = divergent,
               max_treedepth_hits = treedepth, stringsAsFactors = FALSE)
  })
  per_chain <- do.call(rbind, per_chain)
  list(per_chain = per_chain,
       divergences = sum(per_chain$divergences, na.rm = TRUE),
       max_treedepth_hits = sum(per_chain$max_treedepth_hits, na.rm = TRUE))
}


.sft_run_stan_hierarchy <- function(method, effect_hat, precision, ndraws,
                                    burnin, thin, chains, prior_mean,
                                    prior_sd, prior_tau_sd, seed,
                                    adapt_delta, max_treedepth, stan_control,
                                    effect_name = "phi", group = NULL) {
  if (is.null(group)) group <- rep(1L, length(effect_hat))
  group <- as.integer(group)
  n_group <- max(group)
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("method='", method, "' requires the optional 'rstan' package.")
  }
  precision_name <- switch(effect_name, eta = "P", phi = "W", "V")
  code <- .sft_stan_code(method, effect_name, precision_name)
  # prior_tau_sd enters the model as *data*, not as a literal in the generated
  # code, so it must not discriminate the compiled-model cache: including it
  # forced a fresh (identical) compile for every distinct prior width.
  key <- paste0(method, "_", effect_name)
  model <- if (exists(key, envir = .sft_stan_model_cache, inherits = FALSE)) {
    get(key, envir = .sft_stan_model_cache, inherits = FALSE)
  } else {
    fitted_model <- rstan::stan_model(model_code = code,
                                       model_name = paste0("ucip_bayes_", method),
                                       verbose = FALSE)
    assign(key, fitted_model, envir = .sft_stan_model_cache)
    fitted_model
  }
  control <- utils::modifyList(
    list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    stan_control
  )
  stan_seed <- if (is.null(seed)) -1L else seed
  stan_data <- list(N = length(effect_hat), C = n_group, cond = group,
                    prior_mean = prior_mean, prior_sd = prior_sd,
                    prior_tau_sd = prior_tau_sd)
  stan_data[[paste0(effect_name, "_hat")]] <- effect_hat
  stan_data[[precision_name]] <- precision
  fit <- rstan::sampling(
    model,
    data = stan_data,
    chains = chains, iter = burnin + ndraws * thin, warmup = burnin,
    thin = thin, seed = stan_seed, refresh = 0, cores = 1L,
    control = control
  )
  raw <- rstan::extract(fit, pars = c("mu", "tau", effect_name),
                        permuted = FALSE, inc_warmup = FALSE)
  variables <- dimnames(raw)[[3L]]
  effect_idx <- grep(paste0("^", effect_name, "\\["), variables)
  if (length(effect_idx) != length(effect_hat)) {
    stop("Stan did not return the expected participant effects.")
  }
  effect_idx <- effect_idx[order(as.integer(gsub("[^0-9]", "", variables[effect_idx])))]
  mu_idx <- grep("^mu(\\[|$)", variables)
  mu_idx <- mu_idx[order(as.integer(gsub("[^0-9]", "", paste0(variables[mu_idx], "1"))))]
  if (length(mu_idx) != n_group) {
    stop("Stan did not return the expected condition means.")
  }
  chain_array <- array(NA_real_,
                       dim = c(dim(raw)[1L], dim(raw)[2L],
                               n_group + 1L + length(effect_hat)))
  for (j in seq_len(n_group)) chain_array[, , j] <- raw[, , mu_idx[j]]
  chain_array[, , n_group + 1L] <- raw[, , match("tau", variables)]
  for (j in seq_along(effect_idx)) {
    chain_array[, , n_group + 1L + j] <- raw[, , effect_idx[j]]
  }
  dimnames(chain_array) <- list(
    NULL, NULL,
    c(if (n_group == 1L) "mu" else paste0("mu[", seq_len(n_group), "]"), "tau",
      paste0(effect_name, "[", seq_along(effect_hat), "]")))
  list(chain_array = chain_array, fit = fit,
       sampler_diagnostics = .sft_stan_sampler_diagnostics(fit, max_treedepth))
}
