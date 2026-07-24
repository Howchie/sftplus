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
# effect being pooled (UCIP score/Cz, log-capacity, or the mean interaction
# contrast): every caller supplies a per-subject effect estimate and its
# sampling precision, so the updates below are unchanged across models.
.sft_hierarchical_normal_chain <- function(theta_hat, precision, ndraws, burnin,
                                           thin, prior_mean, prior_sd,
                                           tau2_prior_shape, tau2_prior_rate) {
  n <- length(theta_hat)
  niter <- burnin + ndraws * thin
  mu_draws <- numeric(ndraws)
  tau2_draws <- numeric(ndraws)
  theta_draws <- matrix(NA_real_, nrow = ndraws, ncol = n)
  # Over-dispersed per-chain starting values so split-R-hat reflects genuine
  # mixing rather than a shared deterministic start. Each chain advances the
  # global RNG, so tau2 (drawn from its prior) and the jittered mu differ across
  # chains; theta is resampled at the top of the loop, so its start is immaterial.
  tau2 <- 1 / stats::rgamma(1L, shape = tau2_prior_shape, rate = tau2_prior_rate)
  mu <- stats::rnorm(1L, stats::weighted.mean(theta_hat, precision), sqrt(tau2))
  theta <- theta_hat
  save_i <- 0L

  for (iter in seq_len(niter)) {
    theta_var <- 1 / (precision + 1 / tau2)
    theta_mean <- theta_var * (precision * theta_hat + mu / tau2)
    theta <- stats::rnorm(n, theta_mean, sqrt(theta_var))

    mu_var <- 1 / (n / tau2 + 1 / prior_sd^2)
    mu_mean <- mu_var * (sum(theta) / tau2 + prior_mean / prior_sd^2)
    mu <- stats::rnorm(1L, mu_mean, sqrt(mu_var))

    # tau^2 has an inverse-Gamma full conditional. Sampling its reciprocal
    # as Gamma keeps the update stable and avoids a parameter at zero.
    tau2 <- 1 / stats::rgamma(
      1L, shape = tau2_prior_shape + n / 2,
      rate = tau2_prior_rate + sum((theta - mu)^2) / 2
    )

    if (iter > burnin && (iter - burnin) %% thin == 0L) {
      save_i <- save_i + 1L
      mu_draws[save_i] <- mu
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
                                      stan_control, hdi, rope) {
  prior_predictive <- .sft_prior_predictive(
    effect_hat, precision, ndraws, method, prior_shape, prior_rate,
    prior_mean, prior_sd, prior_tau_sd, effect_name = effect_name
  )
  n <- length(effect_hat)
  stan_fit <- NULL
  sampler_diagnostics <- NULL
  if (method == "InvGamma") {
    chain_draws <- lapply(seq_len(chains), function(chain) {
      .sft_hierarchical_normal_chain(effect_hat, precision, ndraws, burnin, thin,
                                     prior_mean, prior_sd, prior_shape, prior_rate)
    })
    chain_array <- array(NA_real_, dim = c(ndraws, chains, 2L + n))
    for (j in seq_len(chains)) {
      chain_array[, j, 1L] <- chain_draws[[j]]$mu
      chain_array[, j, 2L] <- chain_draws[[j]]$tau
      chain_array[, j, seq_len(n) + 2L] <- chain_draws[[j]]$theta
    }
  } else {
    stan <- .sft_run_stan_hierarchy(
      method, effect_hat, precision, ndraws, burnin, thin, chains,
      prior_mean, prior_sd, prior_tau_sd, seed, adapt_delta,
      max_treedepth, stan_control, effect_name = effect_name
    )
    chain_array <- stan$chain_array
    stan_fit <- stan$fit
    sampler_diagnostics <- stan$sampler_diagnostics
  }
  parameter_names <- c("mu", "tau", paste0(effect_name, "[", subject, "]"))
  dimnames(chain_array) <- list(NULL, paste0("chain", seq_len(chains)), parameter_names)
  mu_chain <- chain_array[, , 1L, drop = TRUE]
  tau_chain <- chain_array[, , 2L, drop = TRUE]
  effect_array <- chain_array[, , seq_len(n) + 2L, drop = FALSE]
  mu_draws <- as.vector(mu_chain)
  tau_draws <- as.vector(tau_chain)
  effect_draws <- matrix(effect_array, nrow = ndraws * chains, ncol = n)
  parameter_draws <- cbind(mu = mu_draws, tau = tau_draws, effect_draws)
  summary <- .sft_parameter_summary(parameter_draws, parameter_names, hdi, rope)
  population <- summary[summary$parameter %in% c("mu", "tau"), , drop = FALSE]
  subject_summary <- summary[grepl(paste0("^", effect_name, "\\["), summary$parameter), , drop = FALSE]
  rownames(population) <- NULL
  rownames(subject_summary) <- subject
  posterior_predictive <- .sft_posterior_predictive(effect_draws, effect_hat, precision)
  shrinkage_draws <- .sft_shrinkage_draws(tau_draws, precision)
  shrinkage_summary <- .sft_shrinkage_summary(shrinkage_draws, subject, hdi)
  draws <- data.frame(mu = mu_draws, tau = tau_draws, effect_draws,
                      check.names = FALSE)
  names(draws)[-(1:2)] <- paste0(effect_name, "[", subject, "]")
  diagnostics <- .sft_mcmc_diagnostics(chain_array, parameter_names)
  list(mu_draws = mu_draws, tau_draws = tau_draws, effect_draws = effect_draws,
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
# point null suits the raw score and the multiplicative/eta effect; an interval
# null suits the dimensionless relative MIC and the eta effect). The generalised
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
                               effect_draws = NULL, tau_draws = NULL) {
  prior_density <- stats::dnorm(null, prior_mean, prior_sd)
  posterior_density <- if (!is.null(post_mean)) {
    stats::dnorm(null, post_mean, sqrt(post_var))
  } else {
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
.sft_interval_null_bf <- function(effect_draws, delta, prior_mean, prior_sd) {
  posterior_null <- mean(abs(effect_draws) < delta)
  prior_null <- stats::pnorm(delta, prior_mean, prior_sd) -
    stats::pnorm(-delta, prior_mean, prior_sd)
  posterior_odds <- posterior_null / (1 - posterior_null)
  prior_odds <- prior_null / (1 - prior_null)
  bf01 <- posterior_odds / prior_odds
  list(delta = delta,
       posterior_null_probability = posterior_null,
       prior_null_probability = prior_null,
       BF01 = bf01, BF10 = 1 / bf01)
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
  if (!is.finite(prior_shape) || prior_shape <= 0 ||
      !is.finite(prior_rate) || prior_rate <= 0) {
    stop("prior_shape and prior_rate must be positive finite values.")
  }
  if (!is.finite(prior_mean) || !is.finite(prior_sd) || prior_sd <= 0) {
    stop("prior_mean must be finite and prior_sd must be positive.")
  }
  if (!is.finite(prior_tau_sd) || prior_tau_sd <= 0) {
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


# The effect a single UCIP bootstrap replicate contributes: the score-effect
# U/V for method = "score" or the log-capacity ratio for the "multiplicative"
# (eta = log(A/B)) method.
.sft_ucip_boot_effect <- function(rt, CR, stopping.rule, method) {
  s <- .sft_ucip_score(rt, CR, stopping.rule = stopping.rule)
  if (method == "multiplicative") s$log_capacity else s$numer / s$variance
}


# Within-condition, with-replacement bootstrap variance of the UCIP effect for
# one participant, returned as a precision (1 / variance). Trials are resampled
# within each condition, keeping each RT paired with its correctness, and the
# effect is recomputed on every replicate. Non-finite replicates (e.g. a
# resample with no correct trials in a channel for log-capacity) are dropped.
.sft_ucip_boot_precision <- function(rt, CR, stopping.rule, method, n_boot) {
  reps <- vapply(seq_len(n_boot), function(b) {
    rb <- vector("list", length(rt))
    cb <- vector("list", length(rt))
    for (i in seq_along(rt)) {
      ni <- length(rt[[i]])
      idx <- if (ni > 0L) sample.int(ni, replace = TRUE) else integer(0)
      rb[[i]] <- rt[[i]][idx]
      cb[[i]] <- CR[[i]][idx]
    }
    tryCatch(.sft_ucip_boot_effect(rb, cb, stopping.rule, method),
             error = function(e) NA_real_)
  }, numeric(1))
  reps <- reps[is.finite(reps)]
  if (length(reps) < 2L) return(NA_real_)
  v <- stats::var(reps)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  1 / v
}


.sft_ucip_score_data <- function(input, stopping.rule,
                                 method = c("score", "multiplicative"),
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

  effect_hat <- if (method == "multiplicative") eta_hat else theta_hat
  precision <- if (method == "multiplicative") eta_precision else theta_precision
  effect_name <- if (method == "multiplicative") "eta" else "theta"

  # An optional within-subject bootstrap variance for the active effect scale.
  # It replaces the closed-form precision used by the hierarchy; the analytic
  # value is retained for any subject whose bootstrap variance is undefined.
  analytic_se <- 1 / sqrt(precision)
  used_se <- analytic_se
  if (var_method == "bootstrap") {
    boot_precision <- vapply(seq_along(input$RT), function(i) {
      cr_i <- .sft_cr_list(input$RT[[i]], input$CR[[i]])
      .sft_ucip_boot_precision(lapply(input$RT[[i]], as.numeric), cr_i,
                               stopping.rule, method, n_boot)
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
       theta_hat = theta_hat, eta_hat = eta_hat,
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
    rope_probability <- if (is.null(rope)) NA_real_ else mean(abs(x) <= rope)
    data.frame(parameter = parameter_names[j], mean = mean(x),
               lower = interval[["lower"]], upper = interval[["upper"]],
               posterior_positive = mean(x > 0),
               posterior_negative = mean(x < 0),
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
  keep <- which(gamma > 0)
  if (!length(keep)) return(as.numeric(n))
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
                                  effect_name = "theta") {
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


.sft_stan_code <- function(method, effect_name = "theta",
                           precision_name = "V") {
  effect_hat_name <- paste0(effect_name, "_hat")
  if (method == "HalfNormal") {
    return(paste0("data {
  int<lower=1> N;
  vector[N] ", effect_hat_name, ";
  vector<lower=0>[N] ", precision_name, ";
  real prior_mean;
  real<lower=0> prior_sd;
  real<lower=0> prior_tau_sd;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[N] z;
}
transformed parameters {
  vector[N] ", effect_name, ";
  ", effect_name, " = mu + tau * z;
}
model {
  mu ~ normal(prior_mean, prior_sd);
  tau ~ normal(0, prior_tau_sd);
  z ~ normal(0, 1);
  for (i in 1:N)
    ", effect_hat_name, "[i] ~ normal(", effect_name, "[i], 1 / sqrt(", precision_name, "[i]));", "}"))
  }
  paste0("data {
  int<lower=1> N;
  vector[N] ", effect_hat_name, ";
  vector<lower=0>[N] ", precision_name, ";
  real prior_mean;
  real<lower=0> prior_sd;
  real<lower=0> prior_tau_sd;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[N] ", effect_name, ";
}
model {
  mu ~ normal(prior_mean, prior_sd);
  tau ~ normal(0, prior_tau_sd);
  ", effect_name, " ~ normal(mu, tau);
  for (i in 1:N)
    ", effect_hat_name, "[i] ~ normal(", effect_name, "[i], 1 / sqrt(", precision_name, "[i]));", "}")
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
                                    effect_name = "theta") {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("method='", method, "' requires the optional 'rstan' package.")
  }
  precision_name <- if (effect_name == "eta") "P" else "V"
  code <- .sft_stan_code(method, effect_name, precision_name)
  key <- paste0(method, "_", as.character(prior_tau_sd), "_", effect_name)
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
  stan_data <- list(N = length(effect_hat),
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
  chain_array <- array(NA_real_,
                       dim = c(dim(raw)[1L], dim(raw)[2L], 2L + length(effect_hat)))
  chain_array[, , 1L] <- raw[, , match("mu", variables)]
  chain_array[, , 2L] <- raw[, , match("tau", variables)]
  for (j in seq_along(effect_idx)) chain_array[, , 2L + j] <- raw[, , effect_idx[j]]
  dimnames(chain_array) <- list(NULL, NULL,
                                c("mu", "tau", paste0(effect_name, "[", seq_along(effect_hat), "]")))
  list(chain_array = chain_array, fit = fit,
       sampler_diagnostics = .sft_stan_sampler_diagnostics(fit, max_treedepth))
}

