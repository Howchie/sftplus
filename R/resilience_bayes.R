# Bayesian companions to the resilience coefficient: resilience.bayes() for one
# participant and resilienceGroup.bayes() for a hierarchy. Both reuse the
# normal-normal engine in ucip_bayes_engine.R.
#
# The frequentist resilience test contrasts Nelson-Aalen cumulative hazards
# against the UCIP-OR benchmark,
#
#   delta(t) = H_AB(t) - H_A(t) - H_B(t),
#
# which is zero for an unlimited-capacity independent parallel OR system. Both
# effect scales below are exactly zero under that benchmark:
#
#   "difference"     delta = H_AB - H_A - H_B at the evaluation horizon, in
#                    cumulative-hazard units, with the Nelson-Aalen variance sum.
#   "multiplicative" psi = log(H_AB / (H_A + H_B)), the log resilience ratio, a
#                    dimensionless quantity with a delta-method variance.
#
# psi is the default for pooling because it is horizon-robust: it compares two
# cumulative hazards to each other rather than reporting an absolute level, so
# participants observed over different time ranges remain comparable. delta is
# on an absolute scale and is only comparable across participants evaluated at
# the same horizon -- see the `at` argument.


.sft_resilience_method <- function(method) {
  key <- tolower(gsub("[^a-z]", "", as.character(method[[1L]])))
  switch(key,
         multiplicative = "multiplicative", logratio = "multiplicative",
         ratio = "multiplicative", psi = "multiplicative",
         difference = "difference", diff = "difference",
         contrast = "difference", delta = "difference",
         stop("score_method must be one of 'multiplicative' or 'difference'."))
}


# One participant's resilience effect at a horizon. `at` selects the evaluation
# time; NULL uses that participant's own largest finite RT, which reproduces
# resilience.test() but is participant-specific.
.sft_resilience_score <- function(RT, CR = NULL, at = NULL) {
  if (!is.list(RT) || length(RT) < 3L) {
    stop("RT must contain AB, A, and B conditions.", call. = FALSE)
  }
  RT <- RT[1:3]
  CR <- .sft_cr_list(RT, CR)
  times <- .sft_finite_times(RT)
  t_star <- if (is.null(at)) times[[length(times)]] else as.numeric(at)[[1L]]
  h <- lapply(seq_len(3L), function(i) estimateNAH(RT[[i]], CR[[i]]))
  H <- vapply(h, function(x) x$H(t_star), numeric(1))
  V <- vapply(h, function(x) x$Var(t_star), numeric(1))
  numerator <- H[[1L]]
  denominator <- H[[2L]] + H[[3L]]
  denominator_var <- V[[2L]] + V[[3L]]
  delta <- numerator - denominator
  delta_var <- sum(V)
  ok_ratio <- is.finite(numerator) && is.finite(denominator) &&
    numerator > 0 && denominator > 0
  psi <- if (ok_ratio) log(numerator / denominator) else NA_real_
  psi_var <- if (ok_ratio) {
    V[[1L]] / numerator^2 + denominator_var / denominator^2
  } else NA_real_
  z <- if (is.finite(delta_var) && delta_var > 0) delta / sqrt(delta_var) else NA_real_
  list(time = t_star, max_time = times[[length(times)]],
       H_AB = H[[1L]], H_A = H[[2L]], H_B = H[[3L]],
       delta = delta, delta_variance = delta_var,
       resilience_ratio = if (ok_ratio) numerator / denominator else NA_real_,
       psi = psi, psi_variance = psi_var, z = z)
}


# The effect one bootstrap replicate contributes, on whichever scale is active.
.sft_resilience_boot_effect <- function(rt, CR, at, method) {
  s <- .sft_resilience_score(rt, CR, at = at)
  if (method == "multiplicative") s$psi else s$delta
}


# Within-condition resampling variance of the resilience effect, as a precision.
# The horizon is held fixed across replicates so every replicate estimates the
# same quantity.
.sft_resilience_boot_precision <- function(rt, CR, at, method, n_boot) {
  reps <- vapply(seq_len(n_boot), function(b) {
    rb <- vector("list", length(rt)); cb <- vector("list", length(rt))
    for (i in seq_along(rt)) {
      ni <- length(rt[[i]])
      idx <- if (ni > 0L) sample.int(ni, replace = TRUE) else integer(0)
      rb[[i]] <- rt[[i]][idx]; cb[[i]] <- CR[[i]][idx]
    }
    tryCatch(.sft_resilience_boot_effect(rb, cb, at, method),
             error = function(e) NA_real_)
  }, numeric(1))
  reps <- reps[is.finite(reps)]
  if (length(reps) < 2L) return(NA_real_)
  v <- stats::var(reps)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  1 / v
}


# Score every participant on a shared horizon, dropping any without an estimable
# effect. When `at` is NULL the horizon is the largest time available to *every*
# participant, so the pooled effects are all evaluated at the same point; using
# each participant's own maximum instead would compare cumulative hazards taken
# over different time ranges.
.sft_resilience_score_data <- function(input, method, at = NULL,
                                       at_quantile = 0.5,
                                       var_method = "analytic", n_boot = 2000L) {
  if (length(at_quantile) != 1L || !is.finite(at_quantile) ||
      at_quantile <= 0 || at_quantile >= 1) {
    stop("`at_quantile` must lie strictly between 0 and 1.", call. = FALSE)
  }
  method <- .sft_resilience_method(method)
  var_method <- .sft_var_method(var_method)
  if (is.null(at)) {
    # Never default to the terminal time. For complete (uncensored) data the
    # Nelson-Aalen cumulative hazard at a sample's own maximum is exactly the
    # harmonic number sum_{k=1}^{n} 1/k whatever the distribution, so all three
    # cumulative hazards collapse to the same constant and the contrast becomes
    # -H_n for every dataset -- carrying no information about the UCIP-OR
    # benchmark. Censoring softens this but the tail stays uninformative.
    # A central quantile of the pooled response times is where the cumulative
    # hazards actually separate.
    per_subject_max <- vapply(input$RT, function(rt) {
      max(.sft_finite_times(rt[1:3]))
    }, numeric(1))
    pooled <- unlist(lapply(input$RT, function(rt)
      unlist(rt[1:3], use.names = FALSE)), use.names = FALSE)
    pooled <- pooled[is.finite(pooled)]
    at <- as.numeric(stats::quantile(pooled, at_quantile, names = FALSE,
                                     type = 8))
    # Stay inside every participant's observed range so each contributes a
    # genuine estimate rather than an extrapolated plateau.
    at <- min(at, min(per_subject_max))
  } else {
    at <- as.numeric(at)[[1L]]
    if (!is.finite(at) || at <= 0) {
      stop("`at` must be a positive finite evaluation time.", call. = FALSE)
    }
  }
  scores <- lapply(seq_along(input$RT), function(i) {
    .sft_resilience_score(input$RT[[i]], input$CR[[i]], at = at)
  })

  usable <- vapply(scores, function(s) {
    if (method == "multiplicative") {
      is.finite(s$psi) && is.finite(s$psi_variance) && s$psi_variance > 0
    } else {
      is.finite(s$delta) && is.finite(s$delta_variance) && s$delta_variance > 0
    }
  }, logical(1))
  if (!any(usable)) {
    stop("No participant has an estimable resilience effect with positive ",
         "variance at the evaluation horizon.", call. = FALSE)
  }
  if (!all(usable)) {
    warning("Dropping ", sum(!usable), " participant(s) without an estimable ",
            "resilience effect: ", paste(input$subject[!usable], collapse = ", "),
            ".", call. = FALSE)
    scores <- scores[usable]; input$RT <- input$RT[usable]
    input$CR <- input$CR[usable]; input$subject <- input$subject[usable]
  }

  effect_hat <- vapply(scores, function(s)
    if (method == "multiplicative") s$psi else s$delta, numeric(1))
  variance <- vapply(scores, function(s)
    if (method == "multiplicative") s$psi_variance else s$delta_variance, numeric(1))
  precision <- 1 / variance
  analytic_se <- sqrt(variance)
  used_se <- analytic_se
  if (var_method == "bootstrap") {
    boot_precision <- vapply(seq_along(input$RT), function(i) {
      rt <- lapply(input$RT[[i]][1:3], as.numeric)
      cr <- .sft_cr_list(input$RT[[i]][1:3], input$CR[[i]])
      .sft_resilience_boot_precision(rt, cr, at, method, n_boot)
    }, numeric(1))
    undefined <- !is.finite(boot_precision)
    boot_precision[undefined] <- precision[undefined]
    precision <- boot_precision
    used_se <- 1 / sqrt(precision)
  }

  score_df <- data.frame(
    subject = input$subject,
    time = vapply(scores, `[[`, numeric(1), "time"),
    H_AB = vapply(scores, `[[`, numeric(1), "H_AB"),
    H_A = vapply(scores, `[[`, numeric(1), "H_A"),
    H_B = vapply(scores, `[[`, numeric(1), "H_B"),
    delta = vapply(scores, `[[`, numeric(1), "delta"),
    resilience_ratio = vapply(scores, `[[`, numeric(1), "resilience_ratio"),
    psi = vapply(scores, `[[`, numeric(1), "psi"),
    se = used_se,
    se_analytic = analytic_se,
    z = vapply(scores, `[[`, numeric(1), "z"),
    stringsAsFactors = FALSE
  )
  if (var_method == "analytic") score_df$se_analytic <- NULL
  rownames(score_df) <- NULL

  list(scores = scores, subject = input$subject, at = at,
       effect_hat = effect_hat, precision = precision,
       effect_name = if (method == "multiplicative") "psi" else "delta",
       method = method, var_method = var_method, score = score_df)
}


.sft_resilience_interpretation <- function(method) {
  if (method == "multiplicative") {
    paste("psi = log(H_AB / (H_A + H_B)) is the log resilience ratio;",
          "psi > 0 is super-resilient, psi < 0 is limited, psi = 0 is UCIP-OR")
  } else {
    paste("delta = H_AB - H_A - H_B is the cumulative-hazard contrast against",
          "the UCIP-OR benchmark at the evaluation horizon")
  }
}


# Default priors. psi is a dimensionless log ratio and takes the same weakly
# informative widths as the other log-ratio scales in the package. delta is on an
# absolute cumulative-hazard scale whose magnitude depends on the horizon, so its
# prior is scaled by the observed spread of the participant estimates rather than
# fixed; a fixed width there would be arbitrary.
.sft_resilience_priors <- function(method, effect_hat, prior_sd = NULL,
                                   prior_tau_sd = NULL, prior_shape = NULL,
                                   prior_rate = NULL) {
  if (method == "multiplicative") {
    base_sd <- 0.35; base_tau <- 0.35
  } else {
    spread <- stats::sd(effect_hat)
    if (!is.finite(spread) || spread <= 0) spread <- max(abs(effect_hat), 1)
    base_sd <- 2 * spread; base_tau <- spread
  }
  list(prior_sd = if (is.null(prior_sd)) base_sd else prior_sd,
       prior_tau_sd = if (is.null(prior_tau_sd)) base_tau else prior_tau_sd,
       prior_shape = if (is.null(prior_shape)) 2 else prior_shape,
       prior_rate = if (is.null(prior_rate)) base_tau^2 else prior_rate)
}


#' Analytic Bayesian resilience posterior for one participant.
#'
#' Bayesian companion to \code{\link{resilience.test}}.  The effect is the
#' resilience contrast against the UCIP-OR benchmark evaluated at a single
#' horizon, with the Nelson-Aalen sampling variance carried through to a
#' conjugate normal posterior.
#'
#' @param inData Canonical SFT data frame or a list of AB, A, and B RT vectors.
#' @param CR Optional correctness indicators.
#' @param at Evaluation horizon.  \code{NULL} selects it from \code{at_quantile}.
#' @param at_quantile Pooled response-time quantile defining the horizon when
#'   \code{at} is \code{NULL} (default the median).  The terminal time is
#'   deliberately not the default: for complete data the Nelson-Aalen cumulative
#'   hazard at a sample's own maximum equals the harmonic number
#'   \eqn{\sum_{k=1}^{n} 1/k} whatever the distribution, so all three cumulative
#'   hazards coincide there and the contrast carries no information about the
#'   UCIP-OR benchmark.
#' @param score_method \code{"multiplicative"} (default, the log resilience
#'   ratio \eqn{\psi}) or \code{"difference"} (the cumulative-hazard contrast
#'   \eqn{\delta}).
#' @param var_method,n_boot \code{"analytic"} Nelson-Aalen variance or a
#'   within-condition \code{"bootstrap"} variance.
#' @param ndraws,seed,hdi,prior_mean,prior_sd,chains,rope Posterior and summary
#'   controls, as in \code{\link{ucip.bayes}}.
#' @param Condition,Subject Optional selectors when \code{inData} is a data frame.
#' @return A list mirroring \code{\link{ucip.bayes}}.
#' @seealso \code{\link{resilienceGroup.bayes}}, \code{\link{resilience.test}}
#' @export
resilience.bayes <- function(inData, CR = NULL, at = NULL, at_quantile = 0.5,
                             ndraws = 10000L, seed = NULL, hdi = .94,
                             prior_mean = 0, prior_sd = NULL, chains = 4L,
                             rope = NULL, Condition = NULL, Subject = NULL,
                             score_method = c("multiplicative", "difference"),
                             var_method = c("analytic", "bootstrap"),
                             n_boot = 2000L) {
  score_method <- .sft_resilience_method(score_method)
  var_method <- .sft_var_method(var_method)
  converted <- .sft_as_rt_cr(inData, CR, stopping.rule = "OR",
                             Condition = Condition, Subject = Subject,
                             by_subject = "never")
  .sft_validate_bayes_args(ndraws, 0L, 1L, chains, 1, 1, prior_mean,
                           prior_sd, 1, hdi, rope, .95, 12L)
  .sft_validate_boot(var_method, n_boot)
  input <- .sft_subject_ucip_input(converted$RT, converted$CR)
  if (length(input$RT) != 1L) stop("resilience.bayes() requires exactly one subject.")
  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  scored <- .sft_resilience_score_data(input, score_method, at = at,
                                       at_quantile = at_quantile,
                                       var_method = var_method, n_boot = n_boot)
  priors <- .sft_resilience_priors(score_method, scored$effect_hat,
                                   prior_sd = prior_sd)
  prior_sd <- priors$prior_sd
  effect_name <- scored$effect_name
  fit <- .sft_normal_analytic_fit(scored$effect_hat[[1L]], scored$precision[[1L]],
                                  scored$subject[[1L]], effect_name, ndraws,
                                  chains, prior_mean, prior_sd, hdi, rope)
  bayes_factor <- .sft_analytic_bayes_factor(fit, prior_mean, prior_sd, rope,
                                             point = TRUE, interval = TRUE)
  prior <- list(rope = rope)
  prior[[effect_name]] <- list(mean = prior_mean, sd = prior_sd)
  prior_predictive <- list(score = fit$prior_score)
  prior_predictive[[effect_name]] <- fit$prior_effect
  list(
    statistic = setNames(mean(fit$effect_draws[, 1L]), fit$subject_parameter),
    posterior_probability = list(
      subject_super = setNames(mean(fit$effect_draws[, 1L] > 0), scored$subject),
      subject_limited = setNames(mean(fit$effect_draws[, 1L] < 0), scored$subject),
      subject_rope = setNames(if (is.null(rope)) NA_real_ else
        mean(abs(fit$effect_draws[, 1L]) <= rope), scored$subject)),
    bayes_factor = bayes_factor,
    summary = fit$summary, population_summary = NULL,
    subject_summary = fit$subject_summary, score = scored$score,
    at = scored$at, draws = fit$draws,
    prior_predictive = prior_predictive[c(effect_name, "score")],
    posterior_predictive = fit$posterior_predictive,
    diagnostics = fit$diagnostics,
    shrinkage_summary = fit$shrinkage_summary,
    shrinkage_draws = fit$shrinkage_draws,
    refit_data = list(effect_hat = scored$effect_hat,
                      precision = scored$precision,
                      subject = scored$subject, effect_name = effect_name),
    method = "Analytic single-subject normal-normal posterior",
    method_code = "Analytic",
    alternative = "the participant resilience effect is nonzero",
    prior = prior,
    model = list(
      observation = sprintf("%s_hat ~ Normal(mean = %s, sd = 1 / sqrt(R))",
                            effect_name, effect_name),
      subject = sprintf("%s ~ Normal(mean = prior_mean, sd = prior_sd)", effect_name),
      interpretation = .sft_resilience_interpretation(score_method),
      score_method = score_method, var_method = var_method, hierarchy = FALSE),
    score_method = score_method, var_method = var_method,
    hdi = hdi, ndraws = ndraws, burnin = 0L, thin = 1L, chains = chains,
    seed = seed)
}


#' Hierarchical Bayesian resilience across participants.
#'
#' Pools the resilience contrast against the UCIP-OR benchmark with the same
#' normal-normal machinery as \code{\link{capacityGroup.bayes}}.  A single
#' participant is routed to \code{\link{resilience.bayes}}.
#'
#' All participants are evaluated at one shared horizon, because the effect is a
#' cumulative-hazard quantity and comparing participants observed over different
#' time ranges would not be meaningful.  By default that horizon is the largest
#' time available to every participant.
#'
#' @inheritParams resilience.bayes
#' @param burnin,thin Sampler controls.
#' @param method Hierarchical sampler: \code{"InvGamma"}, \code{"HalfNormal"}, or
#'   \code{"Centered"}.
#' @param prior_shape,prior_rate,prior_tau_sd Between-participant scale priors.
#' @param adapt_delta,max_treedepth,stan_control Stan sampler controls.
#' @param Condition Optional condition selector; supply more than one value to
#'   fit the condition-factor hierarchy and obtain pairwise contrasts.
#' @return A list mirroring \code{\link{capacityGroup.bayes}}.
#' @seealso \code{\link{resilience.bayes}}, \code{\link{capacityGroup.bayes}}
#' @export
resilienceGroup.bayes <- function(inData, CR = NULL, at = NULL,
                                  at_quantile = 0.5, ndraws = 10000L, prior_shape = NULL,
                                  prior_rate = NULL, seed = NULL, hdi = .94,
                                  prior_mean = 0, prior_sd = NULL,
                                  burnin = 1000L, thin = 1L, chains = 4L,
                                  rope = NULL,
                                  method = c("InvGamma", "HalfNormal", "Centered"),
                                  prior_tau_sd = NULL, adapt_delta = .95,
                                  max_treedepth = 12L, stan_control = list(),
                                  Condition = NULL, Subject = NULL,
                                  score_method = c("multiplicative", "difference"),
                                  var_method = c("analytic", "bootstrap"),
                                  n_boot = 2000L) {
  score_method <- .sft_resilience_method(score_method)
  var_method <- .sft_var_method(var_method)
  method <- .sft_bayes_method(method)
  if (!is.list(stan_control)) stop("stan_control must be a list.")
  .sft_validate_bayes_args(ndraws, burnin, thin, chains, prior_shape,
                           prior_rate, prior_mean, prior_sd, prior_tau_sd,
                           hdi, rope, adapt_delta, max_treedepth)
  .sft_validate_boot(var_method, n_boot)

  grouped <- .sft_condition_input(inData, CR, stopping.rule = "OR",
                                  Condition = Condition, Subject = Subject)
  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)

  scored <- .sft_condition_score(
    grouped, function(input) .sft_resilience_score_data(
      input, score_method, at = at, at_quantile = at_quantile,
      var_method = var_method, n_boot = n_boot))

  if (length(scored$subject) == 1L && !scored$multi_condition) {
    return(resilience.bayes(
      inData = grouped$conditions[[1L]]$RT, CR = grouped$conditions[[1L]]$CR,
      at = at, at_quantile = at_quantile, ndraws = ndraws, seed = seed,
      hdi = hdi, prior_mean = prior_mean,
      prior_sd = prior_sd, chains = chains, rope = rope,
      score_method = score_method, var_method = var_method, n_boot = n_boot))
  }
  priors <- .sft_resilience_priors(score_method, scored$effect_hat, prior_sd,
                                   prior_tau_sd, prior_shape, prior_rate)
  .sft_normal_group_result(
    scored = scored, priors = priors, prior_mean = prior_mean, method = method,
    ndraws = ndraws, burnin = burnin, thin = thin, chains = chains, seed = seed,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    stan_control = stan_control, hdi = hdi, rope = rope,
    label = "resilience",
    alternative = "the population resilience effect is nonzero",
    interpretation = .sft_resilience_interpretation(score_method),
    extra = list(at = scored$at, score_method = score_method,
                 var_method = var_method))
}
