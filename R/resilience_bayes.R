# Bayesian companions to the resilience coefficient: resilience.bayes() for one
# participant and resilienceGroup.bayes() for a hierarchy.
#
# Resilience (Little, Eidels, Fific & Wang, 2015) has the *same functional form*
# as the OR capacity coefficient,
#
#   R(t) = H_AB(t) / (H_AY(t) + H_XB(t)),
#
# differing only in which trials enter it: the denominator uses single-target
# displays that still carry distracting information (AY, XB) rather than lone
# targets (A, B). Houpt and Little (2016) accordingly build its statistical test
# on the Houpt-Townsend (2012) weighted-logrank score statistic, exactly as for
# capacity. These functions therefore reuse the UCIP score machinery verbatim,
# with the resilience conditions supplied in place of the capacity conditions,
# and share ucip.bayes()'s three effect scales:
#
#   "score" (default)  theta = U / V, the Peto one-step log hazard ratio for the
#                      resilience contrast, with precision V.
#   "standardized"     phi = theta * sqrt(V_ref), the same effect on the
#                      reference-Cz scale.
#   "multiplicative"   eta = log(A / B), the log weighted-resilience ratio.
#
# All three are exactly zero under the UCIP-OR benchmark. There is no evaluation
# horizon: the score statistic accumulates risk-set-weighted increments over the
# whole observed time course, which is what makes it well behaved where a
# cumulative hazard read off at a single time is not (see resilience.test()).


.sft_resilience_interpretation <- function(score_method) {
  switch(
    score_method,
    multiplicative = "exp(eta) is the participant weighted resilience ratio",
    standardized = paste("phi is the resilience effect on the reference-Cz",
                         "scale; phi = theta * sqrt(V_ref)"),
    paste("theta is the resilience effect on the log-hazard-ratio scale",
          "(the Peto one-step estimate); z_i = theta_i * sqrt(V_i)"))
}


#' Analytic Bayesian resilience posterior for one participant.
#'
#' Bayesian companion to \code{\link{resilience.test}}.  The effect is the
#' Houpt-Townsend score contrast against the UCIP-OR benchmark computed from the
#' resilience conditions, carried through to a conjugate normal posterior.
#'
#' @param inData Canonical SFT data frame, or a list of AB, AY, and XB RT
#'   vectors (double target, and each single target accompanied by a distractor).
#' @param CR Optional correctness indicators.
#' @param score_method \code{"score"} (default, the log-hazard-ratio effect
#'   \eqn{\theta}), \code{"standardized"} (\eqn{\phi}, the same effect on the
#'   reference-Cz scale), or \code{"multiplicative"} (\eqn{\eta = \log(A/B)}).
#' @param var_method,n_boot \code{"analytic"} martingale variance or a
#'   within-condition \code{"bootstrap"} variance.
#' @param ndraws,seed,hdi,prior_mean,prior_sd,chains,rope Posterior and summary
#'   controls, as in \code{\link{ucip.bayes}}.
#' @param Condition,Subject Optional selectors when \code{inData} is a data frame.
#' @return A list mirroring \code{\link{ucip.bayes}}.
#' @seealso \code{\link{resilienceGroup.bayes}}, \code{\link{resilience.test}},
#'   \code{\link{ucip.bayes}}
#' @export
resilience.bayes <- function(inData, CR = NULL,
                             ndraws = 10000L, seed = NULL, hdi = .94,
                             prior_mean = 0, prior_sd = NULL, chains = 4L,
                             rope = NULL, Condition = NULL, Subject = NULL,
                             score_method = c("score", "standardized",
                                              "multiplicative"),
                             var_method = c("analytic", "bootstrap"),
                             n_boot = 2000L) {
  score_method <- .sft_score_method(score_method)
  var_method <- .sft_var_method(var_method)
  converted <- .sft_as_rt_cr(inData, CR, stopping.rule = "OR",
                             Condition = Condition, Subject = Subject,
                             by_subject = "never")
  .sft_validate_bayes_args(ndraws, 0L, 1L, chains, 1, 1, prior_mean,
                           prior_sd, 1, hdi, rope, .95, 12L)
  .sft_validate_boot(var_method, n_boot)
  input <- .sft_subject_ucip_input(converted$RT, converted$CR)
  if (length(input$RT) != 1L) stop("resilience.bayes() requires exactly one subject.")
  if (length(input$RT[[1L]]) < 3L) {
    stop("RT must contain AB, AY, and XB conditions.", call. = FALSE)
  }
  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  scored <- .sft_ucip_score_data(input, "OR", score_method,
                                 var_method = var_method, n_boot = n_boot)
  prior_sd <- .sft_resolve_priors(score_method, scored$V_ref,
                                  prior_sd = prior_sd)$prior_sd
  effect_name <- scored$effect_name
  fit <- .sft_normal_analytic_fit(scored$effect_hat[[1L]], scored$precision[[1L]],
                                  input$subject[[1L]], effect_name, ndraws,
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
      subject_super = setNames(mean(fit$effect_draws[, 1L] > 0), input$subject),
      subject_limited = setNames(mean(fit$effect_draws[, 1L] < 0), input$subject),
      subject_rope = setNames(if (is.null(rope)) NA_real_ else
        mean(abs(fit$effect_draws[, 1L]) <= rope), input$subject)),
    bayes_factor = bayes_factor,
    summary = fit$summary, population_summary = NULL,
    subject_summary = fit$subject_summary, score = scored$score,
    V_ref = scored$V_ref, draws = fit$draws,
    prior_predictive = prior_predictive[c(effect_name, "score")],
    posterior_predictive = fit$posterior_predictive,
    diagnostics = fit$diagnostics,
    shrinkage_summary = fit$shrinkage_summary,
    shrinkage_draws = fit$shrinkage_draws,
    refit_data = list(effect_hat = scored$effect_hat,
                      precision = scored$precision,
                      subject = input$subject, effect_name = effect_name),
    method = "Analytic single-subject normal-normal posterior",
    method_code = "Analytic",
    alternative = "the participant resilience effect is nonzero",
    prior = prior,
    model = list(
      observation = sprintf("%s_hat ~ Normal(mean = %s, sd = 1 / sqrt(%s))",
                            effect_name, effect_name,
                            .sft_precision_name(score_method)),
      subject = sprintf("%s ~ Normal(mean = prior_mean, sd = prior_sd)", effect_name),
      interpretation = .sft_resilience_interpretation(score_method),
      score_method = score_method, var_method = var_method, hierarchy = FALSE),
    score_method = score_method, var_method = var_method,
    hdi = hdi, ndraws = ndraws, burnin = 0L, thin = 1L, chains = chains,
    seed = seed)
}


#' Hierarchical Bayesian resilience across participants.
#'
#' Pools the resilience score contrast against the UCIP-OR benchmark with the
#' same normal-normal machinery as \code{\link{capacityGroup.bayes}}.  A single
#' participant in a single condition is routed to \code{\link{resilience.bayes}}.
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
resilienceGroup.bayes <- function(inData, CR = NULL, ndraws = 10000L,
                                  prior_shape = NULL,
                                  prior_rate = NULL, seed = NULL, hdi = .94,
                                  prior_mean = 0, prior_sd = NULL,
                                  burnin = 1000L, thin = 1L, chains = 4L,
                                  rope = NULL,
                                  method = c("InvGamma", "HalfNormal", "Centered"),
                                  prior_tau_sd = NULL, adapt_delta = .95,
                                  max_treedepth = 12L, stan_control = list(),
                                  Condition = NULL, Subject = NULL,
                                  score_method = c("score", "standardized",
                                                   "multiplicative"),
                                  var_method = c("analytic", "bootstrap"),
                                  n_boot = 2000L) {
  score_method <- .sft_score_method(score_method)
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
    grouped, function(input) .sft_ucip_score_data(
      input, "OR", score_method, var_method = var_method, n_boot = n_boot))

  if (length(scored$effect_hat) == 1L && !scored$multi_condition) {
    one <- grouped$conditions[[1L]]
    return(resilience.bayes(
      inData = one$RT[[1L]], CR = one$CR[[1L]],
      ndraws = ndraws, seed = seed, hdi = hdi, prior_mean = prior_mean,
      prior_sd = prior_sd, chains = chains, rope = rope,
      score_method = score_method, var_method = var_method, n_boot = n_boot))
  }
  if (length(scored$effect_hat) < 2L) {
    stop("Fewer than two estimable participant effects remain; ",
         "use resilience.bayes() for a single participant.", call. = FALSE)
  }

  v_ref <- vapply(scored$per_condition, function(p) p$V_ref, numeric(1))
  if (!scored$multi_condition) v_ref <- unname(v_ref[[1L]])
  priors <- .sft_resolve_priors(score_method, v_ref[[1L]], prior_sd,
                                prior_tau_sd, prior_shape, prior_rate)
  .sft_normal_group_result(
    scored = scored, priors = priors, prior_mean = prior_mean, method = method,
    ndraws = ndraws, burnin = burnin, thin = thin, chains = chains, seed = seed,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    stan_control = stan_control, hdi = hdi, rope = rope,
    label = "resilience score/information",
    alternative = "the population resilience effect is nonzero",
    interpretation = .sft_resilience_interpretation(score_method),
    precision_name = .sft_precision_name(score_method),
    extra = list(V_ref = v_ref, score_method = score_method,
                 var_method = var_method))
}
