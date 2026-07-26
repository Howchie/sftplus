# Public Bayesian UCIP orchestrators: ucip.bayes() (single subject) and
# capacityGroup.bayes() (hierarchical group). They score kernels and wrap the
# estimation engine in ucip_bayes_engine.R.


# Analytic normal-normal posterior for one participant. There is no
# population mean or between-participant scale in this model.
ucip.bayes <- function(inData, CR = NULL,
                       stopping.rule = c("OR", "AND", "STST"),
                       ndraws = 10000L, seed = NULL,
                       hdi = .94, prior_mean = 0, prior_sd = NULL,
                       chains = 4L, rope = NULL,
                       Condition = NULL, Subject = NULL,
                       score_method = c("score", "standardized", "multiplicative"),
                       var_method = c("analytic", "bootstrap"), n_boot = 2000L) {
  stopping.rule <- match.arg(stopping.rule)
  score_method <- .sft_score_method(score_method)
  var_method <- .sft_var_method(var_method)
  converted <- .sft_as_rt_cr(inData, CR, stopping.rule = stopping.rule,
                             Condition = Condition, Subject = Subject,
                             by_subject = "never")
  RT <- converted$RT; CR <- converted$CR
  .sft_validate_bayes_args(ndraws, 0L, 1L, chains, 1, 1, prior_mean,
                           prior_sd, 1, hdi, rope, .95, 12L)
  .sft_validate_boot(var_method, n_boot)
  input <- .sft_subject_ucip_input(RT, CR)
  if (length(input$RT) != 1L) {
    stop("ucip.bayes() requires exactly one subject.")
  }
  # The seed is set before scoring so an optional bootstrap variance and the
  # posterior draws share one reproducible RNG stream; analytic scoring is
  # deterministic, so this does not change the analytic posterior.
  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  scored <- .sft_ucip_score_data(input, stopping.rule, score_method,
                                 var_method = var_method, n_boot = n_boot)
  # A NULL prior_sd resolves to the scale-aware default for this score_method
  # (anchored on the theta log-hazard-ratio effect scale and mapped onto
  # whichever scale is active); an explicit value is kept as supplied.
  prior_sd <- .sft_resolve_priors(score_method, scored$V_ref,
                                  prior_sd = prior_sd)$prior_sd
  effect_hat <- scored$effect_hat[[1L]]
  precision <- scored$precision[[1L]]
  effect_name <- scored$effect_name
  observed_name <- paste0(effect_name, "_hat")
  precision_name <- .sft_precision_name(score_method)

  fit <- .sft_normal_analytic_fit(effect_hat, precision, input$subject[[1L]],
                                  effect_name, ndraws, chains, prior_mean,
                                  prior_sd, hdi, rope)
  # An exact point null (effect = 0 is the UCIP prediction) applies on every
  # scale, and each also carries a coherent interval null: theta is a log hazard
  # ratio, phi is a reference-Cz shift, and eta = log(A/B) is a log capacity
  # ratio, so a ROPE is interpretable in all three.
  bayes_factor <- .sft_analytic_bayes_factor(
    fit, prior_mean, prior_sd, rope, point = TRUE, interval = TRUE)
  effect_draws <- fit$effect_draws
  subject_parameter <- fit$subject_parameter
  summary <- fit$summary
  subject_summary <- fit$subject_summary
  posterior_predictive <- fit$posterior_predictive
  prior_predictive <- list(score = fit$prior_score)
  prior_predictive[[effect_name]] <- fit$prior_effect
  prior_predictive <- prior_predictive[c(effect_name, "score")]
  draws <- fit$draws
  prior <- list(rope = rope)
  prior[[effect_name]] <- list(mean = prior_mean, sd = prior_sd)
  list(
    statistic = setNames(mean(effect_draws[, 1L]), subject_parameter),
    posterior_probability = list(
      subject_super = setNames(mean(effect_draws[, 1L] > 0), input$subject[[1L]]),
      subject_limited = setNames(mean(effect_draws[, 1L] < 0), input$subject[[1L]]),
      subject_rope = setNames(if (is.null(rope)) NA_real_ else
                                mean(abs(effect_draws[, 1L]) <= rope),
                              input$subject[[1L]])
    ),
    bayes_factor = bayes_factor,
    summary = summary, population_summary = NULL,
    subject_summary = subject_summary, score = scored$score,
    V_ref = scored$V_ref, draws = draws,
    prior_predictive = prior_predictive,
    posterior_predictive = posterior_predictive,
    diagnostics = fit$diagnostics,
    shrinkage_summary = fit$shrinkage_summary,
    shrinkage_draws = fit$shrinkage_draws,
    # Enough for prior_sensitivity() to refit without rescoring.
    refit_data = list(effect_hat = scored$effect_hat,
                      precision = scored$precision,
                      subject = input$subject, effect_name = effect_name),
    method = "Analytic single-subject normal-normal posterior",
    method_code = "Analytic",
    alternative = paste("the participant UCIP", stopping.rule, "effect is nonzero"),
    prior = prior,
    model = list(
      observation = sprintf("%s ~ Normal(mean = %s, sd = 1 / sqrt(%s))",
                            observed_name, effect_name, precision_name),
      subject = sprintf("%s ~ Normal(mean = prior_mean, sd = prior_sd)",
                        effect_name),
      interpretation = .sft_effect_interpretation(score_method),
      score_method = score_method,
      var_method = var_method,
      hierarchy = FALSE
    ),
    score_method = score_method, var_method = var_method,
    stopping.rule = stopping.rule, hdi = hdi, ndraws = ndraws,
    burnin = 0L, thin = 1L, chains = chains, seed = seed
  )
}


# Hierarchical Bayesian companion to the UCIP/Cz score test.
#
# The default score_method = "score" pools theta = U / V, the Peto one-step
# estimate of the log hazard ratio for the UCIP contrast, with precision V.
# theta is the effect the Cz statistic tests: Cz_i = theta_i * sqrt(V_i)
# reproduces ucip.test() exactly. Its location is stable as trials accumulate
# while V grows linearly with them, which is what makes a fixed prior width a
# genuine effect-size prior that the data can outweigh.
#
# score_method = "standardized" reports the same fit on the reference-Cz scale
# phi = theta * sqrt(V_ref); since that is a deterministic linear map given the
# data, it returns the identical posterior under the default priors, just in Cz
# units. score_method = "multiplicative" instead uses eta = log(A/B) with its
# delta-method precision (formerly "capacity", still accepted as an alias).
#
# The InvGamma method is the exact conjugate Gibbs sampler; HalfNormal is a
# non-centred Stan model and Centered is its centred Stan counterpart, both with
# a half-normal prior on tau.
#
# prior_mean stays at 0: the population score-effect prior is centred on the
# UCIP benchmark (mu = 0 == unlimited-capacity independent parallel). The
# remaining prior widths (prior_sd on mu, prior_tau_sd / prior_shape / prior_rate
# on the between-subject scale tau) default to NULL and are resolved by
# .sft_default_priors() once the data are scored: mu ~ N(0, 0.5^2) and
# tau ~ HalfNormal(0, 0.5) on the theta log-hazard-ratio scale, pushed through
# the sqrt(V_ref) map for "standardized". Supplying any of these arguments
# explicitly overrides its default.
capacityGroup.bayes <- function(inData, CR = NULL, stopping.rule = c("OR", "AND", "STST"),
                            ndraws = 10000L, prior_shape = NULL, prior_rate = NULL,
                            seed = NULL, hdi = .94, prior_mean = 0,
                            prior_sd = NULL, burnin = 1000L, thin = 1L,
                            chains = 4L, rope = NULL,
                            method = c("InvGamma", "HalfNormal", "Centered"),
                            prior_tau_sd = NULL, adapt_delta = .95,
                            max_treedepth = 12L, stan_control = list(),
                            tau2_prior_shape = NULL, tau2_prior_rate = NULL,
                            Condition = NULL, Subject = NULL,
                            score_method = c("score", "standardized", "multiplicative"),
                            var_method = c("analytic", "bootstrap"), n_boot = 2000L) {
  stopping.rule <- match.arg(stopping.rule)
  score_method <- .sft_score_method(score_method)
  var_method <- .sft_var_method(var_method)
  method <- .sft_bayes_method(method)
  if (!is.list(stan_control)) stop("stan_control must be a list.")
  if (!is.null(tau2_prior_shape)) prior_shape <- tau2_prior_shape
  if (!is.null(tau2_prior_rate)) prior_rate <- tau2_prior_rate
  .sft_validate_bayes_args(ndraws, burnin, thin, chains, prior_shape,
                           prior_rate, prior_mean, prior_sd, prior_tau_sd,
                           hdi, rope, adapt_delta, max_treedepth)
  .sft_validate_boot(var_method, n_boot)

  grouped <- .sft_condition_input(inData, CR, stopping.rule = stopping.rule,
                                  Condition = Condition, Subject = Subject)
  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  scored <- .sft_condition_score(
    grouped, function(input) .sft_ucip_score_data(
      input, stopping.rule, score_method, var_method = var_method,
      n_boot = n_boot))

  # A one-subject, one-condition input is deliberately routed to the analytic
  # posterior: the hierarchy cannot identify a population mean and a
  # between-subject scale from one participant.
  if (length(scored$effect_hat) == 1L && !scored$multi_condition) {
    one <- grouped$conditions[[1L]]
    return(ucip.bayes(
      inData = one$RT[[1L]], CR = one$CR[[1L]], stopping.rule = stopping.rule,
      ndraws = ndraws, seed = seed, hdi = hdi, prior_mean = prior_mean,
      prior_sd = prior_sd, chains = chains, rope = rope,
      score_method = score_method, var_method = var_method, n_boot = n_boot))
  }
  if (length(scored$effect_hat) < 2L) {
    stop("Fewer than two estimable participant effects remain; ",
         "use ucip.bayes() for a single participant.", call. = FALSE)
  }

  # V_ref is condition-specific, so report one per condition when several are fit.
  v_ref <- vapply(scored$per_condition, function(p) p$V_ref, numeric(1))
  if (!scored$multi_condition) v_ref <- unname(v_ref[[1L]])
  priors <- .sft_resolve_priors(score_method, v_ref[[1L]], prior_sd,
                                prior_tau_sd, prior_shape, prior_rate)
  .sft_normal_group_result(
    scored = scored, priors = priors, prior_mean = prior_mean, method = method,
    ndraws = ndraws, burnin = burnin, thin = thin, chains = chains, seed = seed,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    stan_control = stan_control, hdi = hdi, rope = rope,
    label = "UCIP score/information",
    alternative = paste("the population UCIP", stopping.rule, "effect is nonzero"),
    interpretation = .sft_effect_interpretation(score_method),
    precision_name = .sft_precision_name(score_method),
    extra = list(V_ref = v_ref, score_method = score_method,
                 var_method = var_method, stopping.rule = stopping.rule))
}
