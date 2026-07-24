# Public Bayesian UCIP orchestrators: ucip.bayes() (single subject) and
# capacityGroup.bayes() (hierarchical group). They score kernels and wrap the
# estimation engine in ucip_bayes_engine.R.


# Analytic normal-normal posterior for one participant. There is no
# population mean or between-participant scale in this model.
ucip.bayes <- function(inData, CR = NULL,
                       stopping.rule = c("OR", "AND", "STST"),
                       ndraws = 10000L, seed = NULL,
                       hdi = .94, prior_mean = 0, prior_sd = 1,
                       chains = 4L, rope = NULL,
                       Condition = NULL, Subject = NULL,
                       score_method = c("score", "multiplicative"),
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
  effect_hat <- scored$effect_hat[[1L]]
  precision <- scored$precision[[1L]]
  effect_name <- scored$effect_name
  observed_name <- paste0(effect_name, "_hat")
  precision_name <- if (score_method == "multiplicative") "P" else "V"

  fit <- .sft_normal_analytic_fit(effect_hat, precision, input$subject[[1L]],
                                  effect_name, ndraws, chains, prior_mean,
                                  prior_sd, hdi, rope)
  # An exact point null (effect = 0 is the UCIP prediction) applies on either
  # scale; an interval null only suits the dimensionless multiplicative/eta
  # effect, not the information-referenced raw score.
  bayes_factor <- .sft_analytic_bayes_factor(
    fit, prior_mean, prior_sd, rope,
    point = TRUE, interval = identical(score_method, "multiplicative"))
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
    subject_summary = subject_summary, score = scored$score, draws = draws,
    prior_predictive = prior_predictive,
    posterior_predictive = posterior_predictive,
    diagnostics = fit$diagnostics,
    shrinkage_summary = fit$shrinkage_summary,
    shrinkage_draws = fit$shrinkage_draws,
    method = "Analytic single-subject normal-normal posterior",
    method_code = "Analytic",
    alternative = paste("the participant UCIP", stopping.rule, "effect is nonzero"),
    prior = prior,
    model = list(
      observation = sprintf("%s ~ Normal(mean = %s, sd = 1 / sqrt(%s))",
                            observed_name, effect_name, precision_name),
      subject = sprintf("%s ~ Normal(mean = prior_mean, sd = prior_sd)",
                        effect_name),
      interpretation = if (score_method == "multiplicative") {
        "exp(eta) is the participant weighted capacity ratio"
      } else {
        "theta is the unstandardised score-effect parameter"
      },
      score_method = score_method,
      var_method = var_method,
      hierarchy = FALSE
    ),
    score_method = score_method, var_method = var_method,
    stopping.rule = stopping.rule, hdi = hdi, ndraws = ndraws,
    burnin = 0L, thin = 1L, chains = chains, seed = seed
  )
}


# Hierarchical Bayesian companion to the UCIP/Cz score test. The default
# score_method = "score" preserves the original U/V score-effect model;
# score_method = "multiplicative" uses eta = log(A/B) with its delta-method
# precision (formerly "capacity", still accepted as an alias). The InvGamma
# method is the exact conjugate Gibbs sampler.
# HalfNormal is a non-centred Stan model and Centered is its centred Stan
# counterpart; both use a half-normal prior on tau.
#
# prior_mean stays at 0: the population score-effect prior is centred on the
# UCIP benchmark (mu = 0 == unlimited-capacity independent parallel). The
# InvGamma(prior_shape, prior_rate) prior on the between-subject variance tau^2
# defaults to shape 2, rate 0.5 (prior mean tau^2 = 0.5), a weakly informative
# choice that keeps more mass at small tau than the old rate = 1 and better
# matches the half-Normal(0, 1) tau prior used by the Stan variants.
capacityGroup.bayes <- function(inData, CR = NULL, stopping.rule = c("OR", "AND", "STST"),
                            ndraws = 10000L, prior_shape = 2, prior_rate = 0.5,
                            seed = NULL, hdi = .94, prior_mean = 0,
                            prior_sd = 1, burnin = 1000L, thin = 1L,
                            chains = 4L, rope = NULL,
                            method = c("InvGamma", "HalfNormal", "Centered"),
                            prior_tau_sd = 1, adapt_delta = .95,
                            max_treedepth = 12L, stan_control = list(),
                            tau2_prior_shape = NULL, tau2_prior_rate = NULL,
                            Condition = NULL, Subject = NULL,
                            score_method = c("score", "multiplicative"),
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
  converted <- .sft_as_rt_cr(inData, CR, stopping.rule = stopping.rule,
                             Condition = Condition, Subject = Subject,
                             by_subject = "always")
  RT <- converted$RT; CR <- converted$CR
  input <- .sft_subject_ucip_input(RT, CR)

  # A one-subject input is deliberately routed to the analytic posterior.
  # The hierarchy cannot identify a population mean and between-subject scale
  # from one participant, even though priors would make the Gibbs model proper.
  if (length(input$RT) == 1L) {
    return(ucip.bayes(
      inData = RT, CR = CR, stopping.rule = stopping.rule, ndraws = ndraws,
      seed = seed, hdi = hdi, prior_mean = prior_mean, prior_sd = prior_sd,
      chains = chains, rope = rope, score_method = score_method,
      var_method = var_method, n_boot = n_boot
    ))
  }

  # The seed is set before scoring so an optional bootstrap variance and the
  # sampler share one reproducible RNG stream; analytic scoring is deterministic.
  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  scored <- .sft_ucip_score_data(input, stopping.rule, score_method,
                                 var_method = var_method, n_boot = n_boot)
  # Adopt the surviving participants (degenerate subjects may have been dropped).
  input$subject <- scored$subject
  if (length(input$subject) < 2L) {
    stop("Fewer than two participants remain after dropping degenerate subjects; ",
         "use ucip.bayes() for a single participant.",
         call. = FALSE)
  }
  effect_hat <- scored$effect_hat
  precision <- scored$precision
  effect_name <- scored$effect_name
  observed_name <- paste0(effect_name, "_hat")
  precision_name <- if (score_method == "multiplicative") "P" else "V"

  fit <- .sft_normal_hierarchy_fit(effect_hat, precision, input$subject,
                                   effect_name, method, ndraws, burnin, thin,
                                   chains, prior_mean, prior_sd, prior_shape,
                                   prior_rate, prior_tau_sd, seed, adapt_delta,
                                   max_treedepth, stan_control, hdi, rope)
  mu_draws <- fit$mu_draws
  tau_draws <- fit$tau_draws
  effect_draws <- fit$effect_draws
  summary <- fit$summary
  population <- fit$population
  subject_summary <- fit$subject_summary
  posterior_predictive <- fit$posterior_predictive
  shrinkage_draws <- fit$shrinkage_draws
  shrinkage_summary <- fit$shrinkage_summary
  draws <- fit$draws
  diagnostics <- fit$diagnostics
  prior_predictive <- fit$prior_predictive
  stan_fit <- fit$stan_fit
  sampler_diagnostics <- fit$sampler_diagnostics
  # UCIP is an exact theoretical prediction, so the population mean carries a
  # point null (mu = 0); the interval null is reserved for the dimensionless
  # multiplicative/eta effect, whose ROPE is interpretable.
  bayes_factor <- .sft_hierarchy_bayes_factor(
    fit, prior_mean, prior_sd, rope,
    point = TRUE, interval = identical(score_method, "multiplicative"))
  prior <- if (method == "InvGamma") {
    list(mu = list(mean = prior_mean, sd = prior_sd),
         tau2 = list(family = "inverse-Gamma", shape = prior_shape,
                     rate = prior_rate), rope = rope)
  } else {
    list(mu = list(mean = prior_mean, sd = prior_sd),
         tau = list(family = "half-Normal", sd = prior_tau_sd), rope = rope)
  }
  model <- if (method == "InvGamma") {
    list(observation = sprintf("%s[i] ~ Normal(mean = %s[i], sd = 1 / sqrt(%s[i]))",
                               observed_name, effect_name, precision_name),
         subject = sprintf("%s[i] ~ Normal(mean = mu, sd = tau)", effect_name),
         population = "mu ~ Normal(mean = prior_mean, sd = prior_sd)",
         interpretation = if (score_method == "multiplicative") {
           "exp(eta[i]) is the participant weighted capacity ratio"
         } else {
           "theta[i] is the participant unstandardised score-effect parameter"
         },
         score_method = score_method, var_method = var_method,
         parameterization = "centered", sampler = "conjugate Gibbs")
  } else {
    list(observation = sprintf("%s[i] ~ Normal(mean = %s[i], sd = 1 / sqrt(%s[i]))",
                               observed_name, effect_name, precision_name),
         subject = sprintf("%s[i] ~ Normal(mean = mu, sd = tau)", effect_name),
         population = "mu ~ Normal(mean = prior_mean, sd = prior_sd)",
         interpretation = if (score_method == "multiplicative") {
           "exp(eta[i]) is the participant weighted capacity ratio"
         } else {
           "theta[i] is the participant unstandardised score-effect parameter"
         },
         score_method = score_method, var_method = var_method,
         tau = "tau ~ HalfNormal(prior_tau_sd)",
         parameterization = if (method == "Centered") "centered" else "non-centred",
         sampler = "Stan NUTS")
  }
  list(
    statistic = setNames(mean(mu_draws), "mu"),
    posterior_probability = list(
      population_super = mean(mu_draws > 0),
      population_limited = mean(mu_draws < 0),
      population_rope = if (is.null(rope)) NA_real_ else mean(abs(mu_draws) <= rope),
      subject_super = setNames(colMeans(effect_draws > 0), input$subject),
      subject_limited = setNames(colMeans(effect_draws < 0), input$subject),
      subject_rope = setNames(if (is.null(rope)) rep(NA_real_, length(effect_hat)) else
                                colMeans(abs(effect_draws) <= rope), input$subject)
    ),
    bayes_factor = bayes_factor,
    summary = summary, population_summary = population,
    subject_summary = subject_summary, score = scored$score, draws = draws,
    prior_predictive = prior_predictive,
    posterior_predictive = posterior_predictive,
    diagnostics = diagnostics,
    sampler_diagnostics = sampler_diagnostics,
    shrinkage_summary = shrinkage_summary,
    shrinkage_draws = shrinkage_draws,
    stan_fit = stan_fit,
    method = if (method == "InvGamma") {
      "Hierarchical Bayesian UCIP score/information model"
    } else {
      paste(method, "Bayesian UCIP score/information model")
    },
    method_code = method,
    score_method = score_method, var_method = var_method,
    alternative = paste("the population UCIP", stopping.rule, "effect is nonzero"),
    prior = prior, model = model,
    stopping.rule = stopping.rule, hdi = hdi, ndraws = ndraws,
    burnin = burnin, thin = thin, chains = chains, seed = seed
  )
}

