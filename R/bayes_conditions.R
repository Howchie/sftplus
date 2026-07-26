# Multi-condition support shared by the hierarchical .bayes group models.
#
# When more than one Condition is in play the hierarchy gains one population
# mean per condition over a shared between-subject scale:
#
#   theta_ic ~ Normal(mu_c, tau),   mu_c ~ Normal(prior_mean, prior_sd)
#
# and every pairwise contrast mu_c - mu_c' is reported draw by draw. Subjects
# are treated as independent given (mu, tau), so a participant appearing in
# several conditions is not linked across them: the contrasts are
# between-condition population comparisons, not within-subject differences.
# Fitting each condition separately and comparing the two posteriors would
# instead give each condition its own tau and no shared shrinkage.


# Split canonical input into one RT/CR subject-list per condition.
#
# `Condition = NULL` uses every condition present, which is what makes the
# multi-condition fit the default when the data carry one. Supplying a single
# value reproduces the one-condition fit exactly.
.sft_condition_input <- function(inData, CR = NULL, stopping.rule = "OR",
                                 Condition = NULL, Subject = NULL) {
  if (!is.data.frame(inData)) {
    # List input carries no condition column; treat it as a single condition.
    converted <- .sft_as_rt_cr(inData, CR, stopping.rule = stopping.rule,
                               Condition = Condition, Subject = Subject,
                               by_subject = "always")
    input <- .sft_subject_ucip_input(converted$RT, converted$CR)
    return(list(conditions = list(default = input), levels = "default",
                multi_condition = FALSE))
  }
  data <- .sft_normalize_columns(inData)
  present <- if ("Condition" %in% names(data)) {
    as.character(data$Condition)
  } else rep("default", nrow(data))
  if (anyNA(present)) {
    stop("Condition must not contain missing values.", call. = FALSE)
  }
  levels <- if (is.null(Condition)) {
    sort(unique(present))
  } else {
    chosen <- as.character(Condition)
    unknown <- setdiff(chosen, present)
    if (length(unknown)) {
      stop("Condition value(s) not present in data: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    unique(chosen)
  }
  conditions <- lapply(levels, function(level) {
    converted <- .sft_as_rt_cr(data, NULL, stopping.rule = stopping.rule,
                               Condition = level, Subject = Subject,
                               by_subject = "always")
    .sft_subject_ucip_input(converted$RT, converted$CR)
  })
  names(conditions) <- levels
  list(conditions = conditions, levels = levels,
       multi_condition = length(levels) > 1L)
}


# Apply a per-condition scorer and stack the results, carrying a condition index
# and disambiguating subject labels that recur across conditions.
.sft_condition_score <- function(grouped, scorer) {
  parts <- lapply(grouped$levels, function(level) scorer(grouped$conditions[[level]]))
  names(parts) <- grouped$levels
  keep <- vapply(parts, function(p) length(p$effect_hat) > 0L, logical(1))
  if (!any(keep)) stop("No condition yielded an estimable effect.", call. = FALSE)
  if (!all(keep)) {
    warning("Dropping condition(s) without an estimable effect: ",
            paste(grouped$levels[!keep], collapse = ", "), ".", call. = FALSE)
    parts <- parts[keep]
  }
  levels <- names(parts)
  group <- rep(seq_along(parts), vapply(parts, function(p) length(p$effect_hat), integer(1)))
  subject <- unlist(lapply(parts, `[[`, "subject"), use.names = FALSE)
  label <- if (length(levels) > 1L) {
    paste0(subject, "@", levels[group])
  } else subject
  score <- do.call(rbind, lapply(seq_along(parts), function(k) {
    df <- parts[[k]]$score
    if (length(levels) > 1L) df <- cbind(Condition = levels[[k]], df,
                                         stringsAsFactors = FALSE)
    df
  }))
  rownames(score) <- NULL
  first <- parts[[1L]]
  list(effect_hat = unlist(lapply(parts, `[[`, "effect_hat"), use.names = FALSE),
       precision = unlist(lapply(parts, `[[`, "precision"), use.names = FALSE),
       subject = label, raw_subject = subject,
       group = group, levels = levels,
       multi_condition = length(levels) > 1L,
       effect_name = first$effect_name, score = score,
       at = first$at, per_condition = parts)
}


# Pairwise population contrasts mu_c - mu_c', with sign probabilities, an HDI,
# and (given a rope) an interval-null Bayes factor. The contrast prior is the
# difference of two independent Normal(prior_mean, prior_sd) means, hence
# Normal(0, sqrt(2) * prior_sd).
.sft_condition_contrasts <- function(mu_matrix, levels, hdi, rope, prior_sd) {
  if (ncol(mu_matrix) < 2L) return(NULL)
  pairs <- utils::combn(seq_along(levels), 2L, simplify = FALSE)
  rows <- lapply(pairs, function(idx) {
    d <- mu_matrix[, idx[[1L]]] - mu_matrix[, idx[[2L]]]
    interval <- .sft_hdi(d, hdi)
    bf <- if (is.null(rope)) NULL else
      .sft_interval_null_bf(d, rope, 0, sqrt(2) * prior_sd)
    data.frame(
      contrast = paste(levels[[idx[[1L]]]], "-", levels[[idx[[2L]]]]),
      mean = mean(d), lower = interval[["lower"]], upper = interval[["upper"]],
      posterior_positive = mean(d > 0), posterior_negative = mean(d < 0),
      posterior_rope = if (is.null(rope)) NA_real_ else mean(abs(d) <= rope),
      BF01_interval = if (is.null(bf)) NA_real_ else bf$BF01,
      BF10_interval = if (is.null(bf)) NA_real_ else bf$BF10,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "draws") <- vapply(pairs, function(idx)
    mean(mu_matrix[, idx[[1L]]] - mu_matrix[, idx[[2L]]]), numeric(1))
  out
}


# Assemble the standard hierarchical result list shared by the group models.
.sft_normal_group_result <- function(scored, priors, prior_mean, method, ndraws,
                                     burnin, thin, chains, seed, adapt_delta,
                                     max_treedepth, stan_control, hdi, rope,
                                     label, alternative, interpretation,
                                     point_null = TRUE, interval_null = TRUE,
                                     precision_name = "R", extra = list()) {
  effect_name <- scored$effect_name
  fit <- .sft_normal_hierarchy_fit(
    scored$effect_hat, scored$precision, scored$subject, effect_name, method,
    ndraws, burnin, thin, chains, prior_mean, priors$prior_sd,
    priors$prior_shape, priors$prior_rate, priors$prior_tau_sd, seed,
    adapt_delta, max_treedepth, stan_control, hdi, rope,
    group = scored$group, group_levels = scored$levels)

  mu_matrix <- fit$mu_matrix
  n_group <- ncol(mu_matrix)
  # One null test per condition mean, each using only that condition's units.
  bayes_factor <- if (n_group == 1L) {
    .sft_hierarchy_bayes_factor(fit, prior_mean, priors$prior_sd, rope,
                                point = point_null, interval = interval_null)
  } else {
    setNames(lapply(seq_len(n_group), function(k) {
      out <- list()
      if (point_null) {
        out$point_null <- .sft_point_null_bf(
          prior_mean, priors$prior_sd, effect_draws = fit$effect_draws,
          tau_draws = fit$tau_draws, group = scored$group == k)
      }
      if (interval_null && !is.null(rope)) {
        out$interval_null <- .sft_interval_null_bf(mu_matrix[, k], rope,
                                                   prior_mean, priors$prior_sd)
      }
      if (length(out)) out else NULL
    }), scored$levels)
  }
  contrasts <- .sft_condition_contrasts(mu_matrix, scored$levels, hdi, rope,
                                        priors$prior_sd)

  prior <- if (method == "InvGamma") {
    list(mu = list(mean = prior_mean, sd = priors$prior_sd),
         tau2 = list(family = "inverse-Gamma", shape = priors$prior_shape,
                     rate = priors$prior_rate), rope = rope)
  } else {
    list(mu = list(mean = prior_mean, sd = priors$prior_sd),
         tau = list(family = "half-Normal", sd = priors$prior_tau_sd), rope = rope)
  }
  probs <- list(
    population_super = if (n_group == 1L) mean(mu_matrix[, 1L] > 0) else
      setNames(colMeans(mu_matrix > 0), scored$levels),
    population_limited = if (n_group == 1L) mean(mu_matrix[, 1L] < 0) else
      setNames(colMeans(mu_matrix < 0), scored$levels),
    population_rope = if (is.null(rope)) NA_real_ else if (n_group == 1L)
      mean(abs(mu_matrix[, 1L]) <= rope) else
        setNames(colMeans(abs(mu_matrix) <= rope), scored$levels),
    subject_super = setNames(colMeans(fit$effect_draws > 0), scored$subject),
    subject_limited = setNames(colMeans(fit$effect_draws < 0), scored$subject),
    subject_rope = setNames(if (is.null(rope))
      rep(NA_real_, length(scored$subject)) else
        colMeans(abs(fit$effect_draws) <= rope), scored$subject))

  out <- list(
    statistic = if (n_group == 1L) setNames(mean(mu_matrix[, 1L]), "mu") else
      setNames(colMeans(mu_matrix), fit$mu_names),
    posterior_probability = probs,
    bayes_factor = bayes_factor,
    condition_contrasts = contrasts,
    conditions = scored$levels,
    mu_matrix = mu_matrix,
    summary = fit$summary, population_summary = fit$population,
    subject_summary = fit$subject_summary, score = scored$score,
    draws = fit$draws, prior_predictive = fit$prior_predictive,
    posterior_predictive = fit$posterior_predictive,
    diagnostics = fit$diagnostics, sampler_diagnostics = fit$sampler_diagnostics,
    shrinkage_summary = fit$shrinkage_summary,
    shrinkage_draws = fit$shrinkage_draws, stan_fit = fit$stan_fit,
    # Enough to refit the same hierarchy without rescoring, which is what
    # prior_sensitivity() and spike_slab() reuse.
    refit_data = list(effect_hat = scored$effect_hat,
                      precision = scored$precision, subject = scored$subject,
                      group = scored$group, levels = scored$levels,
                      effect_name = effect_name),
    method = if (method == "InvGamma")
      paste("Hierarchical Bayesian", label, "model") else
        paste(method, "Bayesian", label, "model"),
    method_code = method,
    alternative = alternative,
    prior = prior,
    model = c(list(
      observation = sprintf("%s_hat[i] ~ Normal(mean = %s[i], sd = 1 / sqrt(%s[i]))",
                            effect_name, effect_name, precision_name),
      subject = if (n_group == 1L)
        sprintf("%s[i] ~ Normal(mean = mu, sd = tau)", effect_name) else
          sprintf("%s[i] ~ Normal(mean = mu[condition[i]], sd = tau)", effect_name),
      population = "mu ~ Normal(mean = prior_mean, sd = prior_sd)",
      interpretation = interpretation,
      parameterization = if (method == "Centered" || method == "InvGamma")
        "centered" else "non-centred",
      sampler = if (method == "InvGamma") "conjugate Gibbs" else "Stan NUTS",
      conditions = scored$levels), extra),
    hdi = hdi, ndraws = ndraws, burnin = burnin, thin = thin,
    chains = chains, seed = seed)
  out[names(extra)] <- extra
  out
}
