# Semiparametric hierarchical SFT model: public API (semiparametricSFT.bayes(),
# print method, and the sic.bayes()/mic.bayes() extractors).


#' Fit an OR-centred semiparametric hierarchical Bayesian SFT model.
#'
#' Fits a hierarchical piecewise-exponential model of the underlying survival
#' and hazard functions and derives capacity from the posterior; with a
#' `salience_split` it fits the full double-factorial design and additionally
#' derives the survivor interaction contrast (see [sic.bayes()]) and the mean
#' interaction contrast (see [mic.bayes()]) hierarchically at subject and group
#' levels.
#'
#' The input must contain Subject, Condition, RT, Correct, Channel1, and
#' Channel2.  By default positive salience values are pooled by channel
#' presence into AB, A, and B.  With `salience_split`, the model retains
#' A_L/A_H, B_L/B_H, and AB_LL/AB_LH/AB_HL/AB_HH, then hierarchically pools
#' matched capacity deviations across those cells.  Correct finite trials are
#' events; incorrect finite trials are censored exposure.  Fitting requires
#' the optional rstan package, but `sample = FALSE` still performs and returns
#' the complete preparation stage.
#'
#' @param inData Canonical SFT data frame. `data` and `sftData` are accepted
#'   aliases for backward compatibility. Within the data frame, `rt`,
#'   `subjects`, and `LogicalRule` are accepted as aliases for `RT`, `Subject`,
#'   and `Condition`.
#' @param Condition Optional condition value or values to retain.
#' @param n_bins Requested number of pooled RT intervals.
#' @param rt_units Input RT units, either seconds or milliseconds.
#' @param report_quantiles Central pooled RT quantiles retained in tidy curves.
#' @param priors Named prior overrides; see the returned `prior_metadata`.
#' @param smoothness Named smoothness overrides, including `basis_dim`, `A`,
#'   `B`, and `delta`.
#' @param require_complete Require every subject to have AB, A, and B.
#' @param cell_mapping Optional function of Channel1 and Channel2 returning an
#'   OR series label per trial; it may return pooled labels or the salience-split
#'   labels documented below.
#' @param salience_split Retain low/high salience cells.  Use `TRUE` or
#'   `"auto"` when each channel has exactly two positive levels, a numeric
#'   vector of the two levels, a list with `Channel1` and `Channel2` mappings,
#'   or a function returning `L`/`H` labels.  `NULL` preserves pooled mode.
#' @param sample If TRUE, run the optional rstan sampler.
#' @param posterior_draws Optional precomputed log-hazard arrays for testing or
#'   downstream shared-engine use; bypasses rstan when supplied.
#' @param ... Compatibility aliases: `bins`, `nbin`, `grid_bins`, `rt.units`,
#'   `central_quantiles`, `prior`, `smooth`, `fit`, and `salience`.
#' @return An object of class `sft_bayes`.
#' @export
semiparametricSFT.bayes <- function(inData = NULL, Condition = NULL, n_bins = 25L,
                           rt_units = c("seconds", "milliseconds"),
                           report_quantiles = c(0.05, 0.95), priors = NULL,
                           smoothness = NULL, require_complete = TRUE,
                           cell_mapping = NULL, salience_split = NULL,
                           sample = TRUE, posterior_draws = NULL,
                           chains = 4L, iter = 2000L, warmup = NULL,
                           seed = NULL, cores = 1L, refresh = 0L,
                           control = list(adapt_delta = 0.95, max_treedepth = 12L),
                           hdi = 0.94, rope = 0, prior_predictive_draws = 200L,
                           posterior_predictive_draws = 200L, return_fit = TRUE,
                           sftData = NULL, data = NULL, ...) {
  dots <- list(...)
  aliases <- c("bins", "nbin", "grid_bins")
  for (nm in aliases) if (!is.null(dots[[nm]])) n_bins <- dots[[nm]]
  if (!is.null(dots$rt.units)) rt_units <- dots$rt.units
  if (!is.null(dots$central_quantiles)) report_quantiles <- dots$central_quantiles
  if (!is.null(dots$prior)) priors <- dots$prior
  if (!is.null(dots$smooth)) smoothness <- dots$smooth
  if (!is.null(dots$fit)) sample <- dots$fit
  if (!is.null(dots$salience)) salience_split <- dots$salience
  if (!is.null(dots$salience_mapping)) salience_split <- dots$salience_mapping
  unknown <- setdiff(names(dots), c(aliases, "rt.units", "central_quantiles", "prior", "smooth",
                                    "fit", "salience", "salience_mapping"))
  if (length(unknown)) stop("Unused argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  # inData is the canonical argument; data and sftData are accepted aliases.
  supplied <- list(inData = inData, data = data, sftData = sftData)
  supplied <- supplied[!vapply(supplied, is.null, logical(1))]
  if (length(supplied) > 1L) {
    stop("Supply the data once, via inData (data and sftData are accepted aliases).",
         call. = FALSE)
  }
  if (!length(supplied)) stop("inData must be supplied.", call. = FALSE)
  data <- supplied[[1L]]
  rt_units <- .sft_bayes_units(rt_units)
  if (length(hdi) != 1L || !is.finite(hdi) || hdi <= 0 || hdi >= 1)
    stop("hdi must lie strictly between 0 and 1.", call. = FALSE)
  if (length(rope) != 1L || !is.finite(rope) || rope < 0)
    stop("rope must be a non-negative finite value.", call. = FALSE)
  seed <- .sft_capacity_bayes_seed(seed)
  prepared <- .sft_bayes_prepared(data, Condition, n_bins, report_quantiles,
                                  rt_units, priors, smoothness, require_complete,
                                  cell_mapping, salience_split)
  ppc_seed <- seed + 101L
  prior_ppc <- if (prepared$split)
    .sft_bayes_prior_predictive_salience(prepared, prior_predictive_draws, ppc_seed) else
    .sft_bayes_prior_predictive(prepared, prior_predictive_draws, ppc_seed)
  fit <- NULL; raw <- NULL; effects <- NULL; criteria <- NULL
  transformed <- .sft_bayes_empty_transformed()
  diagnostics <- list(available = FALSE, reason = "No posterior fit requested.")
  posterior_ppc <- NULL
  if (!is.null(posterior_draws)) {
    if (!is.list(posterior_draws)) stop("posterior_draws must be a named list of log-hazard arrays.",
                                       call. = FALSE)
    if (prepared$split) {
      if (!all(prepared$series %in% names(posterior_draws)))
        stop("Salience-split posterior_draws must include all eight requested series.", call. = FALSE)
      raw <- setNames(lapply(prepared$series, function(nm) posterior_draws[[nm]]), prepared$series)
      names(raw) <- paste0("eta_", names(raw))
      # Carry population hyper-parameters through when supplied so the typical
      # subject (population parameter) SIC/MIC estimand is available, mirroring
      # the pooled path below.
      pop_names <- c("population_A_L", "population_A_H", "population_B_L",
                     "population_B_H", "population_delta", "population_gamma",
                     "sigma_speed", "sigma_asymmetry", "sigma_capacity")
      for (nm in intersect(pop_names, names(posterior_draws))) raw[[nm]] <- posterior_draws[[nm]]
      effects <- .sft_bayes_effects_from_salience_posterior(raw, prepared$priors,
                                                             prepared$n_subjects, prepared$n_bins,
                                                             prepared$subjects, seed + 17L,
                                                             new_subject = FALSE)
      transformed <- .sft_bayes_transform_salience(effects$log_hazard,
                                                   prepared$grid$bins$width,
                                                   prepared$grid$bins$upper,
                                                   prepared$subjects)
      posterior_ppc <- .sft_bayes_predictive_from_series(prepared, effects$log_hazard,
                                                         seed + 29L, posterior_predictive_draws)
    } else {
      if (!all(c("A", "B", "AB") %in% names(posterior_draws)))
        stop("posterior_draws must be a list with A, B, and AB log-hazard arrays.", call. = FALSE)
      raw <- list(eta_A = posterior_draws$A, eta_B = posterior_draws$B,
                  eta_AB = posterior_draws$AB)
      if (isTRUE(all(c("population_A", "population_B", "population_delta") %in% names(posterior_draws)))) {
        raw$population_A <- posterior_draws$population_A
        raw$population_B <- posterior_draws$population_B
        raw$population_delta <- posterior_draws$population_delta
      }
      effects <- .sft_bayes_effects_from_posterior(raw, prepared$priors,
                                                     prepared$n_subjects, prepared$n_bins,
                                                     prepared$subjects, prepared$grid$bins$width,
                                                     seed + 17L, new_subject = FALSE)
      transformed <- .sft_bayes_transform(effects$log_hazard, prepared$grid$bins$width,
                                           prepared$grid$bins$upper, prepared$subjects)
      posterior_ppc <- .sft_bayes_posterior_predictive(prepared, transformed, seed + 29L,
                                                       posterior_predictive_draws)
    }
  } else if (isTRUE(sample)) {
    sampling <- .sft_bayes_validate_sampling(iter, warmup, chains, cores, refresh, control)
    model <- .sft_bayes_stan_model(prepared$split)
    fit <- rstan::sampling(model, data = prepared$stan_data, iter = sampling$iter,
                           warmup = sampling$warmup, chains = sampling$chains,
                           cores = sampling$cores, seed = seed, refresh = sampling$refresh,
                           control = sampling$control)
    pars <- if (prepared$split) c(
      "beta_A_L", "beta_A_H", "beta_B_L", "beta_B_H", "beta_delta", "beta_gamma_raw",
      "z_speed", "z_asymmetry", "z_capacity", "speed", "asymmetry", "capacity_shift",
      "eta_A_L", "eta_A_H", "eta_B_L", "eta_B_H", "eta_AB_LL", "eta_AB_LH",
      "eta_AB_HL", "eta_AB_HH", "population_A_L", "population_A_H",
      "population_B_L", "population_B_H", "population_delta", "population_gamma",
      "sigma_speed", "sigma_asymmetry", "sigma_capacity", "sigma_smooth_A",
      "sigma_smooth_B", "sigma_smooth_delta", "sigma_gamma", "sigma_smooth_gamma") else c(
        "beta_A", "beta_B", "beta_delta", "z_speed", "z_asymmetry",
        "z_capacity", "speed", "asymmetry", "capacity_shift", "eta_A",
        "eta_B", "eta_AB", "population_A", "population_B",
        "population_delta", "sigma_speed", "sigma_asymmetry", "sigma_capacity",
        "sigma_smooth_A", "sigma_smooth_B", "sigma_smooth_delta")
    raw <- rstan::extract(fit, pars = pars, permuted = TRUE, inc_warmup = FALSE)
    # log_lik is pulled separately: it is large, is not part of the effect
    # machinery, and is dropped again once the criteria are computed.
    log_lik <- tryCatch(
      rstan::extract(fit, pars = "log_lik", permuted = TRUE)$log_lik,
      error = function(e) NULL)
    if (prepared$split) {
      effects <- .sft_bayes_effects_from_salience_posterior(raw, prepared$priors,
                                                             prepared$n_subjects, prepared$n_bins,
                                                             prepared$subjects, seed + 17L,
                                                             new_subject = TRUE)
      transformed <- .sft_bayes_transform_salience(effects$log_hazard,
                                                   prepared$grid$bins$width,
                                                   prepared$grid$bins$upper,
                                                   prepared$subjects, effects$new_log_hazard)
      posterior_ppc <- .sft_bayes_predictive_from_series(prepared, effects$log_hazard,
                                                         seed + 29L, posterior_predictive_draws)
    } else {
      effects <- .sft_bayes_effects_from_posterior(raw, prepared$priors,
                                                     prepared$n_subjects, prepared$n_bins,
                                                     prepared$subjects, prepared$grid$bins$width,
                                                     seed + 17L, new_subject = TRUE)
      transformed <- .sft_bayes_transform(effects$log_hazard, prepared$grid$bins$width,
                                           prepared$grid$bins$upper, prepared$subjects,
                                           effects$new_log_hazard)
      posterior_ppc <- .sft_bayes_posterior_predictive(prepared, transformed, seed + 29L,
                                                       posterior_predictive_draws)
    }
    diagnostics <- .sft_bayes_diagnostics(fit)
    diagnostics$available <- TRUE
    if (!is.null(log_lik)) {
      criteria <- list(waic = .sft_bayes_waic(log_lik),
                       loo = .sft_bayes_loo(log_lik))
      criteria$note <- paste(
        "Cell-wise criteria: predictive accuracy for a new time bin of an",
        "observed participant, not for a new participant. Compare only fits on",
        "the same pooled grid and observation set.")
    }
  }
  if (!is.null(transformed$D)) {
    if (prepared$split) {
      transformed$population <- list(
        D = .sft_bayes_mean_over(transformed$D, 2L, na.rm = FALSE),
        C = exp(.sft_bayes_mean_over(log(transformed$C), 2L)),
        logC = .sft_bayes_mean_over(transformed$logC_average, 2L),
        D_cells = .sft_bayes_mean_over(transformed$D_cells, 2L, na.rm = FALSE),
        C_cells = exp(.sft_bayes_mean_over(log(transformed$C_cells), 2L)),
        SIC = .sft_bayes_mean_over(transformed$SIC, 2L, na.rm = FALSE)
      )
      summaries <- .sft_bayes_split_summaries(transformed, prepared$grid,
                                               prepared$subjects, hdi, rope)
    } else {
      transformed$population <- list(
        D = .sft_bayes_mean_over(transformed$D, 2L, na.rm = FALSE),
        C = .sft_bayes_mean_over(transformed$C, 2L, na.rm = FALSE),
        H_A = .sft_bayes_mean_over(transformed$H_A, 2L, na.rm = FALSE),
        H_B = .sft_bayes_mean_over(transformed$H_B, 2L, na.rm = FALSE),
        H_AB = .sft_bayes_mean_over(transformed$H_AB, 2L, na.rm = FALSE)
      )
      summaries <- .sft_bayes_summaries(transformed, prepared$grid, prepared$subjects, hdi, rope)
    }
  } else {
    empty <- .sft_bayes_empty_curve_data()
    summaries <- list(tidy = empty, subject = empty, population = empty,
                      new_subject = empty)
  }
  # Hierarchical SIC(t) and MIC are only identified from the salience split; the
  # pooled OR fit has no HL/LH/HH/LL cells to contrast.
  mic <- if (prepared$split)
    .sft_bayes_mic_result(transformed, raw, prepared$grid$bins$width,
                          prepared$subjects, hdi, rope) else NULL
  sic <- if (prepared$split)
    .sft_bayes_sic_result(transformed, raw, prepared$grid,
                          prepared$subjects, hdi, rope) else NULL
  posterior <- if (is.null(raw)) NULL else list(
    draws = raw,
    log_hazard = if (!is.null(effects)) effects$log_hazard else NULL,
    n_draws = if (!is.null(transformed$D)) dim(transformed$D)[1L] else 0L
  )
  # The direct log-hazard arrays are more useful and exact than attempting to
  # reconstruct them from cumulative hazards.  Keep them in a compact alias.
  if (!is.null(effects)) posterior$log_hazard <- effects$log_hazard
  if (!is.null(effects)) {
    transformed$log_hazard <- effects$log_hazard
    transformed$hazard <- lapply(effects$log_hazard, exp)
  }
  tidy <- summaries$tidy
  out <- list(
    call = match.call(), method = "OR-centred hierarchical piecewise-exponential",
    stopping_rule = "OR", data = prepared$trials, prepared = prepared,
    grid = prepared$grid, posterior = posterior, fit = if (isTRUE(return_fit)) fit else NULL,
    transformed = transformed, curves = summaries, tidy_curves = tidy,
    curve_data = tidy, sic = sic, mic = mic, diagnostics = diagnostics,
    waic = criteria,
    predictive_checks = list(prior = prior_ppc, posterior = posterior_ppc),
    prior_metadata = list(priors = prepared$priors,
                          units = "log hazards are calibrated in seconds",
                          smoothness = list(type = "second-difference random walk",
                                            basis_dimension = ncol(prepared$basis),
                                            salience_effect = if (prepared$split)
                                              "sum-to-zero gamma cells with strong shrinkage" else NULL),
                          central_interval = prepared$grid$report_range,
                          hdi = hdi, rope = rope),
    salience_pooling = prepared$pooling,
    metadata = list(series = prepared$series, subjects = prepared$subjects,
                    selected_conditions = prepared$selected_conditions,
                    salience_pooling = prepared$pooling,
                    incorrect_trials = "finite censored exposure with no event",
                    complete_subject_series = prepared$complete,
                    missing_subject_series = prepared$missing_series,
                    report_quantiles = report_quantiles,
                    rt_units = prepared$rt_units,
                    posterior_summary = "all capacity summaries are draw-by-draw",
                    sic = if (prepared$split)
                      "posterior SIC = S_AB_LL - S_AB_LH - S_AB_HL + S_AB_HH" else
                      "not identified without an explicit salience split"),
    stan = list(path = tryCatch(.sft_bayes_stan_path(prepared$split), error = function(e) NA_character_),
                data = prepared$stan_data),
    settings = list(n_bins = n_bins, iter = iter, warmup = warmup,
                    chains = chains, seed = seed, sample = isTRUE(sample),
                    return_fit = isTRUE(return_fit))
  )
  class(out) <- c("sft_bayes", "list")
  out
}


print.sft_bayes <- function(x, ...) {
  cat("OR-centred hierarchical Bayesian SFT model (capacity",
      if (!is.null(x$sic)) "+ SIC + MIC" else "", ")\n", sep = "")
  cat("Subjects:", length(x$prepared$subjects), " Series:", paste(x$prepared$series, collapse = ", "),
      " Bins:", x$prepared$n_bins, "\n")
  cat("RT units:", x$prepared$rt_units, "  Conditions:",
      paste(x$prepared$selected_conditions, collapse = ", "), "\n")
  if (isTRUE(x$diagnostics$available)) {
    cat("Posterior draws:", x$posterior$n_draws, "  Divergences:", x$diagnostics$divergences, "\n")
  } else cat("No posterior fit stored; preparation/prior predictive results are available.\n")
  if (!is.null(x$mic)) {
    g <- x$mic$summary[x$mic$summary$Level == "population", , drop = FALSE][1L, ]
    cat(sprintf("Group MIC (%s): mean %.4g  P(MIC>0) = %.2f  P(MIC<0) = %.2f\n",
                x$mic$group_estimand, g$Mean, g$Prob_positive, g$Prob_negative))
  }
  invisible(x)
}


#' Extract the hierarchical posterior MIC from an OR capacity model.
#'
#' The mean interaction contrast MIC = mean(LL) - mean(LH) - mean(HL) + mean(HH)
#' is the signed area under the survivor interaction contrast; it is computed
#' draw by draw at the subject, population (typical-subject parameter),
#' population_finite_mean, and new_subject levels, together with subject-versus-
#' population contrasts, so its sign can be tested against zero at every level.
#' It is only identified from a salience-split fit.
#'
#' @param object A fitted `sft_bayes` object, or a canonical SFT data
#'   frame to fit with [semiparametricSFT.bayes()].
#' @param ... Passed to [semiparametricSFT.bayes()] when `object` is a data frame.  A
#'   `salience_split` is required to identify the cells and defaults to `TRUE`.
#' @return A list with `summary`, `contrasts`, per-level `draws`, and metadata.
#' @export
mic.bayes <- function(object, ...) {
  if (inherits(object, "sft_bayes")) {
    if (is.null(object$mic)) {
      stop("This semiparametricSFT.bayes fit has no MIC; refit with salience_split to ",
           "retain the HH/HL/LH/LL cells the interaction contrast needs.", call. = FALSE)
    }
    return(object$mic)
  }
  if (is.data.frame(object)) {
    dots <- list(...)
    if (is.null(dots$salience_split)) dots$salience_split <- TRUE
    fit <- do.call(semiparametricSFT.bayes, c(list(object), dots))
    if (is.null(fit$mic)) {
      stop("The fit did not identify a MIC; supply a salience_split that yields ",
           "the four HH/HL/LH/LL cells.", call. = FALSE)
    }
    return(fit$mic)
  }
  stop("mic.bayes() needs a fitted sft_bayes object or a canonical SFT data frame.",
       call. = FALSE)
}


#' Extract the hierarchical posterior SIC curves from an OR capacity model.
#'
#' Returns the survivor interaction contrast SIC(t) at the subject (partially
#' pooled), population (typical-subject parameter), population_finite_mean, and
#' new_subject (posterior-predictive) levels, with pointwise posterior sign
#' probabilities.  It is only identified from a salience-split fit.
#'
#' @inheritParams mic.bayes
#' @return A list with a tidy `summary` and per-level views plus metadata.
#' @export
sic.bayes <- function(object, ...) {
  if (inherits(object, "sft_bayes")) {
    if (is.null(object$sic)) {
      stop("This semiparametricSFT.bayes fit has no SIC; refit with salience_split to ",
           "retain the HH/HL/LH/LL cells the interaction contrast needs.", call. = FALSE)
    }
    return(object$sic)
  }
  if (is.data.frame(object)) {
    dots <- list(...)
    if (is.null(dots$salience_split)) dots$salience_split <- TRUE
    fit <- do.call(semiparametricSFT.bayes, c(list(object), dots))
    if (is.null(fit$sic)) {
      stop("The fit did not identify a SIC; supply a salience_split that yields ",
           "the four HH/HL/LH/LL cells.", call. = FALSE)
    }
    return(fit$sic)
  }
  stop("sic.bayes() needs a fitted sft_bayes object or a canonical SFT data frame.",
       call. = FALSE)
}
