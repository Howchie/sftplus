# Hierarchical Bayesian relative mean interaction contrast (MIC).
#
# micGroup.bayes() estimates the mean interaction contrast from the salience
# 2x2 (Channel1, Channel2 in {2 = High, 1 = Low}) cells and pools it across
# participants with the shared normal-normal machinery in ucip_bayes_engine.R.
#
# The estimand is the *relative* MIC,
#
#   rho_i = MIC_i / M_i,
#   MIC_i = mean(HH) - mean(HL) - mean(LH) + mean(LL),
#   M_i   = (mean(HH) + mean(HL) + mean(LH) + mean(LL)) / 4  (equal-weight mean),
#
# a dimensionless proportion of the participant's own RT scale. This removes
# multiplicative speed differences between participants (a generally slower
# participant would otherwise have both a larger raw MIC and a larger raw MIC
# standard error), while preserving the sign, the exact-null boundary
# (rho = 0 iff MIC = 0), and the SIC connection (rho is the normalised signed
# area under the survivor interaction contrast). It also gives a proportional
# ROPE: rope = 0.025 treats interactions below 2.5% of typical RT as negligible.


# The four salience cells (HH, HL, LH, LL) of one or more participants, taken
# from a canonical SFT data frame or supplied directly as RT lists. Cell RTs
# are the finite, correct-by-default response times; 2 = High and 1 = Low on
# each channel, and trials with any other channel level are dropped.
.sft_mic_cells_from_data <- function(data, Condition = NULL, Subject = NULL,
                                     correct = TRUE) {
  data <- .sft_normalize_columns(data)
  required <- c("Subject", "RT", "Correct", "Channel1", "Channel2")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("data is missing required column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!nrow(data)) stop("data must contain at least one trial.", call. = FALSE)
  subject <- as.character(data$Subject)
  if (anyNA(subject) || any(!nzchar(subject))) {
    stop("Subject must not contain missing or empty values.", call. = FALSE)
  }
  condition <- if ("Condition" %in% names(data)) as.character(data$Condition) else
    rep("default", nrow(data))
  if (anyNA(condition)) stop("Condition must not contain missing values.", call. = FALSE)
  rt <- .sft_data_numeric(data$RT, "RT")
  correct_flag <- .sft_data_correct(data$Correct)
  channel1 <- .sft_data_numeric(data$Channel1, "Channel1")
  channel2 <- .sft_data_numeric(data$Channel2, "Channel2")

  if (is.null(Condition)) {
    selected_condition <- unique(condition)
    if (length(selected_condition) != 1L) {
      stop("data contains multiple Condition values; supply Condition to select one.",
           call. = FALSE)
    }
  } else {
    selected_condition <- as.character(Condition)
    if (length(selected_condition) != 1L || is.na(selected_condition) ||
        !selected_condition %in% condition) {
      stop("Condition must select exactly one value present in data.", call. = FALSE)
    }
  }
  keep <- condition == selected_condition & is.finite(rt) &
    channel1 %in% c(1, 2) & channel2 %in% c(1, 2)
  if (isTRUE(correct)) keep <- keep & correct_flag
  if (!any(keep)) {
    stop("No usable High/Low (Channel = 1 or 2) trials remain after selection.",
         call. = FALSE)
  }

  if (is.null(Subject)) {
    selected_subject <- unique(subject[keep])
  } else {
    selected_subject <- as.character(Subject)
    if (!length(selected_subject) || anyNA(selected_subject) ||
        any(!selected_subject %in% subject[keep])) {
      stop("Subject must select values present in the selected condition.", call. = FALSE)
    }
    selected_subject <- unique(selected_subject)
  }

  label <- ifelse(channel1 == 2 & channel2 == 2, "HH",
                  ifelse(channel1 == 2 & channel2 == 1, "HL",
                         ifelse(channel1 == 1 & channel2 == 2, "LH", "LL")))
  cells <- c("HH", "HL", "LH", "LL")
  make_one <- function(id) {
    take <- keep & subject == id
    setNames(lapply(cells, function(cell) rt[take & label == cell]), cells)
  }
  out <- lapply(selected_subject, make_one)
  names(out) <- selected_subject
  list(cells = out, subject = selected_subject, Condition = selected_condition)
}


# Condition levels available to micGroup.bayes(). A data frame with several
# conditions and no explicit selection fits them all; anything else is one level.
.sft_mic_condition_levels <- function(inData, Condition = NULL) {
  if (!is.null(Condition)) return(unique(as.character(Condition)))
  if (!is.data.frame(inData)) return("default")
  data <- .sft_normalize_columns(inData)
  if (!"Condition" %in% names(data)) return("default")
  present <- as.character(data$Condition)
  if (anyNA(present)) stop("Condition must not contain missing values.", call. = FALSE)
  sort(unique(present))
}


# Normalise micGroup.bayes() input to a named list of per-subject cell lists.
# Accepts a canonical data frame, a nested list (one HH/HL/LH/LL list per
# participant), or a flat single-subject HH/HL/LH/LL list.
.sft_mic_input <- function(inData, Condition = NULL, Subject = NULL,
                           correct = TRUE) {
  if (is.data.frame(inData)) {
    ext <- .sft_mic_cells_from_data(inData, Condition = Condition,
                                    Subject = Subject, correct = correct)
    return(list(cells = ext$cells, subject = ext$subject))
  }
  if (!is.list(inData) || !length(inData)) {
    stop("inData must be a data frame or a non-empty list of HH/HL/LH/LL cells.",
         call. = FALSE)
  }
  required <- c("HH", "HL", "LH", "LL")
  flat <- all(vapply(inData, is.numeric, logical(1)))
  if (flat) {
    if (!all(required %in% names(inData))) {
      stop("A single-subject cell list must be named HH, HL, LH, and LL.",
           call. = FALSE)
    }
    return(list(cells = list(`1` = inData[required]), subject = "1"))
  }
  ok <- vapply(inData, function(s) is.list(s) && all(required %in% names(s)),
               logical(1))
  if (!all(ok)) {
    stop("For nested input, each participant entry must be a list named ",
         "HH, HL, LH, and LL.", call. = FALSE)
  }
  subject <- names(inData)
  if (is.null(subject) || any(!nzchar(subject))) subject <- as.character(seq_along(inData))
  cells <- lapply(inData, function(s) s[required])
  names(cells) <- subject
  list(cells = cells, subject = subject)
}


# One participant's relative MIC and its delta-method sampling variance. The
# gradient of rho = (a' m) / (q' m) is (a - rho q) / M, so with independent
# cells and Sigma = diag(s_c^2 / n_c) the variance is grad' Sigma grad.
.sft_mic_relative_score <- function(cell_rt) {
  required <- c("HH", "HL", "LH", "LL")
  x <- lapply(cell_rt[required], function(v) {
    v <- as.numeric(v)
    v[is.finite(v)]
  })
  n <- lengths(x)
  a <- c(HH = 1, HL = -1, LH = -1, LL = 1)
  q <- rep(0.25, 4L)
  if (any(n < 2L)) {
    return(list(mic = NA_real_, mean_rt = NA_real_, rho = NA_real_,
                variance = NA_real_, n = n, ok = FALSE))
  }
  cell_mean <- vapply(x, mean, numeric(1))
  cell_var <- vapply(x, stats::var, numeric(1))
  mic <- sum(a * cell_mean)
  scale <- sum(q * cell_mean)
  rho <- mic / scale
  gradient <- (a - rho * q) / scale
  sigma <- diag(cell_var / n)
  variance <- drop(t(gradient) %*% sigma %*% gradient)
  ok <- all(is.finite(cell_mean)) && all(is.finite(cell_var)) &&
    is.finite(scale) && scale > 0 && is.finite(rho) &&
    is.finite(variance) && variance > 0
  list(mic = mic, mean_rt = scale, rho = rho, variance = variance,
       cell_mean = cell_mean, cell_var = cell_var, n = n, ok = ok)
}


# Within-cell, with-replacement bootstrap variance of the relative MIC, returned
# as a precision. Resampling each cell independently and recomputing rho on every
# replicate captures numerator-denominator dependence and RT skew that the delta
# approximation ignores.
.sft_mic_boot_precision <- function(cell_rt, n_boot) {
  required <- c("HH", "HL", "LH", "LL")
  x <- lapply(cell_rt[required], function(v) {
    v <- as.numeric(v)
    v[is.finite(v)]
  })
  if (any(lengths(x) < 2L)) return(NA_real_)
  reps <- vapply(seq_len(n_boot), function(b) {
    m <- vapply(x, function(v) mean(v[sample.int(length(v), replace = TRUE)]),
                numeric(1))
    (m[[1L]] - m[[2L]] - m[[3L]] + m[[4L]]) / mean(m)
  }, numeric(1))
  reps <- reps[is.finite(reps)]
  if (length(reps) < 2L) return(NA_real_)
  v <- stats::var(reps)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  1 / v
}


# Score every participant, drop those without an estimable relative MIC (with a
# warning), and optionally replace the delta-method precision with a bootstrap
# precision on the same relative scale.
.sft_mic_score_data <- function(cells, subject, var_method = "analytic",
                                n_boot = 2000L) {
  var_method <- .sft_var_method(var_method)
  scores <- lapply(cells, .sft_mic_relative_score)
  ok <- vapply(scores, `[[`, logical(1), "ok")
  if (!any(ok)) {
    stop("No participant has an estimable relative MIC: each of the four ",
         "HH/HL/LH/LL cells needs at least two finite correct RTs and a ",
         "positive grand-mean RT.", call. = FALSE)
  }
  if (!all(ok)) {
    warning("Dropping ", sum(!ok), " participant(s) without an estimable ",
            "relative MIC: ", paste(subject[!ok], collapse = ", "), ".",
            call. = FALSE)
    scores <- scores[ok]; cells <- cells[ok]; subject <- subject[ok]
  }

  rho <- vapply(scores, `[[`, numeric(1), "rho")
  variance <- vapply(scores, `[[`, numeric(1), "variance")
  precision <- 1 / variance
  analytic_se <- sqrt(variance)
  used_se <- analytic_se
  if (var_method == "bootstrap") {
    boot_precision <- vapply(cells, .sft_mic_boot_precision, numeric(1),
                             n_boot = n_boot)
    undefined <- !is.finite(boot_precision)
    boot_precision[undefined] <- precision[undefined]
    precision <- boot_precision
    used_se <- 1 / sqrt(precision)
  }

  cell_n <- do.call(rbind, lapply(scores, function(s) s$n))
  score_df <- data.frame(
    subject = subject,
    MIC = vapply(scores, `[[`, numeric(1), "mic"),
    mean_RT = vapply(scores, `[[`, numeric(1), "mean_rt"),
    relative_MIC = rho,
    se = used_se,
    z = rho * sqrt(precision),
    n_HH = cell_n[, "HH"], n_HL = cell_n[, "HL"],
    n_LH = cell_n[, "LH"], n_LL = cell_n[, "LL"],
    stringsAsFactors = FALSE
  )
  if (var_method == "bootstrap") score_df$se_analytic <- analytic_se
  list(subject = subject, effect_hat = unname(rho), precision = unname(precision),
       effect_name = "rho", var_method = var_method, score = score_df)
}


.sft_mic_subject_probs <- function(effect_draws, subject, rope) {
  list(
    subject_overadditive = setNames(colMeans(effect_draws > 0), subject),
    subject_underadditive = setNames(colMeans(effect_draws < 0), subject),
    subject_additive = setNames(
      if (is.null(rope)) rep(NA_real_, length(subject)) else
        colMeans(abs(effect_draws) <= rope), subject)
  )
}


.sft_mic_model <- function(observed_name, effect_name, precision_name, method,
                           hierarchy, var_method) {
  interpretation <- paste0(
    effect_name, " is the relative MIC (MIC divided by the equal-weight ",
    "grand-mean RT): > 0 is over-additive (parallel self-terminating or ",
    "coactive), < 0 is under-additive (parallel exhaustive), and near 0 is ",
    "additive (serial).")
  if (!hierarchy) {
    return(list(
      observation = sprintf("%s ~ Normal(mean = %s, sd = 1 / sqrt(%s))",
                            observed_name, effect_name, precision_name),
      subject = sprintf("%s ~ Normal(mean = prior_mean, sd = prior_sd)", effect_name),
      interpretation = interpretation, var_method = var_method, hierarchy = FALSE))
  }
  base <- list(
    observation = sprintf("%s[i] ~ Normal(mean = %s[i], sd = 1 / sqrt(%s[i]))",
                          observed_name, effect_name, precision_name),
    subject = sprintf("%s[i] ~ Normal(mean = mu, sd = tau)", effect_name),
    population = "mu ~ Normal(mean = prior_mean, sd = prior_sd)",
    interpretation = interpretation, var_method = var_method)
  if (method == "InvGamma") {
    c(base, list(parameterization = "centered", sampler = "conjugate Gibbs"))
  } else {
    c(base, list(tau = "tau ~ HalfNormal(prior_tau_sd)",
                 parameterization = if (method == "Centered") "centered" else "non-centred",
                 sampler = "Stan NUTS"))
  }
}


#' Hierarchical Bayesian relative mean interaction contrast (MIC).
#'
#' Estimates the salience mean interaction contrast and pools it across
#' participants.  The pooled estimand is the *relative* MIC,
#' \eqn{\rho_i = \mathrm{MIC}_i / M_i}, where \eqn{\mathrm{MIC}_i =
#' \bar T_{HH} - \bar T_{HL} - \bar T_{LH} + \bar T_{LL}} and \eqn{M_i} is the
#' equal-weight grand-mean RT of the four cells.  Dividing by \eqn{M_i} removes
#' multiplicative speed differences between participants while preserving the
#' sign, the exact-null boundary, and the survivor-interaction-contrast
#' interpretation, and it yields a proportional region of practical equivalence.
#'
#' Cells are taken from the salience 2x2 defined by \code{Channel1} and
#' \code{Channel2} (2 = High, 1 = Low).  A single participant is fitted with the
#' analytic normal-normal posterior; two or more are pooled hierarchically with
#' the same conjugate Gibbs (\code{"InvGamma"}) or Stan samplers as
#' \code{\link{capacityGroup.bayes}}.
#'
#' @param inData A canonical SFT data frame (\code{Subject}, \code{RT},
#'   \code{Correct}, \code{Channel1}, \code{Channel2}, optional \code{Condition}),
#'   with `subjects`, `rt`, and `LogicalRule` accepted as aliases for the first,
#'   second, and optional condition columns,
#'   a nested list of per-subject \code{HH}/\code{HL}/\code{LH}/\code{LL} RT
#'   lists, or a flat single-subject \code{HH}/\code{HL}/\code{LH}/\code{LL} list.
#' @param Condition,Subject Optional single condition and subject selectors used
#'   when \code{inData} is a data frame.
#' @param correct Keep only correct trials (the conditional-correct MIC).
#' @param var_method \code{"analytic"} (delta-method, the default) or
#'   \code{"bootstrap"} (within-cell resampling) sampling variance of \eqn{\rho}.
#' @param n_boot Bootstrap replicates when \code{var_method = "bootstrap"}.
#' @param prior_mean,prior_sd Normal prior on the population relative MIC
#'   \eqn{\mu_\rho}.  Because \eqn{\rho} is dimensionless, \code{prior_sd = 0.10}
#'   encodes population MICs within roughly \eqn{\pm 20\%} of mean RT.
#' @param prior_shape,prior_rate Inverse-Gamma prior on \eqn{\tau_\rho^2} for the
#'   \code{"InvGamma"} sampler.
#' @param prior_tau_sd Half-Normal prior scale on \eqn{\tau_\rho} for the Stan
#'   samplers.
#' @param rope Optional proportional region of practical equivalence on
#'   \eqn{\rho}; e.g. \code{0.025} treats interactions under 2.5\% of mean RT as
#'   negligible.
#' @param ndraws,burnin,thin,chains,seed,hdi Sampler and summary controls.
#' @param method One of \code{"InvGamma"}, \code{"HalfNormal"}, or
#'   \code{"Centered"}.
#' @param adapt_delta,max_treedepth,stan_control Stan sampler controls.
#' @param tau2_prior_shape,tau2_prior_rate Back-compatible aliases for
#'   \code{prior_shape} and \code{prior_rate}.
#' @return A list mirroring \code{\link{capacityGroup.bayes}} with an added
#'   \code{estimand} label.  \code{posterior_probability} reports
#'   over-additive / under-additive / additive probabilities for the population
#'   (hierarchical fits) and each subject.
#' @seealso \code{\link{capacityGroup.bayes}}, \code{\link{mic.test}},
#'   \code{\link{mic.bayes}}
#' @export
micGroup.bayes <- function(inData, Condition = NULL, Subject = NULL,
                           correct = TRUE,
                           var_method = c("analytic", "bootstrap"), n_boot = 2000L,
                           ndraws = 10000L, prior_shape = 2, prior_rate = 0.01,
                           seed = NULL, hdi = .94, prior_mean = 0, prior_sd = 0.10,
                           burnin = 1000L, thin = 1L, chains = 4L, rope = NULL,
                           method = c("InvGamma", "HalfNormal", "Centered"),
                           prior_tau_sd = 0.10, adapt_delta = .95,
                           max_treedepth = 12L, stan_control = list(),
                           tau2_prior_shape = NULL, tau2_prior_rate = NULL) {
  var_method <- .sft_var_method(var_method)
  method <- .sft_bayes_method(method)
  if (!is.list(stan_control)) stop("stan_control must be a list.")
  if (!is.null(tau2_prior_shape)) prior_shape <- tau2_prior_shape
  if (!is.null(tau2_prior_rate)) prior_rate <- tau2_prior_rate
  .sft_validate_bayes_args(ndraws, burnin, thin, chains, prior_shape,
                           prior_rate, prior_mean, prior_sd, prior_tau_sd,
                           hdi, rope, adapt_delta, max_treedepth)
  .sft_validate_boot(var_method, n_boot)

  levels <- .sft_mic_condition_levels(inData, Condition)
  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  grouped <- list(
    levels = levels,
    multi_condition = length(levels) > 1L,
    conditions = setNames(lapply(levels, function(level) {
      .sft_mic_input(inData, Condition = if (identical(level, "default")) NULL else level,
                     Subject = Subject, correct = correct)
    }), levels))
  scored <- .sft_condition_score(grouped, function(input)
    .sft_mic_score_data(input$cells, input$subject, var_method, n_boot))

  effect_name <- "rho"
  estimand <- "relative MIC = MIC / equal-weight grand-mean RT"

  if (length(scored$effect_hat) == 1L && !scored$multi_condition) {
    fit <- .sft_normal_analytic_fit(scored$effect_hat[[1L]], scored$precision[[1L]],
                                    scored$subject[[1L]], effect_name, ndraws, chains,
                                    prior_mean, prior_sd, hdi, rope)
    # The relative MIC is dimensionless with a proportional ROPE, so the
    # interval null (negligible interaction) is the meaningful test; an exact
    # point null at rho = 0 is a knife-edge and is not reported.
    bayes_factor <- .sft_analytic_bayes_factor(fit, prior_mean, prior_sd, rope,
                                               point = FALSE, interval = TRUE)
    prior_predictive <- list(score = fit$prior_score)
    prior_predictive[[effect_name]] <- fit$prior_effect
    prior <- list(rope = rope)
    prior[[effect_name]] <- list(mean = prior_mean, sd = prior_sd)
    return(list(
      statistic = setNames(mean(fit$effect_draws[, 1L]), fit$subject_parameter),
      posterior_probability = .sft_mic_subject_probs(fit$effect_draws,
                                                     scored$subject, rope),
      bayes_factor = bayes_factor,
      summary = fit$summary, population_summary = NULL,
      subject_summary = fit$subject_summary, score = scored$score,
      draws = fit$draws,
      prior_predictive = prior_predictive[c(effect_name, "score")],
      posterior_predictive = fit$posterior_predictive,
      diagnostics = fit$diagnostics, shrinkage_summary = fit$shrinkage_summary,
      shrinkage_draws = fit$shrinkage_draws,
      refit_data = list(effect_hat = scored$effect_hat,
                        precision = scored$precision, subject = scored$subject,
                        effect_name = effect_name),
      method = "Analytic single-subject normal-normal posterior",
      method_code = "Analytic",
      alternative = "the participant relative mean interaction contrast is nonzero",
      prior = prior,
      model = .sft_mic_model("rho_hat", effect_name, "V", method = "Analytic",
                             hierarchy = FALSE, var_method = var_method),
      estimand = estimand, var_method = var_method,
      hdi = hdi, ndraws = ndraws, burnin = 0L, thin = 1L, chains = chains,
      seed = seed))
  }

  priors <- list(prior_sd = prior_sd, prior_tau_sd = prior_tau_sd,
                 prior_shape = prior_shape, prior_rate = prior_rate)
  out <- .sft_normal_group_result(
    scored = scored, priors = priors, prior_mean = prior_mean, method = method,
    ndraws = ndraws, burnin = burnin, thin = thin, chains = chains, seed = seed,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    stan_control = stan_control, hdi = hdi, rope = rope,
    label = "relative-MIC",
    alternative = "the population relative mean interaction contrast is nonzero",
    interpretation = .sft_mic_model("rho_hat", effect_name, "V", method = method,
                                    hierarchy = TRUE,
                                    var_method = var_method)$interpretation,
    # rho = 0 is a knife-edge on a dimensionless scale: the interval null is the
    # meaningful architecture test, so the exact point null is not reported.
    point_null = FALSE, interval_null = TRUE, precision_name = "V",
    extra = list(estimand = estimand, var_method = var_method))
  # Report the architecture-facing names alongside the generic super/limited ones.
  out$posterior_probability <- c(
    list(population_overadditive = out$posterior_probability$population_super,
         population_underadditive = out$posterior_probability$population_limited,
         population_additive = out$posterior_probability$population_rope),
    .sft_mic_subject_probs(
      as.matrix(out$draws[, grepl("^rho\\[", names(out$draws)), drop = FALSE]),
      scored$subject, rope))
  out$model <- c(.sft_mic_model("rho_hat", effect_name, "V", method = method,
                                hierarchy = TRUE, var_method = var_method),
                 list(conditions = scored$levels))
  out
}
