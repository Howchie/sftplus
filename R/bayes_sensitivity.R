# Prior sensitivity and spike-and-slab extensions for the normal-normal .bayes
# models. Both refit from the stored per-subject effect estimates and precisions
# (`refit_data`), so neither rescores response times.


.sft_refit_data <- function(fit) {
  rd <- fit$refit_data
  if (is.null(rd) || is.null(rd$effect_hat) || is.null(rd$precision)) {
    stop("This object carries no refit data; prior_sensitivity() and ",
         "spike_slab() need a fit from capacityGroup.bayes(), micGroup.bayes(), ",
         "resilienceGroup.bayes(), ucip.bayes(), or resilience.bayes().",
         call. = FALSE)
  }
  rd
}


#' Prior sensitivity of a Bayesian SFT fit.
#'
#' Refits the model over a grid of alternative-prior widths and reports how the
#' population estimate and the null Bayes factors move.  This is the
#' Jeffreys-Lindley check: a point-null Bayes factor is not prior-free, and a
#' more diffuse alternative mechanically favours the null, so the honest summary
#' of a point null is a curve rather than a number.
#'
#' The interval null is prior-sensitive too, and not necessarily less so:
#' widening the alternative lowers the prior mass inside the region, which
#' raises \code{BF01_interval} directly.  Both columns should be read as curves.
#' What they do give you is two different dependencies -- the point null moves
#' with the prior density at zero, the interval null with the prior mass in the
#' ROPE -- so agreement between them across the grid is meaningful in a way that
#' either number alone is not.
#'
#' @param fit A fitted object from \code{\link{capacityGroup.bayes}},
#'   \code{\link{micGroup.bayes}}, \code{\link{resilienceGroup.bayes}},
#'   \code{\link{ucip.bayes}}, or \code{\link{resilience.bayes}}.
#' @param prior_sd Grid of alternative-prior standard deviations.  \code{NULL}
#'   spans an eightfold range either side of the fitted value.
#' @param ndraws,burnin,chains Sampler controls for each refit.  Defaults are
#'   lighter than a headline fit because only summaries are retained.
#' @param seed Base random seed; each grid point uses a distinct offset.
#' @param rope Practical-equivalence threshold for the interval null.  Defaults
#'   to the one used by the original fit.
#' @return A data frame with one row per grid point, carrying the population
#'   mean and interval, both Bayes factors, and a flag marking the fitted value.
#'   \code{attr(x, "reference_prior_sd")} records that value.
#' @seealso \code{\link{spike_slab}}
#' @export
prior_sensitivity <- function(fit, prior_sd = NULL, ndraws = 4000L,
                              burnin = 500L, chains = 2L, seed = NULL,
                              rope = NULL) {
  rd <- .sft_refit_data(fit)
  reference <- fit$prior$mu$sd
  if (is.null(reference)) {
    reference <- fit$prior[[rd$effect_name]]$sd
  }
  if (is.null(reference) || !is.finite(reference)) {
    stop("Could not read the fitted prior width from this object.", call. = FALSE)
  }
  if (is.null(prior_sd)) {
    prior_sd <- sort(unique(c(reference * c(1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8))))
  }
  prior_sd <- sort(unique(as.numeric(prior_sd)))
  if (!length(prior_sd) || any(!is.finite(prior_sd)) || any(prior_sd <= 0)) {
    stop("prior_sd must be positive finite values.", call. = FALSE)
  }
  if (is.null(rope)) rope <- fit$prior$rope
  prior_mean <- 0
  hdi <- if (is.null(fit$hdi)) .94 else fit$hdi
  hierarchical <- !identical(fit$method_code, "Analytic")
  shape <- fit$prior$tau2$shape; rate <- fit$prior$tau2$rate
  tau_sd <- fit$prior$tau$sd
  method <- if (identical(fit$method_code, "Analytic")) "InvGamma" else fit$method_code

  rows <- lapply(seq_along(prior_sd), function(k) {
    ps <- prior_sd[[k]]
    s <- if (is.null(seed)) NULL else seed + k
    if (hierarchical) {
      # tau's prior is held fixed so the column isolates the effect of the
      # alternative width on mu, which is what the Jeffreys-Lindley concern is
      # about.
      f <- .sft_normal_hierarchy_fit(
        rd$effect_hat, rd$precision, rd$subject, rd$effect_name, method,
        ndraws, burnin, 1L, chains, prior_mean, ps,
        if (is.null(shape)) 2 else shape,
        if (is.null(rate)) 0.25 else rate,
        if (is.null(tau_sd)) 0.5 else tau_sd,
        s, .95, 12L, list(), hdi, rope,
        group = rd$group, group_levels = rd$levels)
      mu <- f$mu_matrix[, 1L]
      point <- .sft_point_null_bf(prior_mean, ps, effect_draws = f$effect_draws,
                                  tau_draws = f$tau_draws,
                                  group = if (is.null(rd$group)) NULL else rd$group == 1L)
    } else {
      f <- .sft_normal_analytic_fit(rd$effect_hat[[1L]], rd$precision[[1L]],
                                    rd$subject[[1L]], rd$effect_name, ndraws,
                                    chains, prior_mean, ps, hdi, rope)
      mu <- f$effect_draws[, 1L]
      point <- .sft_point_null_bf(prior_mean, ps, post_mean = f$post_mean,
                                  post_var = f$post_var)
    }
    interval <- if (is.null(rope)) NULL else
      suppressWarnings(.sft_interval_null_bf(mu, rope, prior_mean, ps))
    ci <- .sft_hdi(mu, hdi)
    data.frame(prior_sd = ps,
               mean = mean(mu), lower = ci[["lower"]], upper = ci[["upper"]],
               posterior_positive = mean(mu > 0),
               BF10_point = point$BF10, BF01_point = point$BF01,
               BF10_interval = if (is.null(interval)) NA_real_ else interval$BF10,
               BF01_interval = if (is.null(interval)) NA_real_ else interval$BF01,
               is_fitted = isTRUE(all.equal(ps, reference)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "reference_prior_sd") <- reference
  attr(out, "rope") <- rope
  attr(out, "note") <- paste(
    "Both Bayes factors fall as prior_sd grows (Jeffreys-Lindley): a diffuse",
    "alternative spreads its mass and mechanically favours the null. The point",
    "null moves with the prior density at zero, the interval null with the",
    "prior mass in the ROPE, so they are not redundant -- but neither is",
    "prior-free. Read the curve, not one row.")
  class(out) <- c("sft_prior_sensitivity", "data.frame")
  out
}


print.sft_prior_sensitivity <- function(x, ...) {
  cat("Prior sensitivity (fitted prior_sd =",
      format(attr(x, "reference_prior_sd"), digits = 3), ")\n")
  if (!is.null(attr(x, "rope"))) cat("ROPE:", attr(x, "rope"), "\n")
  print(as.data.frame(x), row.names = FALSE, digits = 4)
  cat("\n", strwrap(attr(x, "note"), width = 78), sep = "\n")
  invisible(x)
}


# ---- Spike-and-slab: which individual participants depart from the null? ----
#
# The hierarchy's point null answers "is the population mean zero, with tau
# free?". It cannot say which participants individually depart from the
# benchmark, and a large tau can absorb exactly the signal one wants to detect.
# The spike-and-slab model answers the per-participant question directly by
# giving each effect a two-component prior:
#
#   z_i ~ Bernoulli(w)
#   theta_i = 0                         if z_i = 0   (the spike: exactly UCIP)
#   theta_i ~ Normal(mu, tau)           if z_i = 1   (the slab)
#   w ~ Beta(a, b), mu ~ Normal(prior_mean, prior_sd), tau^2 ~ InvGamma(shape, rate)
#
# Every full conditional is available in closed form. Marginalising theta_i
# analytically gives the indicator update
#
#   P(z_i = 1 | .) proportional to w * N(theta_hat_i; mu, tau^2 + 1/precision_i)
#   P(z_i = 0 | .) proportional to (1 - w) * N(theta_hat_i; 0, 1/precision_i)
#
# which mixes far better than sampling theta_i and z_i jointly. mu and tau are
# updated from the slab members only; when none are in the slab they are drawn
# from their priors.
.sft_spike_slab_chain <- function(effect_hat, precision, ndraws, burnin, thin,
                                  prior_mean, prior_sd, prior_shape, prior_rate,
                                  w_shape1, w_shape2) {
  n <- length(effect_hat)
  niter <- burnin + ndraws * thin
  z_draws <- matrix(NA_real_, nrow = ndraws, ncol = n)
  theta_draws <- matrix(NA_real_, nrow = ndraws, ncol = n)
  mu_draws <- numeric(ndraws); tau_draws <- numeric(ndraws); w_draws <- numeric(ndraws)

  tau2 <- 1 / stats::rgamma(1L, shape = prior_shape, rate = prior_rate)
  mu <- stats::rnorm(1L, prior_mean, prior_sd)
  w <- stats::rbeta(1L, w_shape1, w_shape2)
  z <- rep(1L, n)
  theta <- effect_hat
  save_i <- 0L

  for (iter in seq_len(niter)) {
    # Indicators, with theta_i marginalised out.
    log_slab <- log(max(w, .Machine$double.xmin)) +
      stats::dnorm(effect_hat, mu, sqrt(tau2 + 1 / precision), log = TRUE)
    log_spike <- log(max(1 - w, .Machine$double.xmin)) +
      stats::dnorm(effect_hat, 0, sqrt(1 / precision), log = TRUE)
    p_slab <- 1 / (1 + exp(log_spike - log_slab))
    p_slab[!is.finite(p_slab)] <- as.numeric(log_slab[!is.finite(p_slab)] > log_spike[!is.finite(p_slab)])
    z <- stats::rbinom(n, 1L, p_slab)

    # Slab effects from their conjugate normal; spike effects are exactly zero.
    theta <- numeric(n)
    slab <- z == 1L
    if (any(slab)) {
      v <- 1 / (precision[slab] + 1 / tau2)
      m <- v * (precision[slab] * effect_hat[slab] + mu / tau2)
      theta[slab] <- stats::rnorm(sum(slab), m, sqrt(v))
    }

    n_slab <- sum(slab)
    if (n_slab > 0L) {
      mu_var <- 1 / (n_slab / tau2 + 1 / prior_sd^2)
      mu_mean <- mu_var * (sum(theta[slab]) / tau2 + prior_mean / prior_sd^2)
      mu <- stats::rnorm(1L, mu_mean, sqrt(mu_var))
      tau2 <- 1 / stats::rgamma(
        1L, shape = prior_shape + n_slab / 2,
        rate = prior_rate + sum((theta[slab] - mu)^2) / 2)
    } else {
      mu <- stats::rnorm(1L, prior_mean, prior_sd)
      tau2 <- 1 / stats::rgamma(1L, shape = prior_shape, rate = prior_rate)
    }
    w <- stats::rbeta(1L, w_shape1 + n_slab, w_shape2 + n - n_slab)

    if (iter > burnin && (iter - burnin) %% thin == 0L) {
      save_i <- save_i + 1L
      z_draws[save_i, ] <- z
      theta_draws[save_i, ] <- theta
      mu_draws[save_i] <- mu; tau_draws[save_i] <- sqrt(tau2); w_draws[save_i] <- w
      if (save_i == ndraws) break
    }
  }
  list(z = z_draws, theta = theta_draws, mu = mu_draws, tau = tau_draws, w = w_draws)
}


#' Per-participant spike-and-slab test of the SFT null.
#'
#' Asks which individual participants depart from the benchmark, rather than
#' whether the population mean does.  Each participant's effect is either
#' exactly zero (the spike) or drawn from the population slab, and the posterior
#' inclusion probability \eqn{P(z_i = 1)} is that participant's evidence for a
#' genuine effect.
#'
#' Each participant also gets a Bayes factor from the change in inclusion odds
#' relative to the prior inclusion probability, which is the prior mean of
#' \eqn{w}.  Because \eqn{w} is estimated, the participants shrink toward one
#' another: an isolated departure in an otherwise null sample is discounted,
#' which is the multiplicity control an independent per-participant test lacks.
#'
#' @param fit A fitted object carrying \code{refit_data}; see
#'   \code{\link{prior_sensitivity}}.
#' @param w_shape1,w_shape2 Beta prior on the inclusion rate \eqn{w}.  The
#'   default \code{Beta(1, 1)} is uniform, giving a prior inclusion probability
#'   of one half.
#' @param ndraws,burnin,thin,chains,seed Sampler controls.
#' @param hdi Interval mass for the slab-effect summaries.
#' @return A list with a per-participant \code{summary} (inclusion probability,
#'   Bayes factor, and posterior effect), a \code{population} summary of
#'   \eqn{w}, \eqn{\mu}, and \eqn{\tau}, the \code{draws}, and MCMC
#'   \code{diagnostics}.
#' @seealso \code{\link{prior_sensitivity}}, \code{\link{capacityGroup.bayes}}
#' @export
spike_slab <- function(fit, w_shape1 = 1, w_shape2 = 1, ndraws = 10000L,
                       burnin = 1000L, thin = 1L, chains = 4L, seed = NULL,
                       hdi = NULL) {
  rd <- .sft_refit_data(fit)
  n <- length(rd$effect_hat)
  if (n < 2L) {
    stop("spike_slab() needs at least two participants; with one there is ",
         "nothing to share an inclusion rate with.", call. = FALSE)
  }
  if (chains < 2L || chains != as.integer(chains)) stop("chains must be an integer >= 2.")
  if (any(c(w_shape1, w_shape2) <= 0) || any(!is.finite(c(w_shape1, w_shape2)))) {
    stop("w_shape1 and w_shape2 must be positive finite values.", call. = FALSE)
  }
  if (is.null(hdi)) hdi <- if (is.null(fit$hdi)) .94 else fit$hdi
  prior_mean <- 0
  prior_sd <- fit$prior$mu$sd
  if (is.null(prior_sd)) prior_sd <- fit$prior[[rd$effect_name]]$sd
  shape <- fit$prior$tau2$shape; rate <- fit$prior$tau2$rate
  if (is.null(shape)) shape <- 2
  if (is.null(rate)) rate <- (if (is.null(fit$prior$tau$sd)) 0.5 else fit$prior$tau$sd)^2

  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  chain_draws <- lapply(seq_len(chains), function(k)
    .sft_spike_slab_chain(rd$effect_hat, rd$precision, ndraws, burnin, thin,
                          prior_mean, prior_sd, shape, rate, w_shape1, w_shape2))

  z <- do.call(rbind, lapply(chain_draws, `[[`, "z"))
  theta <- do.call(rbind, lapply(chain_draws, `[[`, "theta"))
  mu <- unlist(lapply(chain_draws, `[[`, "mu"), use.names = FALSE)
  tau <- unlist(lapply(chain_draws, `[[`, "tau"), use.names = FALSE)
  w <- unlist(lapply(chain_draws, `[[`, "w"), use.names = FALSE)

  # Prior inclusion odds come from the Beta prior mean, so the per-participant
  # Bayes factor is the posterior-to-prior shift in inclusion odds.
  prior_inclusion <- w_shape1 / (w_shape1 + w_shape2)
  prior_odds <- prior_inclusion / (1 - prior_inclusion)
  pip <- colMeans(z)
  total <- nrow(z)
  # Jeffreys smoothing keeps a saturated participant's Bayes factor finite, as
  # in the interval-null helper.
  smoothed <- (colSums(z) + 0.5) / (total + 1)
  bf10 <- (smoothed / (1 - smoothed)) / prior_odds

  effect_summary <- t(vapply(seq_len(n), function(i) {
    slab <- z[, i] == 1L
    v <- if (any(slab)) theta[slab, i] else NA_real_
    ci <- if (any(slab)) .sft_hdi(v, hdi) else c(lower = NA_real_, upper = NA_real_)
    c(mean(theta[, i]), if (any(slab)) mean(v) else NA_real_,
      ci[["lower"]], ci[["upper"]])
  }, numeric(4)))

  summary <- data.frame(
    subject = rd$subject,
    effect_hat = rd$effect_hat,
    se = 1 / sqrt(rd$precision),
    inclusion_probability = pip,
    BF10_inclusion = bf10,
    BF01_inclusion = 1 / bf10,
    mean_effect = effect_summary[, 1L],
    slab_mean = effect_summary[, 2L],
    slab_lower = effect_summary[, 3L],
    slab_upper = effect_summary[, 4L],
    mc_resolution_limited = pip == 0 | pip == 1,
    stringsAsFactors = FALSE)
  if (!is.null(rd$group) && !is.null(rd$levels) && length(rd$levels) > 1L) {
    summary <- cbind(Condition = rd$levels[rd$group], summary,
                     stringsAsFactors = FALSE)
  }
  rownames(summary) <- NULL

  pop <- do.call(rbind, lapply(list(w = w, mu = mu, tau = tau), function(v) {
    ci <- .sft_hdi(v, hdi)
    data.frame(mean = mean(v), lower = ci[["lower"]], upper = ci[["upper"]])
  }))
  pop <- cbind(parameter = c("w", "mu", "tau"), pop, stringsAsFactors = FALSE)
  rownames(pop) <- NULL

  chain_array <- array(NA_real_, dim = c(ndraws, chains, 3L))
  for (k in seq_len(chains)) {
    chain_array[, k, 1L] <- chain_draws[[k]]$w
    chain_array[, k, 2L] <- chain_draws[[k]]$mu
    chain_array[, k, 3L] <- chain_draws[[k]]$tau
  }
  list(
    summary = summary,
    population = pop,
    expected_nonnull = mean(rowSums(z)),
    prior_inclusion_probability = prior_inclusion,
    draws = list(z = z, theta = theta, mu = mu, tau = tau, w = w),
    diagnostics = .sft_mcmc_diagnostics(chain_array, c("w", "mu", "tau")),
    prior = list(mu = list(mean = prior_mean, sd = prior_sd),
                 tau2 = list(shape = shape, rate = rate),
                 w = list(family = "Beta", shape1 = w_shape1, shape2 = w_shape2)),
    effect_name = rd$effect_name,
    model = paste("Spike-and-slab:", rd$effect_name,
                  "[i] = 0 with probability 1 - w, else Normal(mu, tau)"),
    interpretation = paste(
      "inclusion_probability is P(participant departs from the null).",
      "BF10_inclusion is the posterior-to-prior shift in inclusion odds.",
      "Because w is estimated, participants shrink toward one another, so an",
      "isolated departure in an otherwise null sample is discounted."),
    hdi = hdi, ndraws = ndraws, burnin = burnin, thin = thin, chains = chains,
    seed = seed)
}
