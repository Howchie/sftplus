# WAIC and PSIS-LOO for the semiparametric hierarchical SFT model.
#
# The likelihood is Poisson over subject-by-series-by-bin cells of the pooled
# time grid, so the pointwise log-likelihood is exactly the per-cell
# poisson_log_lpmf the Stan models now emit. Cells with no exposure contribute
# an identical zero on every draw and are dropped here: they add nothing to
# either lppd or p_waic, but leaving them in would inflate the observation count
# and the standard errors that depend on it.
#
# A caveat that matters for interpretation: these are *cell-wise* criteria, so
# they estimate out-of-sample predictive accuracy for a new bin of an observed
# participant, not for a new participant. For a participant-level comparison the
# folds would have to be whole subjects, which cross-validating this model
# pointwise does not do.


.sft_log_mean_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(mean(exp(x - m)))
}


# WAIC from a draws-by-observation pointwise log-likelihood matrix.
.sft_bayes_waic <- function(log_lik) {
  log_lik <- as.matrix(log_lik)
  keep <- apply(log_lik, 2L, function(col) any(col != 0) && all(is.finite(col)))
  dropped <- sum(!keep)
  log_lik <- log_lik[, keep, drop = FALSE]
  n <- ncol(log_lik)
  if (!n) stop("No usable pointwise log-likelihood columns.", call. = FALSE)
  lppd_i <- apply(log_lik, 2L, .sft_log_mean_exp)
  p_i <- apply(log_lik, 2L, stats::var)
  elpd_i <- lppd_i - p_i
  lppd <- sum(lppd_i); p_waic <- sum(p_i)
  elpd <- lppd - p_waic
  estimates <- data.frame(
    estimate = c(elpd, p_waic, -2 * elpd),
    se = c(sqrt(n) * stats::sd(elpd_i), sqrt(n) * stats::sd(p_i),
           2 * sqrt(n) * stats::sd(elpd_i)),
    row.names = c("elpd_waic", "p_waic", "waic"))
  list(estimates = estimates, waic = -2 * elpd, elpd_waic = elpd,
       p_waic = p_waic, lppd = lppd, n_observations = n,
       dropped_zero_exposure = dropped,
       pointwise = data.frame(elpd_waic = elpd_i, p_waic = p_i,
                              lppd = lppd_i),
       # p_waic per observation above ~0.4 is the usual flag that the pointwise
       # normal approximation behind WAIC is straining.
       n_high_p_waic = sum(p_i > 0.4),
       scale = "deviance (lower is better); elpd_waic is on the log scale (higher is better)")
}


# PSIS-LOO when the optional loo package is available. Its Pareto-k diagnostic
# is the reason to prefer it: WAIC has no comparable warning when the pointwise
# approximation fails.
.sft_bayes_loo <- function(log_lik) {
  if (!requireNamespace("loo", quietly = TRUE)) return(NULL)
  log_lik <- as.matrix(log_lik)
  keep <- apply(log_lik, 2L, function(col) any(col != 0) && all(is.finite(col)))
  tryCatch({
    fit <- loo::loo(log_lik[, keep, drop = FALSE])
    k <- tryCatch(loo::pareto_k_values(fit), error = function(e) NULL)
    list(estimates = fit$estimates, loo = fit,
         n_observations = sum(keep),
         pareto_k_above_0.7 = if (is.null(k)) NA_integer_ else sum(k > 0.7),
         max_pareto_k = if (is.null(k)) NA_real_ else max(k))
  }, error = function(e) NULL)
}


#' Model-comparison criteria for a semiparametric SFT fit.
#'
#' Returns WAIC, and PSIS-LOO when the optional \pkg{loo} package is installed,
#' from the pointwise Poisson log-likelihood of the fitted hazard model. The
#' criteria can compare pooled and salience-split fits or fits with different
#' bin counts, basis dimensions, or smoothness priors.
#'
#' The criteria are cell-wise: they estimate predictive accuracy for a new time
#' bin of an observed participant, not for a new participant. Compare only fits
#' built on the same pooled grid and the same set of observations because the
#' observation set defines the quantity being summed.
#'
#' @param object A fitted \code{sft_bayes} object from
#'   \code{\link{semiparametricSFT.bayes}}.
#' @return A list with \code{waic}, optional \code{loo}, and metadata, or
#'   \code{NULL} when the fit carries no posterior draws.
#' @seealso \code{\link{semiparametricSFT.bayes}}
#' @export
sft_waic <- function(object) {
  if (!inherits(object, "sft_bayes")) {
    stop("sft_waic() needs a fitted sft_bayes object.", call. = FALSE)
  }
  if (is.null(object$waic)) {
    stop("This fit carries no pointwise log-likelihood. Refit with ",
         "semiparametricSFT.bayes(sample = TRUE); the criteria come from the ",
         "Stan generated quantities block and are unavailable for ",
         "posterior_draws input.", call. = FALSE)
  }
  object$waic
}
