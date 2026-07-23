\name{capacity.bayes}
\alias{capacity.bayes}
\title{OR-centred hierarchical Bayesian capacity model}
\description{
  Fits a hierarchical piecewise-exponential Bayesian capacity model centred on
  the OR stopping-rule UCIP identity. By default positive salience levels are
  pooled by channel presence into \code{AB}, \code{A}, and \code{B}; an
  explicit \code{salience_split} retains the four matched dual-target cells
  and the corresponding single-target levels. Correct finite response times
  are events and incorrect finite response times contribute censored exposure
  without an event.
}
\usage{
capacity.bayes(data = NULL, Condition = NULL, n_bins = 25L,
  rt_units = c("seconds", "milliseconds"),
  report_quantiles = c(0.05, 0.95), priors = NULL, smoothness = NULL,
  require_complete = TRUE, cell_mapping = NULL, salience_split = NULL,
  sample = TRUE, posterior_draws = NULL,
  chains = 4L, iter = 2000L, warmup = NULL, seed = NULL, cores = 1L,
  refresh = 0L, control = list(adapt_delta = 0.95, max_treedepth = 12L),
  hdi = 0.94, rope = 0, prior_predictive_draws = 200L,
  posterior_predictive_draws = 200L, return_fit = TRUE, sftData = NULL,
  ...)
}
\arguments{
  \item{data}{A data frame with \code{Subject}, \code{Condition}, \code{RT},
    \code{Correct}, \code{Channel1}, and \code{Channel2}.}
  \item{Condition}{Optional condition value or values to retain.}
  \item{n_bins}{Requested number of pooled response-time intervals.}
  \item{rt_units}{Input RT units. Hazards and priors are calibrated in seconds.}
  \item{report_quantiles}{Central pooled RT quantiles retained in tidy curves.}
  \item{priors}{Named prior overrides. See \code{prior_metadata} in the result.}
  \item{smoothness}{Named smoothness overrides, including \code{basis_dim},
    \code{A}, \code{B}, and \code{delta}.}
  \item{require_complete}{Require every subject to have every requested series;
    this is AB/A/B in pooled mode and all eight salience cells in split mode.}
  \item{cell_mapping}{Optional function of \code{Channel1} and \code{Channel2}
    returning pooled labels or salience-split labels per trial.}
  \item{salience_split}{Optional low/high salience split. Use \code{TRUE} or
    \code{"auto"} for exactly two positive levels per channel, a numeric
    two-value mapping, a \code{Channel1}/\code{Channel2} mapping list, or a
    function returning \code{L}/\code{H}. \code{NULL} retains pooled mode.}
  \item{sample}{If \code{TRUE}, fit with the optional \pkg{rstan} package.}
  \item{posterior_draws}{Optional list of precomputed log-hazard arrays, each
    with dimensions draw by subject by bin. The required names are A, B, AB in
    pooled mode and the eight salience series in split mode.}
  \item{chains,iter,warmup,seed,cores,refresh,control}{\pkg{rstan} sampling
    controls.}
  \item{hdi}{Posterior interval mass; defaults to 0.94.}
  \item{rope}{Optional practical-equivalence half-width around 0 for difference
    capacity and 1 for ratio capacity.}
  \item{prior_predictive_draws,posterior_predictive_draws}{Predictive-check
    draw controls.}
  \item{return_fit}{Store the \pkg{rstan} fit in the returned object.}
  \item{sftData}{Alias for \code{data}.}
  \item{...}{Compatibility aliases including \code{bins}, \code{nbin},
    \code{grid_bins}, \code{rt.units}, \code{central_quantiles},
    \code{prior}, \code{smooth}, \code{fit}, and \code{salience}.}
}
\value{
  An object of class \code{sft_bayes_capacity}.  It contains the prepared
  trials, pooled grid and exposure/event arrays, posterior draws and optional
  Stan fit, cumulative hazards, draw-by-draw difference and ratio capacity,
  subject/population/new-subject summaries, predictive checks, diagnostics,
  prior metadata, and exact salience pooling rules.  Tidy curve rows are in
  \code{tidy_curves} and \code{curve_data}; the columns \code{Time}, \code{Ct},
  and \code{Series} make the result compatible with existing plotting helpers.
}
\details{
  The fitted hierarchy is
  \code{log h_A = alpha_A(t) + speed + asymmetry/2},
  \code{log h_B = alpha_B(t) + speed - asymmetry/2}, and
  \code{log h_AB = log_sum_exp(log h_A, log h_B) + delta_population(t) +
  capacity_shift}. In salience-split mode the single-target curves are indexed
  by L/H, the dual-target hazard uses the matched single-target hazards, and a
  sum-to-zero, strongly shrunk \code{gamma} curve permits salience-specific
  departures. Capacity is transformed within each salience cell and averaged
  over cells draw-by-draw. SIC is transformed from the same survivor draws as
  \code{S_AB_LL - S_AB_LH - S_AB_HL + S_AB_HH}. The smooth population
  functions use second-difference random-walk priors and subject effects are
  non-centred.
}
\note{
  The legacy \code{sicDPtest()} and \code{sicGroupBF()} functions are
  unchanged. Split mode supplies posterior SIC curves and uncertainty; model
  classification remains a separate downstream decision. Fitting requires
  the optional dependency \pkg{rstan}; use \code{sample = FALSE} to inspect
  preparation without it.
}
\examples{
\dontrun{
  fit <- capacity.bayes(sftData, Condition = "OR", n_bins = 20, seed = 801)
  head(fit$tidy_curves)
}
}
