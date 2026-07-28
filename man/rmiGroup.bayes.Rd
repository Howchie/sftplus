\name{rmiGroup.bayes}
\alias{rmiGroup.bayes}
\title{Hierarchical Bayesian race model inequality test}
\description{
  Fits the race-model inequality on the probability scale, where it is free of
  participant base time, and reports a population violation curve indexed by
  percentile level together with a single min-over-grid global summary.
}
\usage{
rmiGroup.bayes(inData, CR = NULL, bound = c("miller", "grice"),
               probs = seq(.10, .25, by = .05), ndraws = 10000L,
               prior_shape = NULL, prior_rate = NULL, seed = NULL,
               hdi = .94, prior_mean = 0, prior_sd = NULL,
               burnin = 1000L, thin = 1L, chains = 4L, rope = NULL,
               method = c("InvGamma", "HalfNormal", "Centered"),
               prior_tau_sd = NULL, adapt_delta = .95,
               max_treedepth = 12L, stan_control = list(),
               Condition = NULL, Subject = NULL,
               var_method = c("analytic", "bootstrap"), n_boot = 2000L,
               min_n = 20L, by_block = FALSE, qtype = 5,
               errors = c("discard", "defective"),
               cdf_method = c("polygon", "kernel"),
               kernel_bw = "nrd0")
}
\arguments{
  \item{inData}{Canonical trial data frame with \code{Subject}, \code{RT},
    \code{Correct}, \code{Channel1}, \code{Channel2} and optionally
    \code{Condition}, or a pre-split RT/CR list.}
  \item{CR}{Optional correct-response list when \code{inData} is a list.}
  \item{bound}{\code{"miller"} (the upper bound, i.e. the race model
    inequality) or \code{"grice"} (the lower bound). Unlike
    \code{\link{rmi.test}} the two are fitted separately.}
  \item{probs}{Percentile levels to evaluate. The default is the 10\%--25\%
    range recommended by Kiesel, Miller and Ulrich (2007); see
    \code{\link{rmi.test}}. The multiplicity argument for restricting it is
    weaker here, since the global summary is a single min-over-grid statement
    rather than a family of tests, but the bias argument is not: estimates
    outside that range are the ones most affected by the systematic bias they
    document, and levels above the median are vacuous for the Miller bound.}
  \item{ndraws, burnin, thin, chains}{Sampler settings.}
  \item{seed}{Optional RNG seed.}
  \item{hdi}{Mass of the reported highest-density intervals.}
  \item{prior_mean, prior_sd}{Normal prior on each population mean
    \eqn{\mu_k}. \code{prior_sd} defaults to 0.25: \eqn{\delta} is a difference
    of probabilities, so the default places most prior mass inside
    \eqn{\pm 0.5}, which is already a very large violation. Unlike the UCIP
    score priors, nothing here keys off a data-dependent reference information,
    because the scale of the estimand is known a priori.}
  \item{prior_shape, prior_rate}{Inverse-Gamma prior on \eqn{\tau^2} for
    \code{method = "InvGamma"}.}
  \item{prior_tau_sd}{Half-Normal prior scale on \eqn{\tau} for the Stan
    methods.}
  \item{rope}{Optional half-width of a region of practical equivalence on the
    probability scale.}
  \item{method}{\code{"InvGamma"} (conjugate Gibbs), or \code{"HalfNormal"} /
    \code{"Centered"} (Stan NUTS).}
  \item{adapt_delta, max_treedepth, stan_control}{Stan sampler controls.}
  \item{Condition, Subject}{Optional selectors. With several conditions the
    hierarchy gains one population curve per condition over a shared
    between-participant scale, and every pairwise contrast is reported at each
    percentile level.}
  \item{var_method}{\code{"analytic"} for the closed-form multinomial
    covariance of the cell CDFs, or \code{"bootstrap"} to resample within
    participant, which also propagates uncertainty in the percentile grid
    itself.}
  \item{n_boot}{Bootstrap replicates when \code{var_method = "bootstrap"}.}
  \item{min_n}{Minimum correct trials required per cell; with
    \code{by_block = TRUE} it applies per participant, block, and cell. The
    default of 20 is the sample size Ulrich et al. (2007) and Kiesel et al.
    (2007) recommend; see \code{\link{rmi.test}}.}
  \item{by_block}{Evaluate the inequality separately within each participant
    \emph{and} block, using the \code{Block} column of \code{inData}, and
    average the per-block \eqn{\delta} curves --- the level at which Miller
    (1982) estimated his quantiles. This matters when overall speed drifts
    across the session: pooling a fast late block with a slow early block
    spreads all three cell distributions by the between-block shift, which
    flattens their CDFs and moves the bound, whereas estimating within block
    absorbs the drift block by block. Blocks are disjoint sets of trials and so
    independent, which is why the sampling covariance of the block-averaged
    delta is the sum of the per-block covariances over the squared number of
    blocks. Blocks short of \code{min_n} in any cell are dropped, and a
    participant with no surviving block is dropped as usual. \code{FALSE} (the
    default) treats the whole session as one block and requires no \code{Block}
    column.}
  \item{qtype}{Plotting-position rule for the accompanying time-scale quantile
    tables; see \code{\link{rmi.test}}. It never affects the probability-scale
    estimand.}
  \item{cdf_method, kernel_bw}{Estimator and bandwidth for the accompanying
    time-scale quantile tables; see \code{\link{rmi.test}}. These never affect
    the Bayesian probability-scale estimand.}
  \item{errors}{How incorrect and omitted trials enter the cell CDFs.
    \code{"discard"} keeps correct trials only, as Miller did and as
    \code{\link{estimate.bounds}} does. \code{"defective"} keeps every trial
    and counts an error in the denominator without ever satisfying
    \eqn{RT \le t}, so each cell CDF is the defective distribution
    \eqn{F(t) = P(\mathrm{correct}, RT \le t)} whose asymptote is the accuracy
    rather than 1. Prefer \code{"defective"} when accuracy differs across
    cells: under \code{"discard"} a redundant-target condition that is both
    faster and more accurate has its extra correct responses renormalised away,
    which understates the violation. The classical time-scale route of
    \code{\link{rmi.test}} is correct-only by construction, since a quantile
    above a defective distribution's asymptote does not exist.}
}
\details{
  For participant \eqn{i} and percentile level \eqn{u_k}, the evaluation time
  \eqn{t_{ik}} is the \eqn{u_k}-th quantile of that participant's own pooled
  correct single-target RTs, and the estimand is
  \deqn{\delta_{ik} = F_A(t_{ik}) + F_B(t_{ik}) - F_{AB}(t_{ik})}
  for the Miller bound, or \eqn{F_{AB}(t_{ik}) - \max(F_A(t_{ik}),
  F_B(t_{ik}))} for the Grice bound. Negative values are violations under
  either bound.

  Two properties make this free of participant base time. The inequality is
  evaluated within a participant, so overall speed cancels exactly and no
  millisecond quantity enters the group statistic; and it is invariant to any
  strictly increasing reparameterisation of time applied to all three of a
  participant's conditions, which is what licenses evaluating each participant
  on its own percentile grid and then pooling on the shared \eqn{u} index.
  Averaging millisecond quantiles across participants instead --- Miller's
  original procedure, which \code{\link{rmi.test}} reproduces --- can cancel
  real violations across participants with different base times.

  The hierarchy is
  \deqn{\hat\delta_{ik} \sim N(\delta_{ik}, R_{ik}^{-1/2}), \quad
        \delta_{ik} \sim N(\mu_k, \tau), \quad
        \mu_k \sim N(\mathrm{prior\_mean}, \mathrm{prior\_sd}),}
  the same normal-normal engine the other group models use, with the percentile
  level playing the role the \code{Condition} factor plays there. Two
  consequences are worth stating. \eqn{\tau} is a single between-participant
  scale assumed common across the grid, and participants are independent given
  \eqn{(\mu, \tau)}, so the within-participant covariance across levels enters
  through the marginal variances only. A participant whose informative range
  covers part of the grid contributes to the levels it covers and is absent
  elsewhere, which is why the levels are indexed by percentile rather than by
  time.

  The global summary is \eqn{\min_k \mu_k}. Its posterior probability of being
  negative is the multiplicity-honest counterpart of the classical family of
  per-percentile t-tests, and concerns the tested grid only: a violation
  confined to percentiles outside \code{probs} is not in its scope.

  Percentile levels where the bound is already vacuous are not fitted. For the
  Miller bound \eqn{F_A + F_B} reaches 1 near the median of the single-target
  distributions, above which \eqn{F_{AB} \le 1} satisfies the inequality for
  free, so a fit over the default grid covers roughly the fast half of the
  distribution --- which is where race-model violations occur.
}
\value{
  A list with the standard hierarchical fields (\code{summary},
  \code{population_summary}, \code{subject_summary}, \code{draws},
  \code{diagnostics}, \code{bayes_factor}, \code{prior},
  \code{shrinkage_summary}, and so on), plus:
  \item{percentile}{The population violation curve: one row per condition and
    percentile level with the posterior \code{mean}, \code{median}, HDI, and
    \code{p_violation} \eqn{= P(\mu_k < 0)}.}
  \item{global}{The min-over-grid summary per condition, with
    \code{p_any_violation}, \code{p_bound_holds}, and the percentile level at
    which the violation is worst.}
  \item{percentile_contrasts}{Between-condition differences evaluated
    separately at each percentile level. Also returned as
    \code{condition_contrasts}.}
  \item{quantiles}{Time-scale quantile tables per condition, consumed by
    \code{\link{build_rmi_bound_df}}.}
  \item{statistic}{\code{p_any_violation} per condition.}
  \item{note}{A reminder about the vacuous levels described above.}
  \item{errors}{The error handling actually used.}
}
\references{
  Miller, J. (1982). Divided attention: Evidence for coactivation with
  redundant signals. \emph{Cognitive Psychology}, 14, 247--279.

  Ulrich, R., Miller, J., & Schroeter, H. (2007). Testing the race model
  inequality: An algorithm and computer programs. \emph{Behavior Research
  Methods}, 39, 291--302.

  Kiesel, A., Miller, J., & Ulrich, R. (2007). Systematic biases and Type I
  error accumulation in tests of the race model inequality. \emph{Behavior
  Research Methods}, 39, 539--551.

  Miller, J. (2016). Statistical facilitation and the redundant signals
  effect: What are race and coactivation models? \emph{Attention, Perception,
  & Psychophysics}, 78, 516--519.
}
\seealso{
  \code{\link{rmi.test}}, \code{\link{build_rmi_bound_df}},
  \code{\link{build_rmi_violation_df}}, \code{\link{capacityGroup.bayes}}.
}
\examples{
# Twelve participants with very different base times; AB is coactive.
set.seed(1)
sim <- function(s) {
  base <- runif(1, 200, 900); n <- 150
  data.frame(Subject = paste0("s", s),
             RT = c(base + rexp(n, 1/150)/4, base + rexp(n, 1/150),
                    base + rexp(n, 1/150)),
             Correct = TRUE,
             Channel1 = rep(c(2, 2, 0), each = n),
             Channel2 = rep(c(2, 0, 2), each = n))
}
data <- do.call(rbind, lapply(1:12, sim))

\donttest{
fit <- rmiGroup.bayes(data, probs = seq(.1, .5, by = .1),
                      ndraws = 2000, chains = 2, seed = 1)
fit$global
fit$percentile
}
}
