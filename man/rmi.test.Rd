\name{rmi.test}
\alias{rmi.test}
\title{Race model inequality tests (Miller 1982; Ulrich, Miller & Schroeter 2007)}
\description{
  Reproduces the classical race-model procedure: per-participant quantiles of
  the bound distribution and of the redundant-target distribution are compared
  at each percentile level with a one-sided t-test across participants. It is
  provided so that a replication can run the classical and the hierarchical
  Bayesian analyses (\code{\link{rmiGroup.bayes}}) on identical data through one
  input path.
}
\usage{
rmi.test(inData, bound = c("miller", "grice", "both"),
         probs = seq(.10, .25, by = .05), p.adjust.method = "holm",
         Condition = NULL, Subject = NULL, min_n = 20L, by_block = FALSE,
         qtype = 5, cdf_method = c("polygon", "kernel"),
         kernel_bw = "nrd0")
}
\arguments{
  \item{inData}{Canonical trial data frame with \code{Subject}, \code{RT},
    \code{Correct}, \code{Channel1}, \code{Channel2} and optionally
    \code{Condition}, or a pre-split RT/CR list supplying the \code{AB},
    \code{A}, and \code{B} cells.}
  \item{bound}{\code{"miller"} (the upper bound \eqn{F_A + F_B}, i.e. the race
    model inequality), \code{"grice"} (the lower bound \eqn{\max(F_A, F_B)}),
    or \code{"both"}.}
  \item{probs}{Percentile levels to test. The default is the 10\%--25\% range
    recommended by Kiesel, Miller and Ulrich (2007): testing across the whole
    distribution accumulates the overall Type I error to roughly 13\%,
    percentiles above the median are vacuous for the Miller bound, and
    10\%--25\% is where violations are in practice found, so restricting the
    range costs little power. Widen it deliberately, not by habit.}
  \item{p.adjust.method}{Multiplicity adjustment passed to
    \code{\link[stats]{p.adjust}}. The classical procedure applies none, which
    is why the default here is \code{"holm"}; the unadjusted \code{p.value} is
    reported alongside. Kiesel et al. (2007) note that Bonferroni-type
    corrections are conservative here, because the tests at adjacent
    percentiles correlate .77 to .95, and prefer restricting \code{probs} as
    the primary control; \code{p.adjust.method = "none"} recovers their
    procedure once the range is restricted.}
  \item{Condition, Subject}{Optional selectors.}
  \item{min_n}{Minimum correct trials required per cell; with
    \code{by_block = TRUE} it applies per participant, block, and cell. The
    default of 20 is the sample size Ulrich et al. (2007) and Kiesel et al.
    (2007) recommend: below it the systematic bias in percentile estimation is
    large enough to manufacture violations, and Kiesel et al. advise treating
    rejections obtained with fewer than 20 trials in any cell with suspicion.
    Participants with a cell below \code{min_n} are dropped with a warning.}
  \item{by_block}{Estimate the quantiles separately within each participant
    \emph{and} block, using the \code{Block} column of \code{inData}, and
    average them over blocks, as Miller (1982) did. It guards against a drift
    in overall speed across the session, which otherwise spreads the
    session-pooled distributions and moves the bound. Blocks short of
    \code{min_n} in any cell are dropped and the number retained is reported as
    \code{n_blocks}. \code{FALSE} (the default) pools the whole session and
    estimates one set of quantiles per participant (and condition).}
  \item{qtype}{Plotting-position rule for the cumulative frequency polygon from
    which both the reported quantiles and the bound are estimated. The default
    5 places the ith of n ordered RTs at \eqn{(i - 0.5)/n}, which is Equation 2
    of Ulrich et al. (2007) and the convention Miller (1982) used, so the
    defaults reproduce the published procedure. Ties are placed at the midpoint
    of the vertical jump they produce, following their Appendix A, whatever the
    rule; this matters because RT is recorded to the nearest millisecond.}
  \item{cdf_method}{\code{"polygon"} (default) for the Ulrich et al.
    cumulative-frequency polygon, or \code{"kernel"} for a fully symmetric
    Gaussian kernel-CDF analysis. In the kernel analysis A, B, and AB, their
    quantiles, and the bound all use the same estimator; \code{qtype} is then
    ignored.}
  \item{kernel_bw}{Bandwidth selector or one positive bandwidth for the kernel
    analysis. The default \code{"nrd0"} applies R's standard rule independently
    in each cell. Other selectors accepted by \code{\link[stats]{density}} may
    also be supplied.}
}
\details{
  The test is on the millisecond scale: at each percentile level the
  redundant-target quantile is compared with the quantile of the bound
  distribution obtained by inverting the bound curve from
  \code{\link{estimate.bounds}}, and the per-participant differences are
  submitted to a one-sided t-test. Negative differences are violations.

  This inherits a known weakness of the classical procedure. The group
  statistic is a mean of millisecond quantities, so a slow participant's 40 ms
  violation outweighs a fast participant's equally severe 15 ms violation, and
  real violations can cancel across participants with different base times.
  \code{\link{rmiGroup.bayes}} tests the same inequality on the probability
  scale, within participant, where overall speed cancels exactly. The
  \code{relative} component returned here divides each difference by the
  participant's grand-mean correct RT and is a cheap partial remedy on the same
  scale as the classical test.

  With \code{cdf_method = "kernel"}, every cell uses
  \eqn{\widehat F_h(t)=n^{-1}\sum_i\Phi((t-x_i)/h)}. The CDF is evaluated
  analytically and all cell and bound quantiles are obtained by numerical
  inversion. This is a symmetric sensitivity analysis, not the hybrid
  kernel-bound/empirical-AB procedure used in the Fox et al. script.
}
\value{
  A list with one element per requested bound (unwrapped when a single bound is
  requested), and within each, one element per condition (unwrapped when a
  single condition is present), containing:
  \item{test}{Per-percentile table with \code{n}, \code{mean_difference},
    \code{t}, \code{df}, \code{p.value}, \code{p.adjusted}, and
    \code{violation}.}
  \item{subject}{Per-participant, per-percentile table with \code{q_A},
    \code{q_B}, \code{q_AB}, \code{q_bound}, \code{difference},
    \code{relative}, \code{mean_RT}, the cell counts \code{n_AB},
    \code{n_A}, \code{n_B}, and \code{n_blocks}, the number of blocks the
    quantiles were averaged over (1 unless \code{by_block = TRUE}).}
  \item{difference, relative}{Group summaries of the raw and the base-time-free
    difference, with t-intervals.}
  \item{quantiles}{Quantile-averaged ("Vincentized") group curves consumed by
    \code{\link{build_rmi_bound_df}}.}
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
}
\seealso{
  \code{\link{rmiGroup.bayes}} for the hierarchical Bayesian analysis on the
  probability scale and \code{\link{build_rmi_bound_df}} for the figure that
  accompanies this test. \code{\link{estimate.bounds}} computes the same
  bounds for general SFT use; this test does not call it, estimating its bound
  from the Ulrich et al. (2007) polygon so that the bound and the
  redundant-target quantile compared against it share one estimator.
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

fit <- rmi.test(data)
head(fit$test)
}
