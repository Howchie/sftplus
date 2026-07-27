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
         probs = seq(.05, .95, by = .05), p.adjust.method = "holm",
         Condition = NULL, Subject = NULL, min_n = 5L, qtype = 8)
}
\arguments{
  \item{inData}{Canonical trial data frame with \code{Subject}, \code{RT},
    \code{Correct}, \code{Channel1}, \code{Channel2} and optionally
    \code{Condition}, or a pre-split RT/CR list supplying the \code{AB},
    \code{A}, and \code{B} cells.}
  \item{bound}{\code{"miller"} (the upper bound \eqn{F_A + F_B}, i.e. the race
    model inequality), \code{"grice"} (the lower bound \eqn{\max(F_A, F_B)}),
    or \code{"both"}.}
  \item{probs}{Percentile levels to test.}
  \item{p.adjust.method}{Multiplicity adjustment passed to
    \code{\link[stats]{p.adjust}}. The classical procedure applies none, which
    is why the default here is \code{"holm"}; the unadjusted \code{p.value} is
    reported alongside.}
  \item{Condition, Subject}{Optional selectors.}
  \item{min_n}{Minimum correct trials required per cell.}
  \item{qtype}{Plotting-position rule for the reported quantiles, passed to
    \code{\link[stats]{quantile}}. The default 8 matches the rest of the
    package; \code{qtype = 5} places the ith of n ordered RTs at
    \eqn{(i - 0.5)/n}, the convention Miller (1982) used, and so reproduces his
    figures exactly.}
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
    \code{relative}, \code{mean_RT}, and the cell counts \code{n_AB},
    \code{n_A}, \code{n_B}.}
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
}
\seealso{
  \code{\link{rmiGroup.bayes}} for the hierarchical Bayesian analysis on the
  probability scale, \code{\link{build_rmi_bound_df}} for the figure that
  accompanies this test, and \code{\link{estimate.bounds}} for the bound curves
  themselves.
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
