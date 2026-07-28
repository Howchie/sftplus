\name{build_rmi_df}
\alias{build_rmi_cdf_df}
\alias{build_rmi_bound_df}
\alias{build_rmi_violation_df}
\title{Tidy Data Frames for Race Model Inequality Figures}
\description{
  \code{build_rmi_bound_df} produces the classical race-model figure: the
  redundant-target cumulative distribution plotted against the bound it must
  not cross, so a violation is visible as the AB curve lying to the left of the
  \eqn{F_A + F_B} bound. It is the visual complement to the per-percentile
  t-tests of \code{\link{rmi.test}}, which report where the gap is reliable but
  not how large it is. \code{build_rmi_cdf_df} gives the three condition CDFs
  behind it (Miller, 1982, Fig. 1), and \code{build_rmi_violation_df} the
  violation curve with its interval.
}
\usage{
build_rmi_cdf_df(x, bound = NULL, condition = NULL, by_block = FALSE,
                 qtype = 5, probs = NULL,
                 cdf_method = c("polygon", "kernel"), kernel_bw = "nrd0")
build_rmi_bound_df(x, bound = NULL, condition = NULL, by_block = FALSE,
                   qtype = 5, probs = NULL,
                   cdf_method = c("polygon", "kernel"), kernel_bw = "nrd0")
build_rmi_violation_df(x, bound = NULL, condition = NULL,
                       scale = c("auto", "probability", "time"),
                       by_block = FALSE, qtype = 5, probs = NULL,
                       cdf_method = c("polygon", "kernel"),
                       kernel_bw = "nrd0")
}
\arguments{
  \item{x}{The result of \code{\link{rmi.test}} or
    \code{\link{rmiGroup.bayes}}, or any single per-bound / per-condition
    element of one, or a raw trial data frame.}
  \item{bound}{Optional bound to keep when the result holds both.}
  \item{condition}{Optional condition(s) to keep.}
  \item{probs}{Optional percentile levels. For raw trial data these are passed
  to \code{rmi.test}; for an existing result they select levels already
  present in that result. If omitted, the \code{rmi.test} default is used for
  raw data and all levels are retained from an existing result.}
  \item{scale}{For \code{build_rmi_violation_df}: \code{"auto"} (the
    probability scale if the fit is Bayesian, otherwise time),
    \code{"probability"}, or \code{"time"}.}
  \item{by_block}{Passed to \code{\link{rmi.test}} when \code{x} is raw
    trial data: estimate the quantiles per participant and block rather than
    per participant.}
  \item{qtype}{Plotting-position rule passed to \code{\link{rmi.test}} when
    \code{x} is raw trial data; the default 5 is Equation 2 of Ulrich, Miller
    and Schroeter (2007).}
  \item{cdf_method, kernel_bw}{CDF estimator and kernel bandwidth passed to
    \code{\link{rmi.test}} when \code{x} is raw trial data.}
}
\details{
  All three return a long data frame ready for \code{ggplot2}, following the
  same convention as \code{\link{build_cdf_df}} and the other builders in this
  package: no plot object is constructed and \pkg{ggplot2} is not a dependency.

  The quantiles in \code{build_rmi_cdf_df} and \code{build_rmi_bound_df} are
  averaged across participants at each percentile level --- Miller's
  "Vincentized" group curves. That average is on the millisecond scale and so
  is base-time dependent. It is a display convention: the quantity actually
  tested by \code{\link{rmiGroup.bayes}} is the within-participant
  probability-scale margin returned by \code{build_rmi_violation_df}. Read the
  first two as the picture and the third as the inference.

  \code{build_rmi_violation_df} reports the probability-scale population curve
  with its highest-density interval and \code{p_violation} when given a
  Bayesian fit, and the time-scale mean difference with its t-interval and
  adjusted p-value when given \code{\link{rmi.test}} output. Negative values
  are violations on both scales.
}
\value{
  A long \code{data.frame}. \code{build_rmi_cdf_df} has columns \code{bound},
  \code{condition}, \code{prob}, \code{Condition} (a factor with levels
  \code{A}, \code{B}, \code{AB}), and \code{RT}. \code{build_rmi_bound_df}
  replaces \code{Condition} with \code{curve}, whose levels are \code{"AB"} and
  the bound. \code{build_rmi_violation_df} carries a \code{scale} column
  identifying which of the two summaries it holds.
}
\references{
  Miller, J. (1982). Divided attention: Evidence for coactivation with
  redundant signals. \emph{Cognitive Psychology}, 14, 247--279.
}
\seealso{
  \code{\link{rmi.test}}, \code{\link{rmiGroup.bayes}},
  \code{\link{build_cdf_df}}.
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

# The race-model figure: AB against the A + B bound.
d <- build_rmi_bound_df(fit)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  ggplot(d, aes(RT, prob, linetype = curve)) +
    geom_line(linewidth = 1.2) + geom_point(size = 2) +
    labs(x = "Response time (ms)", y = "Cumulative probability",
         linetype = NULL) +
    theme_classic()
}
}
