\name{build_ct_prob_df}
\alias{build_ct_prob_df}
\alias{build_cdf_prob_df}
\alias{build_prob_summary_df}
\title{
Probability-Scale Curve Builders
}
\description{
Evaluate each participant's capacity coefficient or cumulative distribution function at that participant's own response-time quantiles, so that curves can be averaged across a group at a common position in the distribution rather than at a common clock time.
}
\usage{
build_ct_prob_df(x, probs = seq(.05, .95, by = .05), condition = NULL,
                 benchmark = c("series", "all"), qtype = 5)

build_cdf_prob_df(x, probs = seq(.05, .95, by = .05), condition = NULL,
                  benchmark = c("all", "series"),
                  smooth_cdf = c("none", "mono", "spline"),
                  smooth_spar = .65, qtype = 5)

build_prob_summary_df(df, value = c("Ct", "CDF"), center = c("median", "mean"),
                      ci_probs = c(.1, .9), by = "Condition")
}
\arguments{
  \item{x}{Either a \code{\link{capacityGroup}} result, or a list of series specifications.  Each specification is a list with \code{label}, a curve \code{fn}, the response times \code{rt} the series is built from and optionally their correctness \code{cr}, optionally \code{Subject} and \code{Condition} labels, and optionally \code{benchmark}, a wider set of response times for that participant.}
  \item{probs}{Percentile levels, strictly inside \eqn{(0, 1)}, as in \code{\link{rmi.test}}.}
  \item{condition}{Optional conditions to keep.}
  \item{benchmark}{Which response times define a participant's percentiles.  \code{"series"} uses the series' own times; \code{"all"} uses \code{spec$benchmark} when present, so that several series belonging to one participant share a time axis.}
  \item{qtype}{Quantile algorithm, as in \code{\link{rmi.test}}.  The default, \code{5}, is the tie handling recommended by Ulrich et al. (2007) and is the same estimator the redundant-signals functions use, so the two probability-scale analyses sit on a common footing.}
  \item{smooth_cdf, smooth_spar}{Optional smoothing of each participant's distribution function, as in \code{\link{smooth_one_cdf}}.}
  \item{df}{A data frame returned by one of the builders.}
  \item{value}{Name of the column to summarise; by default whichever of \code{"Ct"} or \code{"CDF"} is present.}
  \item{center}{Central tendency across series.  The median is the default because capacity coefficients are ratios and occasionally take extreme values.}
  \item{ci_probs}{Two probabilities giving the spread reported across series.}
  \item{by}{Extra grouping columns, used when present in \code{df}.}
}
\details{
Averaging \eqn{C_i(t)} across participants at a fixed \eqn{t} estimates a quantity that few participants contribute to.  At an early \eqn{t} the fast participants are midway through their distribution while the slow ones have not responded at all, and their coefficients are undefined there; at a late \eqn{t} the situation is reversed.  The group curve is then an average over a changing and self-selected subset.

Reading each participant at \eqn{Q_i(p)} instead estimates \eqn{E[C_i(Q_i(p))]}, to which every participant contributes at every \eqn{p}.  This is legitimate to inspect, but it is a different functional and should be reported as such.  It is not the capacity coefficient of the pooled data, and it is not \eqn{C(t)} at any particular time; the \code{Time} column gives the mean of the participants' \eqn{Q_i(p)} only as a rough guide to where on the clock a given percentile sits.

This is the same change of index that \code{\link{rmiGroup.bayes}} makes for the race-model inequality, and the builders share its vocabulary: the percentile column is \code{prob}, the levels are \code{probs}, the estimator is selected by \code{qtype}, and \code{build_prob_summary_df} returns \code{n}, \code{mean}, \code{lower} and \code{upper} with a \code{scale} of \code{"probability"}, exactly as \code{\link{build_rmi_violation_df}} does.  Plotting code written for one serves the other.

Two features make the alignment better behaved for capacity coefficients than the corresponding Vincent average of raw response times.  Each \eqn{C_i} has already been normalised against that participant's own UCIP prediction, so the reference value of 1 is common to everyone and is preserved by any pointwise average.  And the shape distortion that Vincent averaging introduces acts on the time axis, which is exactly the axis being replaced.  What the alignment cannot repair is heterogeneity in architecture: averaging a limited-capacity and a super-capacity participant returns something near 1 at every percentile, just as it would at every \eqn{t}.

For distribution functions the default \code{benchmark = "all"} takes the percentile grid from a participant's pooled trials, so the conditions of one participant are read on one time axis and the ordering between them is preserved.  Using each condition's own percentiles instead would map every condition onto the identity in \code{prob} and destroy the comparison.
}
\value{
For the builders, a long data frame with \code{prob}, \code{Time}, the curve value (\code{Ct} or \code{CDF}), \code{Series}, \code{scale}, and any \code{Subject} and \code{Condition} labels the specifications carried.  For \code{build_prob_summary_df}, one row per \code{prob} (and per level of \code{by}) with \code{n}, the mean \code{Time}, \code{mean}, \code{lower}, \code{upper}, \code{center} and \code{scale}.
}
\references{
Ulrich, R., Miller, J. & Schroter, H. (2007). Testing the race model inequality: An algorithm and computer programs. \emph{Behavior Research Methods, 39}, 291-302.
}
\seealso{
\code{\link{build_ct_df}}
\code{\link{build_cdf_df}}
\code{\link{build_rmi_violation_df}}
\code{\link{capacityGroup}}
}
\examples{
\dontrun{
data(dots)
cg <- capacityGroup(dots, plotCt = "group")
q <- build_ct_prob_df(cg)
build_prob_summary_df(q)
}
}
\keyword{ sft }
