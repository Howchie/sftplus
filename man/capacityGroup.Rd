\name{capacityGroup}
\alias{capacityGroup}
\title{
Capacity Analysis
}
\description{
Performs workload capacity analysis on each participant and each condition.  Plots each capacity coefficient individually and returns the results of the nonparametric null-hypothesis test for unlimited capacity independent parallel performance.
}
\usage{
capacityGroup(inData, acc.cutoff=.9, ratio=TRUE, OR=NULL, 
  stopping.rule=c("OR", "AND", "STST"),
  plotCt=c("individual", "group", "both", "none"),
  trim=c(.05, .95), qtype=5, ...)
}
\arguments{
  \item{inData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{acc.cutoff}{Minimum accuracy for each stimulus category used in calculating the capacity coefficient.}
  \item{OR}{Indicates whether to compare performance to an OR or AND processing baseline.  Provided for backwards compatibility for package version < 2. }
  \item{stopping.rule}{Indicates whether to use OR, AND or Single Target Self Terminating (STST) processing baseline to calculate individual capacity functions.}
  \item{plotCt}{Which plots to draw.  \code{"individual"} gives one panel per participant, the historical behaviour and what \code{TRUE} still selects; \code{"group"} gives one panel per condition with every participant's capacity coefficient overlaid and their pointwise median drawn over the top; \code{"both"} gives both; \code{"none"} suppresses plotting, as \code{FALSE} still does.}
  \item{qtype}{Quantile algorithm used for \code{trim}, as in \code{\link{rmi.test}}.  The default, \code{5}, is the tie handling recommended by Ulrich et al. (2007).}
  \item{trim}{Two probabilities giving the quantile window over which each participant's capacity coefficient is reported, or \code{NULL} to report it wherever it is defined.  The window is the intersection of the per-cell windows, so it covers only the times at which every contributing estimate is supported by data.}
  \item{ratio}{Indicates whether to return the standard ratio capacity coefficient or, if FALSE, the difference form.}
  \item{...}{Arguments to be passed to plot function.}
}
\details{
For the details of the capacity coefficients, see \code{\link{capacity.or}}, \code{\link{capacity.and}} and \code{\link{capacity.stst}}.  If accuracy in any of the stimulus categories used to calculate a capacity coefficient falls below the cutoff, NA is retuned for that value in both the statistic and the Ct matrix.

The capacity coefficient is a ratio of cumulative hazard functions, and near the edges of a participant's own response times one of those functions has been estimated from almost no events.  Two things follow.  The coefficient is NA outside the range of a participant's response times, rather than being extrapolated from the nearest time at which it was defined, and \code{trim} further restricts it to the quantile window in which every contributing cell is supported.  \code{times} spans all participants and conditions, so most rows of \code{Ct.fn} are NA over part of the grid.  The grid itself runs over the union of the participants' windows rather than the \code{.001} and \code{.999} quantiles of the pooled response times, so that anticipations contributed by any one participant cannot stretch it below the range where any capacity coefficient is defined.

In the \code{"group"} plot each participant's curve stops where their own support does, so the median towards the edges of the axis is taken over fewer participants than in the middle.  The range of that count is printed above the panel.
}
\value{
  A list containing the following components:
  \item{overview}{A data frame indicating whether the OR and AND capacity coefficients significantly above baseline (super), below baseline (limited) or neither (unlimited) both at the individual level and at the group level in each condition.  NA is returned for any participant that had performance below the accuracy cutoff in a condition.}
  \item{Ct.fn}{Matrix with each row giving the values of the of the estimated capacity coefficient for one participant in one condition for values of times.  The rows match the ordering of statistic.}
  \item{Ct.var}{Matrix with each row giving the values of the of the variance of the estimated capacity coefficient for one participant in one condition for values of times.  Only returned if ratio=FALSE.}
  \item{capacity}{A list with the returned values from the capacity function for each participant and each condition.}
  \item{times}{Times at which the matrix capacity coefficients are calculated in Ct.fn matrix.}
  \item{specs}{A list of series specifications, one per participant and condition that yielded a fit, each with \code{label}, \code{Subject}, \code{Condition}, the capacity function \code{fn}, the redundant-target response times \code{rt} and correctness \code{cr}, and \code{benchmark}, that participant's correct response times pooled over cells.  These can be passed straight to \code{\link{build_ct_df}}; \code{\link{build_ct_prob_df}} also accepts the whole result.}
  \item{trim}{The quantile window used.}
  \item{qtype}{The quantile algorithm used.}
}
\references{
Townsend, J.T. & Wenger, M.J. (2004). A theory of interactive parallel processing: New capacity measures and predictions for a response time inequality series. \emph{Psychological Review, 111}, 1003--1035.

Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.

Houpt, J.W., Blaha, L.M., McIntire, J.P., Havig, P.R. and Townsend, J.T. (2013). Systems Factorial Technology with R. \emph{Behavior Research Methods}.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{capacity.and}}
\code{\link{capacity.or}}
\code{\link{capacity.stst}}
\code{\link{ucip.test}}
}
\examples{
\dontrun{
data(dots)
capacityGroup(subset(dots, Condition=="OR"), 
  stopping.rule="OR")
capacityGroup(subset(dots, Condition=="AND"), 
  stopping.rule="AND")
}
}
\keyword{ sft }
