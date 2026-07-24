\name{mic.bayes}
\alias{mic.bayes}
\title{Hierarchical posterior mean interaction contrast (MIC)}
\description{
  Extracts the posterior mean interaction contrast from an OR-centred
  hierarchical Bayesian capacity model fitted with a salience split. The MIC is
  \eqn{\mathrm{mean}(LL) - \mathrm{mean}(LH) - \mathrm{mean}(HL) +
  \mathrm{mean}(HH)}, equivalently the signed area under the survivor
  interaction contrast \eqn{\int \mathrm{SIC}(t)\,dt}. It is evaluated draw by
  draw at the subject, population (typical-subject parameter),
  \code{population_finite_mean}, and \code{new_subject} levels, so its sign can
  be tested against zero at every level, and each subject is contrasted against
  the group. The MIC is only identified from a salience-split fit.
}
\usage{
mic.bayes(object, ...)
}
\arguments{
  \item{object}{A fitted \code{sft_bayes} object, or a canonical SFT
    data frame to fit with \code{\link{semiparametricSFT.bayes}}.}
  \item{...}{Passed to \code{\link{semiparametricSFT.bayes}} when \code{object} is a data
    frame. A \code{salience_split} is required to identify the four cells and
    defaults to \code{TRUE}.}
}
\value{
  A list with components:
  \item{summary}{A data frame with one row per level (\code{population},
    \code{population_finite_mean}, \code{new_subject}, and each subject) giving
    the posterior \code{Mean}, highest-density interval (\code{Lower},
    \code{Upper}), and sign tests \code{Prob_positive} \eqn{= P(\mathrm{MIC} >
    0)} and \code{Prob_negative} \eqn{= P(\mathrm{MIC} < 0)}.}
  \item{contrasts}{Subject-versus-population MIC differences with the same sign
    tests, identifying subjects whose architecture departs from the group.}
  \item{draws}{The raw posterior MIC draws at each level.}
  \item{group_estimand, definition, interpretation}{Metadata describing the
    reported group estimand and the over-additive / under-additive / additive
    (parallel-OR or coactive / parallel-AND / serial) interpretation.}
}
\details{
  \code{Prob_positive} near one indicates an over-additive MIC (parallel
  self-terminating or coactive processing); \code{Prob_negative} near one
  indicates an under-additive MIC (parallel exhaustive processing); posterior
  mass concentrated near zero is consistent with additive, serial processing.
  The \code{population} level is the typical-subject parameter obtained with the
  subject random effects held at zero, a genuine hierarchical group estimand
  rather than an average of independently fitted subjects.
}
\seealso{\code{\link{semiparametricSFT.bayes}}, \code{\link{sic.bayes}},
  \code{\link{mic.test}}}
