\name{resilience}
\alias{resilience}
\alias{resilience.test}
\title{Capacity Resilience Assessment}
\description{
  Functions to compute and test capacity resilience. Resilience extends traditional capacity measures 
  by quantifying the system's ability to maintain processing efficiency under varied workload or 
  distractor conditions. \code{resilience.test} provides statistical bounds or tests for the estimated curve.
}
\usage{
resilience(RT, CR = NULL, ratio = TRUE, Condition = NULL, Subject = NULL)
resilience.test(RT, CR = NULL, Condition = NULL, Subject = NULL)
}
\arguments{
  \item{RT}{Response times in the standard SFT list format, or a row-wise SFT data frame.}
  \item{CR}{Optional correctness indicators.}
  \item{ratio}{Whether to return the ratio measure (default) or an additive alternative.}
  \item{Condition}{Optional string selecting one condition from the data frame.}
  \item{Subject}{Optional string selecting one subject from the data frame.}
}
\value{
  \code{resilience} returns \code{Rt} (and \code{Var} when \code{ratio = FALSE}) together
  with \code{Rtest}. \code{resilience.test} returns an \code{htest} whose statistic is the
  standard normal score statistic \eqn{z}, with the score components \code{numer},
  \code{variance}, and \code{denom} attached.
}
\details{
  Capacity resilience (Little, Eidels, Fific & Wang, 2015) extends the workload capacity
  coefficient to displays containing distracting information. It has the same functional
  form as the OR capacity coefficient,
  \deqn{R(t) = H_{AB}(t) / (H_{AY}(t) + H_{XB}(t)),}
  differing only in that the denominator conditions are single targets accompanied by a
  distractor (AY, XB) rather than lone targets. \code{RT} is supplied in that order.

  \code{resilience.test} is therefore the Houpt-Townsend (2012) weighted-logrank score
  statistic applied to that contrast, as Houpt and Little (2016) derived; it accumulates
  risk-set-weighted increments across the whole observed time course rather than reading a
  cumulative-hazard difference off at one time. That distinction matters: for complete data
  the Nelson-Aalen cumulative hazard at a sample's own maximum is exactly the harmonic
  number \eqn{\sum_{k=1}^{n} 1/k} whatever the distribution, so all three series coincide
  at the terminal time and a contrast evaluated there is the same constant for every
  dataset.
}
\references{
  Little, D. R., Eidels, A., Fific, M., & Wang, T. S. L. (2015). Understanding the influence
  of distractors on workload capacity. \emph{Journal of Mathematical Psychology},
  \bold{68--69}, 25--36.

  Houpt, J. W., & Little, D. R. (2016). Statistical analysis of the resilience function.
  \emph{Behavior Research Methods}.
}
