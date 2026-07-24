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
  \code{resilience} returns a function that evaluates resilience capacity over time. 
  \code{resilience.test} returns a list with the test statistic, bounds, and inference details.
}
\details{
  Capacity resilience measures a system's robustness when shifting from simple parallel race baselines 
  to complex configurations. It evaluates the degree to which processing speed degrades (or improves) 
  when additional loads are imposed, relative to strict workload-capacity predictions.
}
