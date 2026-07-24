\name{assessment_extensions}
\alias{assessment_ta_or}
\alias{assessment_ta_and}
\alias{assessment_donkin_discrimination}
\title{Timed-assessment and Donkin Discrimination Capacity Extensions}
\description{
  Functions to calculate capacity bounds and assessments for timed-assessment (trial-abundance) models
  and Donkin discrimination designs. \code{assessment_ta_or} and \code{assessment_ta_and} evaluate the log-ratio 
  of predicted to observed cumulative probability metrics (correct and incorrect, fast and slow) under 
  OR and AND logical rules. \code{assessment_donkin_discrimination} implements a capacity assessment 
  for the Donkin discrimination paradigm by integrating the weighted survivor function of one channel 
  with respect to the correct/incorrect latency of the other.
}
\usage{
assessment_ta_or(RT, CR = NULL, times = NULL)
assessment_ta_and(RT, CR = NULL, times = NULL)
assessment_donkin_discrimination(RT, CR = NULL, times = NULL)
}
\arguments{
  \item{RT}{Response times in the standard SFT list format: \code{list(AB, A, B)}.}
  \item{CR}{Optional correctness indicators in the standard SFT list format, same length as \code{RT}.}
  \item{times}{Evaluation times for the resulting assessment curves. If \code{NULL}, times are inferred from finite RTs.}
}
\value{
  A list containing evaluation times and predicted vs. observed functions (CF = correct fast, 
  CS = correct slow, IF = incorrect fast, IS = incorrect slow) as well as the log capacity ratios 
  (\code{A_CF}, \code{A_CS}, \code{A_IF}, \code{A_IS}).
}
\details{
  The timed-assessment (trial-abundance) architecture partitions responses into correct fast, correct slow, 
  incorrect fast, and incorrect slow. The \code{assessment_ta_or} and \code{assessment_ta_and} functions compute 
  these defective CDFs and survivor functions for single targets and use them to predict the double-target performance 
  under independent parallel race models for OR (first-terminating) and AND (exhaustive) tasks. 
  \code{assessment_donkin_discrimination} uses the Donkin framework to predict the joint target conditions 
  based on integrals of the marginal RT distributions.
}
