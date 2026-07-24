\name{capacity_altieri}
\alias{capacity.altieri}
\title{Altieri's Capacity Coefficient for Accuracy and RT}
\description{
  Calculates Altieri's capacity coefficient, which incorporates both response time and accuracy. 
  It provides an omnibus measure of performance by integrating correctness directly into the 
  capacity evaluation.
}
\usage{
capacity.altieri(RT, CR = NULL, ratio = TRUE, Condition = NULL, Subject = NULL)
}
\arguments{
  \item{RT}{Response times in standard SFT list format, or a row-wise SFT data frame.}
  \item{CR}{Correctness indicators (logical/numeric) for each response time.}
  \item{ratio}{Logical; whether to return the standard ratio capacity (default) or additive capacity.}
  \item{Condition}{Optional string selecting one condition if \code{RT} is a data frame.}
  \item{Subject}{Optional string selecting one subject if \code{RT} is a data frame.}
}
\value{
  \code{capacity.altieri} returns a function that evaluates the Altieri capacity coefficient at given times, 
  along with associated metadata.
}
\details{
  Traditional Townsend capacity coefficient measures generally focus on response times of correct 
  responses, discarding incorrect trials or treating them separately. Altieri's extension unifies 
  accuracy and response time into a single measure, penalizing models that trade speed for accuracy.
}
