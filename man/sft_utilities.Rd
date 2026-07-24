\name{sft_utilities}
\alias{correctfun}
\alias{get_time_range}
\alias{smooth_one_cdf}
\title{General SFT Utilities}
\description{
  Miscellaneous utilities to assist in correctness coding, temporal grid construction, 
  and smooth CDF estimation for Systems Factorial Technology workflows.
}
\usage{
correctfun(d)
get_time_range(rt_list, probs = c(.05, .95), by = .001)
smooth_one_cdf(t, y, smooth_cdf = "none", smooth_spar = .65)
}
\arguments{
  \item{d}{A simulated data frame for correctness coding.}
  \item{rt_list}{A list of numeric response-time vectors.}
  \item{probs}{Probability bounds used for trimming the empirical time range.}
  \item{by}{Step size for the generated time grid.}
  \item{t, y}{Time and cumulative distribution values to be smoothed.}
  \item{smooth_cdf}{Smoothing method: \code{"none"}, \code{"mono"}, or \code{"spline"}.}
  \item{smooth_spar}{Smoothing parameter for spline fitting.}
}
\value{
  \code{correctfun} returns a logical vector of correctness. 
  \code{get_time_range} returns a numeric sequence of times. 
  \code{smooth_one_cdf} returns a smooth function approximating the cumulative distribution.
}
\details{
  \code{get_time_range} simplifies defining a common evaluation grid across diverse experimental 
  conditions, eliminating extreme outliers using the \code{probs} quantiles. \code{smooth_one_cdf} 
  applies monotonic or spline constraints to empirical CDF steps to produce a differentiable output 
  frequently used in hazard-rate computations.
}
