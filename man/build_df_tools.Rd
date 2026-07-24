\name{build_df_tools}
\alias{build_at_df}
\alias{build_cdf_df}
\alias{build_ct_df}
\alias{build_sic_df}
\title{Data Frame Builders for Capacity and SIC Curves}
\description{
  Utility functions for transforming raw or computed Systems Factorial Technology (SFT) curves 
  into tidy data frames suitable for plotting (e.g., with ggplot2). They handle smoothing, 
  trimming based on probabilities, and proper formatting of times and values.
}
\usage{
build_at_df(times, series_specs, probs = c(.05, .95))
build_cdf_df(times, series_specs, probs = c(.05, .95), smooth_cdf = c("none", "mono", "spline"), smooth_spar = .65)
build_ct_df(times, series_specs, probs = c(.05, .95))
build_sic_df(times, series_specs, trim_cdf_tails = TRUE, cdf_tail_cut = 5e-4, smooth_cdf = c("none", "mono", "spline"), smooth_spar = .65)
}
\arguments{
  \item{times}{Numeric vector of evaluation times for the curves.}
  \item{series_specs}{A list of named curve/data specifications defining the functions to be evaluated.}
  \item{probs}{Quantile probabilities used to trim plotted curves to avoid extreme long tails.}
  \item{smooth_cdf}{Method for smoothing the CDFs: \code{"none"}, \code{"mono"} (monotonic), or \code{"spline"}.}
  \item{smooth_spar}{Spline smoothness parameter, used if \code{smooth_cdf} is \code{"spline"}.}
  \item{trim_cdf_tails}{Logical; whether to trim common CDF tails for the SIC curves.}
  \item{cdf_tail_cut}{Numeric cut-off for trimming CDF tails aggressively.}
}
\value{
  A tidy \code{data.frame} containing columns for \code{time}, \code{value}, and grouping identifiers 
  derived from \code{series_specs}, ready for visualization.
}
\details{
  These builders are primarily used internally by SFT plotting functions to construct smooth, 
  trimmed paths for capacity coefficients (\code{build_ct_df}, \code{build_at_df}), cumulative distribution 
  functions (\code{build_cdf_df}), and survivor interaction contrasts (\code{build_sic_df}). 
  Trimming avoids visually unappealing long, flat tails in empirical distribution functions.
}
