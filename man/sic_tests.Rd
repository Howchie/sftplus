\name{sic_tests}
\alias{sicGroupBF}
\alias{sicDPtest}
\alias{sictestBayes}
\title{Survivor Interaction Contrast (SIC) Tests}
\description{
  Functions for performing tests on the Survivor Interaction Contrast (SIC). 
  \code{sicGroupBF} provides group-level Bayes Factor analysis for dominance tests.
  \code{sicDPtest} uses a Bayesian Dirichlet Process framework to test SIC model classes, 
  and \code{sictestBayes} routes between Dirichlet Process and Inverse Gamma Bayesian methods 
  for rigorous class inference.
}
\usage{
sicGroupBF(inData, domtest = "ks", plotSIC = TRUE, ...)
sicDPtest(dat, nbin = NULL, nsamp = 10000L, maxn = 500000L, tolSIC = 5e-2, tolMIC = NULL, ci = .95, seed = NULL)
sictestBayes(HH, HL, LH, LL, method = c("DP", "IG"), model = NULL, nbin = NULL, nsamp = 10000L, maxn = 500000L, seed = NULL, ...)
}
\arguments{
  \item{inData}{Standard SFT data frame.}
  \item{domtest}{Dominance test method (e.g., \code{"ks"} for Kolmogorov-Smirnov).}
  \item{plotSIC}{Logical; whether to automatically plot the resulting SIC.}
  \item{dat}{A four-element list of SIC RT samples in HH, HL, LH, LL order for Dirichlet process testing.}
  \item{nbin}{Optional number of bins for Bayesian SIC/dominance calculations.}
  \item{nsamp, maxn}{Monte Carlo size bounds (initial samples and maximum iterations).}
  \item{tolSIC, tolMIC}{Tolerances for model classification of the SIC and Mean Interaction Contrast (MIC).}
  \item{ci}{Monte Carlo interval mass.}
  \item{seed}{Random seed for reproducibility.}
  \item{HH, HL, LH, LL}{Numeric vectors of response times for the four salient factor cells (High/High, High/Low, Low/High, Low/Low).}
  \item{method}{Bayesian method selection: \code{"DP"} for Dirichlet Process (default) or \code{"IG"} for Inverse Gamma.}
  \item{model}{Optional external Bayesian model object.}
  \item{...}{Additional Bayesian-control or plotting arguments.}
}
\value{
  A test object containing Bayes Factors, posterior class probabilities, or dominance test 
  statistics, typically allowing for the identification of the underlying processing architecture 
  (parallel, serial, coactive, etc.).
}
\details{
  The SIC is fundamental in SFT for distinguishing cognitive architectures and stopping rules. 
  The Bayesian Dirichlet Process tests provide a non-parametric approach to estimating the 
  probabilities of different architectural predictions (e.g., all-positive for coactive, 
  positive-then-negative for parallel exhaustive).
}
