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
sicDPtest(dat, nbin = NULL, nsamp = 10000L, maxn = 500000L, tolSIC = NULL,
          tolMIC = NULL, ci = .95, seed = NULL, CR = NULL,
          errors = c("discard", "defective"))
sictestBayes(HH, HL, LH, LL, method = c("DP", "IG"), model = NULL, nbin = NULL,
             nsamp = 10000L, maxn = 500000L, seed = NULL, CR = NULL,
             errors = c("discard", "defective"), ...)
}
\arguments{
  \item{inData}{Standard SFT data frame.}
  \item{domtest}{Dominance test method (e.g., \code{"ks"} for Kolmogorov-Smirnov).}
  \item{plotSIC}{Accepted for backward compatibility and currently unused; no plot is drawn.}
  \item{dat}{A four-element list of SIC RT samples in HH, HL, LH, LL order for Dirichlet process testing.}
  \item{nbin}{Optional number of bins for Bayesian SIC/dominance calculations.}
  \item{nsamp, maxn}{Monte Carlo size bounds (initial samples and maximum iterations).}
  \item{tolSIC, tolMIC}{Tolerances for model classification of the SIC and Mean Interaction Contrast (MIC).  \code{tolSIC} defaults to 0.05 for the correct-only route and to 0.01 for the unconditional route, where the error and accuracy-agnostic contrasts are bounded by accuracy differences and so are an order of magnitude smaller.  \code{tolMIC} is unused when \code{errors="defective"}.}
  \item{CR}{Optional list of four correctness vectors in HH, HL, LH, LL order, one entry per trial of the matching response-time cell.  Required when \code{errors="defective"}.}
  \item{errors}{\code{"discard"} conditions on accuracy by dropping incorrect responses, the default and the historical behaviour.  \code{"defective"} retains them and works with the unconditional subdistributions.}
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

  With \code{errors="defective"} the six architecture classes are not produced.  \code{BF},
  \code{BFp} and \code{posterior} are instead lists with one entry per contrast
  (\code{correct}, \code{error} and \code{all}), each a vector over the four sign classes
  \code{Zero}, \code{Negative}, \code{Positive} and \code{Crosses}.
}
\details{
  The SIC is fundamental in SFT for distinguishing cognitive architectures and stopping rules.

  Under \code{errors="defective"} each cell instead carries a single Dirichlet over
  \eqn{2(nbin + 1)} categories: \code{nbin} bins of correct response times plus a catch-all for
  correct responses never observed, then the same for errors.  Splitting the unit concentration
  evenly between the two halves makes the implied prior on cell accuracy uniform, which a single
  flat Dirichlet over all categories would not.  Each draw yields the correct, error and
  accuracy-agnostic contrasts, whose signs are summarised but not mapped onto an architecture; see
  \code{\link{sic}} for why.
  The Bayesian Dirichlet process test gives each of the four cells its own Dirichlet over a shared
  set of pooled-RT bins and forms the contrast
  \eqn{SIC = F_{LH} + F_{HL} - F_{HH} - F_{LL}} draw by draw, the same definition
  \code{\link{sic.test}} uses. Each draw is classified into one of six shape classes and the Bayes
  factor for a class is its posterior-to-prior odds.

  The classes are \code{Z} (everywhere within \code{tolSIC} of zero: serial-OR),
  \code{N} (all-negative: parallel exhaustive), \code{P} (all-positive: parallel self-terminating,
  and coactive when its early negative dip falls below \code{tolSIC}), and the negative-then-positive
  family \code{np}/\code{nP}/\code{Np}, distinguished by whether the signed areas balance to within
  \code{tolMIC} (\code{np}: serial exhaustive) or favour the positive (\code{nP}) or negative
  (\code{Np}) lobe.

  \code{tolSIC} is applied to every sign test, so \code{N} and \code{P} tolerate excursions of up to
  \code{tolSIC} on the other side rather than demanding literal sign purity; without that, no noisy
  posterior draw could ever reach them. Raising \code{tolSIC} merges shallow lobes into the
  neighbouring pure class, which is the main lever on how coactive curves are reported.
}
