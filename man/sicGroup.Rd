\name{sicGroup}
\alias{sicGroup}
\title{
SIC Analysis for a Group
}
\description{
Calculates the SIC for each individual in each condition of a DFP experiment.  The function will plot each individuals SIC and return the results of the test for stochastic dominance and the statistical test of SIC form.
}
\usage{
sicGroup(inData, sictest=c("ks", "bf"), mictest=c("art", "anova"), domtest=c("ks", "dp"),
         alpha.sic=.05, plotSIC=TRUE, errors=c("discard", "defective"),
         min.trials=10L, perm=TRUE, nperm=1001L,
         stat=c("CvM", "KS", "AD", "Aly"), times=NULL, n.times=1000L,
         seed=NULL, ...)
}
\arguments{
  \item{inData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{sictest}{Which type of hypothesis test to use for SIC form. "ks" is the frequentist route and "bf" uses Bayesian SIC model probabilities.}
  \item{mictest}{Which type of hypothesis test to use for the MIC. The adjusted rank transform (art) and analysis of variance (anova) are the only tests currently implemented.}
  \item{domtest}{Which type of hypothesis test to use for stochastic dominance: series of KS tests ("ks") or a Bayesian Dirichlet-process-style route ("dp").}
  \item{alpha.sic}{Alpha level for determining a difference from zero used by the SIC overview.}
  \item{plotSIC}{Indicates whether or not to generate plots of the survivor interaction contrasts.}
  \item{errors}{\code{"discard"} conditions on accuracy by dropping incorrect responses, the default and the historical behaviour.  \code{"defective"} retains them and returns the unconditional interaction contrasts of \code{\link{sic}}.}
  \item{min.trials}{A participant is included in a condition only if every cell has more than this many usable trials.  Under \code{errors="discard"} that counts correct trials only; under \code{errors="defective"} it counts all trials.}
  \item{perm}{Whether to run the group-level permutation test of the interaction contrast.  Used only when \code{errors="defective"}.}
  \item{nperm}{Number of permutations for that test.}
  \item{stat}{Test statistic for that test: \code{"CvM"}, \code{"KS"}, \code{"AD"} or \code{"Aly"}.}
  \item{times}{Optional vector of times at which to evaluate every participant's SIC.  By default a grid is built from the observed response times.}
  \item{n.times}{Maximum size of that default grid.  When the data contain no more than \code{n.times} distinct response times they are all used; otherwise the grid is thinned to \code{n.times} evenly spaced points spanning the same range.}
  \item{seed}{Optional reproducibility seed, passed to the per-participant fits and the permutation test.}
  \item{...}{Arguments to be passed to plot function.}
}
\details{
%  TO DO:  Add details of sftData format
See the help page for the \code{\link{sic}} function for details of the survivor interaction contrast, and of the unconditional route selected by \code{errors="defective"}.

The default \code{times} grid is derived from the observed response times without rounding, so it is unaffected by the units the times are recorded in.  Seconds and milliseconds give the same curves up to the change of scale.

With \code{errors="defective"} no \code{Model} column is produced: architecture classification is not attempted on the unconditional contrasts.  Instead the interaction contrast itself is tested across participants with \code{\link{sicPermTest}}, separately for the correct, error and accuracy-agnostic subdistributions.  That test permutes whole participants between the two mixtures, so a condition with fewer than two admissible participants is skipped.
}
\value{
  \item{overview}{Data frame summarizing the test outcomes for each participant and condition.}
  \item{Subject}{The participant identifier from inData.}
  \item{Condition}{The condition identifier from inData.}
  \item{Selective.Influence}{The results of the survivor function dominance test for selective influence.  Pass indicates HH < HL, LH and LL > LH, HL, but not HL, LH < HH and not LH, HL > LL,  where A < B indicates that A is significantly faster than B at the level of the distribution.  Ambiguous means neither HL, LH < HH, nor LH, HL > LL, but at least one of HH < HL, LH or LL >HL, LH did not hold.  Fail means that at least one of HL, LH < HH or HL, LH > LL.}
  \item{Positive.SIC}{Indicates whehter the SIC is significantly positive at any time.}
  \item{Negative.SIC}{Indicates whehter the SIC is significantly negative at any time.}
  \item{MIC}{Indicates whether or not the MIC is significantly non-zero.}
  \item{Model}{Indicates which model would predict the pattern of data, assuming selective influence.}
  \item{SICfn}{Matrix with each row giving the values of the of the estimated SIC for one participant in one condition for values of times.  The rows match the ordering of statistic.}
  \item{sic}{List with each element giving the result applying sic() to an individual in a condition.  sic has the same ordering as overview.}
  \item{times}{Times at which the SICs in SICfn are calculated.}

  With \code{errors="defective"} the overview reports per-cell accuracy and the large-\eqn{t} asymptote of the correct and error contrasts in place of the SIC form columns, \code{sic.fn} holds the accuracy-agnostic contrast, and the result gains \code{sic.correct.fn}, \code{sic.error.fn} and \code{perm}, a list of \code{\link{sicPermTest}} results named by condition.
}
\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2010). The statistical properties of the survivor interaction contrast. \emph{Journal of Mathematical Psychology, 54}, 446-453.

Heathcote, A., Brown, S.D., Wagenmakers, E-J. & Eidels, A. (2010) Distribution-free tests of stochastic dominance for small samples. \emph{Journal of Mathematical Psychology, 54}, 454-463. 

Houpt, J.W., Blaha, L.M., McIntire, J.P., Havig, P.R. and Townsend, J.T. (2013). Systems Factorial Technology with R. \emph{Behavior Research Methods}.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{sic}}
\code{\link{sicPermTest}}
\code{\link{capacityGroup}}
}
\examples{
\dontrun{
data(dots)
sicGroup(dots)
}
}
\keyword{ ~sft }
