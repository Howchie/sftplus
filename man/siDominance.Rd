\name{siDominance}
\alias{siDominance}
\title{
Dominance Test for Selective Influence
}
\description{
Function to test for the survivor function ordering predicted by the selective influence of the salience manipulation.
}
    
\usage{
siDominance(HH, HL, LH, LL, method=c("ks", "dp"), nbin=20L, nsamp=10000L, seed=NULL,
            CR=NULL, errors=c("discard", "defective"))
}
\arguments{
  \item{HH}{Response times from the High--High condition.}
  \item{HL}{Response times from the High--Low condition.}
  \item{LH}{Response times from the Low--High condition.}
  \item{LL}{Response times from the Low--Low condition.}
  \item{method}{Which type of hypothesis test to use for testing stochastic dominance relations, either as series of KS tests ("ks") or the Bayesian Dirichlet-process-style test ("dp").}
  \item{nbin}{Number of common RT bins for the Bayesian route.}
  \item{nsamp}{Number of prior/posterior Monte Carlo draws for the Bayesian route.}
  \item{seed}{Optional reproducibility seed.}
  \item{CR}{Optional list of four correctness vectors in HH, HL, LH, LL order (or named with those cells), one entry per trial of the matching response-time cell.  Required when \code{errors="defective"}.}
  \item{errors}{Whether to \code{"discard"} incorrect responses and condition on accuracy (the default, and the historical behaviour) or to retain them and test the \code{"defective"} subdistributions.}
}
\details{
The unconditional route (\code{errors="defective"}) tests the extended assumption set that Gondan and Legner require for their generalisation of the theory to tasks with non-negligible error rates.  Beyond the survivor ordering, higher salience must produce more fast \emph{and correct} responses (8A) and fewer fast \emph{and incorrect} ones (8D).  These are not automatic: they fail for evidence-accumulation models biased towards the incorrect alternative, and for any manipulation that speeds correct and incorrect responses alike.  Assumption 8C, on the error subdensities, is not directly testable here; 8D is its testable consequence and is what the derivations use.

The two-sample Kolmogorov-Smirnov comparisons are computed on defective samples, in which responses of the other type and any omissions are placed at \code{+Inf} so that they remain in the denominator without contributing at any finite time.  The shared atom at \code{+Inf} does not affect the supremum, but it does make the asymptotic p-value conservative.

For an SIC function to distinguish among the processing types, the salience manipulation on each channel must selectively influence its respective channel (although see Eidels, Houpt, Altieri, Pei & Townsend, 2010 for SIC prediction from interactive parallel models).  Although the selective influence assumption cannot be directly tested, one implication is that the distribution the HH response times stochastically dominates the HL and LH distributions which each in turn stochastically dominate the LL response time distribution.  This implication is automatically tested in this function.  The KS dominance test uses eight two-sample Kolmogorov-Smirnov tests:  HH < HL, HH < LH, HL < LL, LH < LL should be significant while HH > HL, HH > LH, HL > LL, LH > LL should not.  The DP uses four tests to determine which relation has the highest Bayes factor assuming a Dirichlet process prior for each of (HH, HL), (HH, LH), (HL, LL) and (LH, LL).  See Heathcote, Brown, Wagenmakers & Eidels, 2010, for more details.

}
\value{
 A data frame with the first column indicating which ordering was tested, the second column indicating the test statistic and the third indicating the p-value for that value of the statistic.

 With \code{errors="defective"} the frame gains an \code{Assumption} column naming the assumption under test and a logical \code{Predicted} column that is \code{TRUE} for the ordering selective influence predicts.  Three blocks are reported: survivor ordering (Gondan & Legner's 8B, the assumption inherited from the original theory), ordering of the correct subdistributions (8A), and the reverse ordering of the error subdistributions (8D).  Under selective influence the \code{Predicted} rows should be significant and their reversals should not.
}
\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2010). The statistical properties of the survivor interaction contrast. \emph{Journal of Mathematical Psychology, 54}, 446-453.

Dzhafarov, E.N., Schweickert, R., & Sung, K. (2004). Mental architectures with selectively influenced but stochastically interdependent components. \emph{Journal of Mathematical Psychology}, 48, 51-64.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{ks.test}}
\code{\link{sic}}
\code{\link{sicGroup}}
\code{\link{mic.test}}
}
\examples{
T1.h <- rexp(50, .2)
T1.l <- rexp(50, .1)
T2.h <- rexp(50, .21)
T2.l <- rexp(50, .11)

HH <- T1.h + T2.h
HL <- T1.h + T2.l
LH <- T1.l + T2.h
LL <- T1.l + T2.l
siDominance(HH, HL, LH, LL)
}
\keyword{ sft }
