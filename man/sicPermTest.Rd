\name{sicPermTest}
\alias{sicPermTest}
\title{
Permutation Test of the Unconditional Interaction Contrast
}
\description{
Applies Angelov and Ekstrom's (2023) repeated-measures test for stochastic dominance to the interaction contrast of the unconditional response-time subdistributions, following Gondan and Legner.  The correct, error and accuracy-agnostic contrasts are tested separately.
}
\usage{
sicPermTest(inData, contrast = c("correct", "error", "agnostic"),
            stat = c("CvM", "KS", "AD", "Aly"), nperm = 1001L,
            alpha = .05, alpha.star = .8, gm = 1, balance = TRUE,
            min.trials = 10L, Condition = NULL, seed = NULL)
}
\arguments{
  \item{inData}{Canonical SFT trial data with \code{Subject}, \code{RT}, \code{Correct}, \code{Channel1} and \code{Channel2} columns.  A \code{Condition} column is used when \code{Condition} is supplied.}
  \item{contrast}{Which subdistributions to test.  Any subset of \code{"correct"}, \code{"error"} and \code{"agnostic"}; all three by default.}
  \item{stat}{Test statistic: Cramer-von Mises (\code{"CvM"}, the default), Kolmogorov-Smirnov (\code{"KS"}) or Anderson-Darling (\code{"AD"}) for the sign of the contrast, or Aly et al.'s slope statistic (\code{"Aly"}) for the monotonicity predictions.  See Details, as the two families report different quantities.}
  \item{nperm}{Number of permutations.}
  \item{alpha}{Significance level.}
  \item{alpha.star}{Tuning parameter separating crossing from dominating distribution functions.}
  \item{gm}{Weight exponent for the Anderson-Darling statistic; \code{gm = 1} gives the usual \eqn{1/(H(1-H))} weight.}
  \item{balance}{Subsample each participant's four cells to a common size before forming the mixtures.  Leave \code{TRUE} unless the design is already balanced.}
  \item{min.trials}{A participant is included only if every cell has more than this many trials.}
  \item{Condition}{Optional condition to subset before testing.}
  \item{seed}{Optional reproducibility seed.}
}
\details{
Each participant contributes two mixtures, \eqn{F_X = (F_{HH} + F_{LL})/2} and \eqn{F_Y = (F_{HL} + F_{LH})/2}, whose difference is half the interaction contrast.  The null distribution is obtained by randomly swapping the two mixture labels within each participant, so this is a group-level test: with a single participant the permutation distribution has only two points and the p-values are uninformative.

Responses of the other type, and any omissions, are set to \code{+Inf} so that they stay in the denominator without contributing at any finite time.  The half-and-half mixture assumes equal trial counts per cell, which \code{balance} enforces by subsampling to the participant's smallest cell.

Signs are reported in this package's survivor convention, \code{SIC = -D2F}, so \code{"SIC positive"} corresponds to Gondan and Legner's "IC negative".

The Cramer-von Mises, Kolmogorov-Smirnov and Anderson-Darling statistics compare the \emph{level} of the two mixtures and so resolve the sign of the interaction contrast.  They report \code{p.SICpos} and \code{p.SICneg} and decide among \code{"n.s."}, \code{"SIC positive"}, \code{"SIC negative"} and \code{"SIC crosses"}.

Aly's slope statistic compares the mixtures' \emph{increments} instead, and so resolves monotonicity rather than sign.  That is the form the subdensity predictions take, for instance the prediction that the error subdensity contrast has one sign at every \eqn{t}, which makes the error contrast monotone.  It reports \code{p.SICinc} and \code{p.SICdec} and decides among \code{"n.s."}, \code{"SIC increasing"}, \code{"SIC decreasing"} and \code{"SIC non-monotone"}.

The four-way decision follows Angelov and Ekstrom's published rule, which is deliberately asymmetric: \code{"SIC positive"} requires \code{p.SICpos <= alpha} together with \code{p.SICneg > alpha}, whereas \code{"SIC negative"} requires \code{p.SICneg <= alpha} together with \code{p.SICpos > alpha.star}.  A contrast with strong one-sided evidence can therefore still be labelled \code{"SIC crosses"} when the opposing p-value falls between \code{alpha} and \code{alpha.star}.  Both one-sided p-values are returned so that a different rule can be applied.
}
\value{
A list with one \code{htest} per requested contrast, each carrying \code{decision} and the two one-sided p-values alongside the usual fields.  The p-values are named \code{p.SICpos} and \code{p.SICneg} for the sign statistics and \code{p.SICinc} and \code{p.SICdec} for \code{stat = "Aly"}.  The list also contains:
 \item{summary}{Data frame of both p-values and the decision for each contrast.  The two p-value columns are named for whichever question the chosen statistic answers.}
 \item{subjects}{Participants contributing to the test.}
 \item{stat}{The statistic used.}
 \item{nperm}{Number of permutations.}
}
\references{
Angelov, A.G. & Ekstrom, M. (2023). Tests of stochastic dominance with repeated measurements data. \emph{AStA Advances in Statistical Analysis, 107}, 443-467.

Aly, E.-E.A.A., Kochar, S.C. & McKeague, I.W. (1994). Some tests for comparing cumulative incidence functions and cause-specific hazard rates. \emph{Journal of the American Statistical Association, 89}, 994-999.

Gondan, M. & Legner, L. Systems factorial technology: unconditional interaction contrasts for correct and incorrect responses.
}
\seealso{
\code{\link{sic}}
\code{\link{sicGroup}}
\code{\link{siDominance}}
}
\examples{
\dontrun{
data(dots)
sicPermTest(dots)
}
}
\keyword{ sft }
