\name{sft_data_to_rt}
\alias{sft_data_to_rt}
\title{Convert row-wise SFT data to RT/CR list input}
\description{
Converts a canonical trial-level data frame into the response-time and
correctness lists accepted by the list-based SFT estimators.  For OR and AND
the returned cells are \code{AB}, \code{A}, and \code{B}; for STST they are
\code{context} and \code{target}.
}
\usage{
sft_data_to_rt(data, Condition = NULL, Subject = NULL,
  stopping.rule = c("OR", "AND", "STST"),
  by_subject = c("auto", "always", "never"), include_nn = FALSE)
}
\arguments{
  \item{data}{A data frame containing \code{Subject}, \code{RT},
    \code{Correct}, \code{Channel1}, and \code{Channel2}.  \code{Condition}
    is optional when the data contain one condition. The aliases \code{subjects},
    \code{rt}, and \code{LogicalRule} are accepted for \code{Subject},
    \code{RT}, and \code{Condition}, respectively.}
  \item{Condition}{Optional single condition value to retain.}
  \item{Subject}{Optional subject value or values to retain.}
  \item{stopping.rule}{One of \code{"OR"}, \code{"AND"}, or \code{"STST"}.}
  \item{by_subject}{Whether to return nested lists by subject.  The default
    \code{"auto"} returns a flat list for one subject and nested lists for
    multiple subjects.}
  \item{include_nn}{For OR/AND conversion, include the both-off \code{NN}
    cell.}
}
\value{
A list with components \code{RT} and \code{CR}.  The RT and CR components
have matching cell lists.  Metadata such as selected subjects, condition,
and dropped rows are stored as attributes.
}
\details{
Positive channel values are treated as target-present for OR and AND.  Rows
with both channels off are omitted unless \code{include_nn = TRUE}.  For STST,
the package's signed-channel convention is used: context rows have one
positive and one negative channel, while target rows have one nonzero
nonnegative channel.
}
\examples{
\dontrun{
  lists <- sft_data_to_rt(sftData, Condition = "OR", Subject = "P1")
  ucip.test(lists$RT, lists$CR, stopping.rule = "OR")
  capacityGroup.bayes(sftData, Condition = "OR")
}
}
\keyword{sft}
