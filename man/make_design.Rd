\name{make_design}
\alias{make_design}
\title{Generate a MaxDiff experimental design}
\description{
  Creates a balanced MaxDiff design using \code{choiceDes::tradeoff.des()}.
  The design is structured so that every item appears approximately the same
  number of times and in each position, and no pair of items appears together
  too frequently. After initial generation, \code{make_design()} applies
  post-processing to repair any forbidden-pair violations and to improve
  position balance.
}
\usage{
make_design(
  n_items,
  n_versions       = 100,
  n_tasks          = 12,
  n_items_per_task = 4,
  forbidden_pairs  = NULL,
  seed             = 42
)
}
\arguments{
  \item{n_items}{Integer. Total number of items in the study.}
  \item{n_versions}{Integer. Number of survey versions (respondent splits).
    Default 100.}
  \item{n_tasks}{Integer. Number of tasks per respondent. Default 12.}
  \item{n_items_per_task}{Integer. Number of items shown per task. Default 4.}
  \item{forbidden_pairs}{Optional list of integer vectors of length 2, each
    giving a pair of item indices that must not appear together in any task.}
  \item{seed}{Integer. RNG seed. Default 42.}
}
\value{
  A data frame with columns:
  \describe{
    \item{Version}{Survey version identifier (integer).}
    \item{Task}{Task number within each version (integer).}
    \item{Item1, \dots, ItemK}{Item indices shown in each task position.}
  }
}
\seealso{\code{\link{fit_hbmnl}}, \code{\link{fit_simple}},
  \code{\link{assemble_data}}}
\examples{
\dontrun{
  design <- make_design(n_items = 10, n_tasks = 12, n_items_per_task = 4)
  head(design)
}
}
