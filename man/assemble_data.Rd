\name{assemble_data}
\alias{assemble_data}
\title{Assemble lgtdata for the HB MNL sampler}
\description{
  Converts survey responses and design into the \code{lgtdata} list format
  consumed by the internal Gibbs sampler. Each respondent becomes one element
  of \code{lgtdata} containing a \code{y} vector and a ragged \code{X} list.

  Best-then-worst encoding: each MaxDiff task produces two MNL observations.
  The best-pick observation uses dummy-coded rows with coefficient \code{+1};
  the worst-pick observation uses rows with coefficient \code{-1} and has
  K-1 alternatives (removing the item that was chosen best).

  You typically call this indirectly via \code{\link{fit_hbmnl}}. Use it
  directly only if you need access to the raw \code{lgtdata} structure.
}
\usage{
assemble_data(
  responses,
  design,
  version_col,
  best_cols,
  worst_cols,
  anchor_cols       = NULL,
  anchor_codes      = c(all = 1, some = 2, none = 99),
  response_encoding = c("item", "position")
)
}
\arguments{
  \item{responses}{Data frame with one row per respondent (wide format).}
  \item{design}{Data frame with columns Version, Task (or Set), and item
    columns.}
  \item{version_col}{Character. Name of the version column in
    \code{responses}.}
  \item{best_cols}{Named character vector: names = task number (as string),
    values = column names containing the best pick.}
  \item{worst_cols}{Same shape as \code{best_cols}, for worst picks.}
  \item{anchor_cols}{Optional named character vector of anchor-task columns.}
  \item{anchor_codes}{Named integer vector mapping data codes to
    \code{"all"}, \code{"some"}, \code{"none"}.
    Default \code{c(all = 1, some = 2, none = 99)}.}
  \item{response_encoding}{\code{"item"} (default) or \code{"position"}.}
}
\value{
  A list with elements:
  \describe{
    \item{lgtdata}{List of length N. Each element has \code{y} (integer
      vector of 2T choice indicators), \code{X} (stacked design matrix with
      ragged per-occasion sizes), and \code{hess} (pre-computed Hessian
      used to scale the Metropolis proposal).}
    \item{p_it_list}{List of length N. Each element is an integer vector of
      per-occasion choice-set sizes for the ragged indexing.}
    \item{N, T, K, M, anchored, respondent_index}{Study metadata.}
  }
}
\seealso{\code{\link{fit_hbmnl}}}
