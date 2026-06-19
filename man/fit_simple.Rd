\name{fit_simple}
\alias{fit_simple}
\alias{print.maxdiff_simple}
\title{Simple Best-Worst scoring for MaxDiff data}
\description{
  Tallies best and worst picks per item and computes a simple score
  \code{(best picks - worst picks) / (total tasks * N)} for each item.
  Returns a tibble with one row per item, sorted by score (descending).
}
\usage{
fit_simple(
  responses,
  design,
  version_col,
  best_cols,
  worst_cols,
  response_encoding = c("item", "position"),
  item_labels       = NULL
)
}
\arguments{
  \item{responses}{Data frame with one row per respondent (wide format).}
  \item{design}{Data frame with columns Version, Task (or Set), and item
    columns (Item1 … ItemK or similar).}
  \item{version_col}{Character. Name of the version column in
    \code{responses}.}
  \item{best_cols}{Named character vector. Names = task number (as string),
    values = column names in \code{responses} containing the best pick.}
  \item{worst_cols}{Same shape as \code{best_cols}, for worst picks.}
  \item{response_encoding}{\code{"item"} (default) if picks are item IDs;
    \code{"position"} if picks are 1-based position indices in the task.}
  \item{item_labels}{Optional length-M character vector of item names.}
}
\value{
  A tibble (also of class \code{"maxdiff_simple"}) with columns:
  \describe{
    \item{item}{Integer item index.}
    \item{item_label}{Item label.}
    \item{best_count}{Number of times chosen as best.}
    \item{worst_count}{Number of times chosen as worst.}
    \item{net_score}{Best minus worst per person per task, as a proportion.}
    \item{rank}{Rank by \code{net_score} (1 = best).}
  }
}
\seealso{\code{\link{fit_hbmnl}}, \code{\link{plot.maxdiff_simple}}}
\examples{
\dontrun{
  design <- make_design(n_items = 8, n_tasks = 8, n_items_per_task = 4)
  fit    <- fit_simple(
    responses   = responses,
    design      = design,
    version_col = "version",
    best_cols   = setNames(paste0("best",  1:8), as.character(1:8)),
    worst_cols  = setNames(paste0("worst", 1:8), as.character(1:8))
  )
  print(fit)
  plot(fit)
}
}
