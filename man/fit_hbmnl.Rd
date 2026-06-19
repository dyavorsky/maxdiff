\name{fit_hbmnl}
\alias{fit_hbmnl}
\alias{print.maxdiff_hbmnl}
\alias{summary.maxdiff_hbmnl}
\alias{coef.maxdiff_hbmnl}
\title{Fit a hierarchical Bayesian MNL to MaxDiff data}
\description{
  Runs the HB MNL Gibbs sampler (Rossi/Allenby/McCulloch algorithm) on
  MaxDiff best-then-worst responses. Each task is encoded as the exact joint
  likelihood: K alternatives for the best pick, K-1 for the worst pick.

  The sampler is self-contained (no external \pkg{bayesm} dependency) and was
  adapted from Peter Rossi's \pkg{bayesm} package with an extension to support
  ragged choice sets. Reference: Rossi, Allenby, and McCulloch (2024),
  \emph{Bayesian Statistics and Marketing} (2nd ed.), Wiley.
}
\usage{
fit_hbmnl(
  responses,
  design,
  version_col,
  best_cols,
  worst_cols,
  anchor_cols       = NULL,
  anchor_codes      = c(all = 1, some = 2, none = 99),
  response_encoding = c("item", "position"),
  item_labels       = NULL,
  R                 = 20000L,
  burnin            = 10000L,
  keep              = 10L,
  nprint            = 1000L,
  ncomp             = 1L,
  seed              = 42L
)
}
\arguments{
  \item{responses}{Data frame with one row per respondent (wide format).}
  \item{design}{Data frame with columns Version, Task (or Set), and item
    columns (Item1 … ItemK or similar).}
  \item{version_col}{Character. Name of the version column in \code{responses}.}
  \item{best_cols}{Named character vector. Names = task number (as string),
    values = column names in \code{responses} containing the best pick.}
  \item{worst_cols}{Same shape as \code{best_cols}, for worst picks.}
  \item{anchor_cols}{Optional named character vector of anchor-task columns.}
  \item{anchor_codes}{Named integer vector mapping data codes to
    \code{"all"}, \code{"some"}, \code{"none"}.
    Default \code{c(all = 1, some = 2, none = 99)}.}
  \item{response_encoding}{\code{"item"} (default) if picks are item IDs;
    \code{"position"} if picks are 1-based position indices in the task.}
  \item{item_labels}{Optional length-M character vector of item names.}
  \item{R}{Integer. Total Gibbs iterations. Default 20000.}
  \item{burnin}{Integer. Iterations to discard; must be a multiple of
    \code{keep}. Default 10000.}
  \item{keep}{Integer. Thinning interval; every \code{keep}-th draw is
    saved. Default 10.}
  \item{nprint}{Integer. Progress-print interval (0 = silent). Default 1000.}
  \item{ncomp}{Integer. Mixture components in the upper-level prior.
    Default 1 (standard HB MNL with a single MVN prior).}
  \item{seed}{Integer. RNG seed. Default 42.}
}
\value{
  An S3 object of class \code{"maxdiff_hbmnl"} with elements:
  \describe{
    \item{betadraw}{N × M × n_kept array of respondent partworth draws.}
    \item{beta_mean}{N × M matrix of respondent posterior means.}
    \item{item_labels}{Length-M character vector of item names.}
    \item{n_kept, iter_R, iter_burnin, iter_keep}{MCMC metadata.}
    \item{runtime_secs}{Wall-clock seconds for the sampler.}
    \item{respondent_index}{Row indices in \code{responses} of included
      respondents.}
    \item{anchored, N, T, K, M}{Study dimensions and anchor flag.}
  }
}
\seealso{
  \code{\link{fit_simple}}, \code{\link{make_design}},
  \code{\link{assemble_data}}, \code{\link{subgroup_summary}},
  \code{\link{compare_groups}}, \code{\link{plot.maxdiff_hbmnl}}
}
\examples{
\dontrun{
  design <- make_design(n_items = 10, n_tasks = 12, n_items_per_task = 4)
  # responses: data frame with version, best1..best12, worst1..worst12 columns
  fit <- fit_hbmnl(
    responses   = responses,
    design      = design,
    version_col = "version",
    best_cols   = setNames(paste0("best",  1:12), as.character(1:12)),
    worst_cols  = setNames(paste0("worst", 1:12), as.character(1:12)),
    item_labels = c("Price", "Quality", "Speed", "Support",
                    "Design", "Battery", "Weight", "Screen",
                    "Camera", "Storage")
  )
  print(fit)
  summary(fit)
  coef(fit)
  plot(fit)
}
}
