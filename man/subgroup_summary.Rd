\name{subgroup_summary}
\alias{subgroup_summary}
\alias{print.maxdiff_subgroup}
\title{Summarize HB partworths by subgroup}
\description{
  For each subgroup, averages respondent posterior means to produce a
  group-level partworth summary. Returns a long tibble with one row per
  subgroup × item combination, of class \code{"maxdiff_subgroup"}.

  Subgroups are defined by a named list of logical vectors (one per group),
  each of length equal to the number of rows in the original \code{responses}
  data frame or the number of fitted respondents.
}
\usage{
subgroup_summary(
  fit_hbmnl,
  subgroups,
  scale  = c("raw", "probability"),
  K_alts = NULL
)
}
\arguments{
  \item{fit_hbmnl}{A \code{maxdiff_hbmnl} object from \code{\link{fit_hbmnl}}.}
  \item{subgroups}{Named list of logical vectors identifying each subgroup.}
  \item{scale}{\code{"raw"} (default) or \code{"probability"}.}
  \item{K_alts}{Integer. Defaults to \code{fit_hbmnl$K}.}
}
\value{
  A tibble of class \code{"maxdiff_subgroup"} with columns:
  \describe{
    \item{group}{Subgroup name.}
    \item{n}{Respondent count in the subgroup.}
    \item{item}{Item index.}
    \item{item_label}{Item label.}
    \item{hb_mean}{Average partworth for that subgroup × item.}
    \item{hb_sd}{Standard deviation across respondents in the subgroup.}
    \item{rank}{Rank within subgroup (1 = highest).}
  }
  Attributes \code{scale}, \code{K_alts}, and \code{anchored} are preserved
  for downstream use by \code{plot()} and \code{summary_table()}.
}
\seealso{
  \code{\link{fit_hbmnl}}, \code{\link{compare_groups}},
  \code{\link{plot.maxdiff_subgroup}}, \code{\link{summary_table}}
}
\examples{
\dontrun{
  sg <- subgroup_summary(
    fit_hbmnl = fit,
    subgroups = list(
      Male   = responses$gender == "Male",
      Female = responses$gender == "Female"
    )
  )
  print(sg)
  plot(sg)
}
}
