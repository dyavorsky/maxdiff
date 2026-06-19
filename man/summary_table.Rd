\name{summary_table}
\alias{summary_table}
\title{Publication-ready gt table of MaxDiff results}
\description{
  Produces a \pkg{gt} table from a \code{maxdiff_hbmnl} fit. Two modes:

  \itemize{
    \item \strong{Overall} (\code{subgroup_results = NULL}): one row per item
      with mean, SD, and 2.5\%–97.5\% heterogeneity range.
    \item \strong{By subgroup}: pivots \code{\link{subgroup_summary}} output
      so each subgroup gets its own mean column.
  }

  Requires the \pkg{gt} package. Subgroup mode additionally requires
  \pkg{tidyr}.
}
\usage{
summary_table(
  fit_hbmnl,
  subgroup_results = NULL,
  scale            = c("raw", "probability"),
  K_alts           = NULL,
  include_footnote = TRUE,
  digits           = 2L
)
}
\arguments{
  \item{fit_hbmnl}{A \code{maxdiff_hbmnl} object from \code{\link{fit_hbmnl}}.}
  \item{subgroup_results}{Optional \code{maxdiff_subgroup} from
    \code{\link{subgroup_summary}}.}
  \item{scale}{\code{"raw"} (default) or \code{"probability"}. Ignored when
    \code{subgroup_results} is supplied (reads scale from its attributes).}
  \item{K_alts}{Integer. Defaults to \code{fit_hbmnl$K}.}
  \item{include_footnote}{Logical. Add an interpretive footnote. Default
    \code{TRUE}.}
  \item{digits}{Integer decimal places. Default 2.}
}
\value{
  A \code{gt_tbl} object.
}
\seealso{
  \code{\link{fit_hbmnl}}, \code{\link{subgroup_summary}}
}
\examples{
\dontrun{
  summary_table(fit)
  summary_table(fit, scale = "probability")
  summary_table(fit, subgroup_results = sg)
}
}
