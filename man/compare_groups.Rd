\name{compare_groups}
\alias{compare_groups}
\title{Per-item significance tests between two subgroups}
\description{
  Tests whether each item's partworth differs significantly between two
  respondent subgroups. Two methods are available:

  \itemize{
    \item \strong{Welch} (\code{method = "welch"}): two-sample Welch t-test
      on respondent-level posterior means. Reports t, df, p-value, Cohen's d,
      and a significance flag.
    \item \strong{Bayesian} (\code{method = "bayesian"}): per-draw posterior
      comparison. Reports the credible interval of the group-mean difference
      and the probability that group A exceeds group B.
  }
}
\usage{
compare_groups(
  fit_hbmnl,
  group_a,
  group_b,
  method = c("welch", "bayesian"),
  level  = 0.95,
  scale  = c("raw", "probability"),
  K_alts = NULL
)
}
\arguments{
  \item{fit_hbmnl}{A \code{maxdiff_hbmnl} object from \code{\link{fit_hbmnl}}.}
  \item{group_a}{Logical vector identifying group A respondents.}
  \item{group_b}{Logical vector identifying group B respondents.}
  \item{method}{\code{"welch"} (default) or \code{"bayesian"}.}
  \item{level}{Confidence / credible interval coverage. Default 0.95.}
  \item{scale}{\code{"raw"} (default) or \code{"probability"}.}
  \item{K_alts}{Integer. Defaults to \code{fit_hbmnl$K}.}
}
\value{
  A tibble with one row per item. Both methods include columns \code{item},
  \code{item_label}, \code{n_a}, \code{n_b}, \code{mean_a}, \code{mean_b},
  \code{mean_diff}, \code{ci_lo}, \code{ci_hi}, and
  \code{sig_\{level*100\}}.

  The Welch method additionally provides \code{t}, \code{df}, \code{p}, and
  \code{cohens_d}. The Bayesian method additionally provides
  \code{prob_a_gt_b}.
}
\seealso{
  \code{\link{subgroup_summary}}, \code{\link{fit_hbmnl}}
}
\examples{
\dontrun{
  cmp <- compare_groups(
    fit_hbmnl = fit,
    group_a   = responses$age >= 35,
    group_b   = responses$age <  35
  )
  cmp[cmp$sig_95, ]
}
}
