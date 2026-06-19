\name{plot.maxdiff_hbmnl}
\alias{plot.maxdiff_hbmnl}
\alias{plot.maxdiff_simple}
\alias{plot.maxdiff_subgroup}
\title{Plot MaxDiff results}
\description{
  \itemize{
    \item \code{plot.maxdiff_hbmnl} — two styles:
      \describe{
        \item{\code{"distribution"}}{For each item: a symmetric histogram of
          respondent-level posterior means, a 2.5\%–97.5\% spine, and a hollow
          circle marking the mean. Items sorted by mean (default).}
        \item{\code{"bar"}}{Simple horizontal bars of mean scores. For anchored
          studies, bars above the reference are blue and bars below are red.}
      }
    \item \code{plot.maxdiff_simple} — tornado chart: blue bars (pct best) and
      gray bars (pct worst) extending from zero in opposite directions.
    \item \code{plot.maxdiff_subgroup} — dodged point-range chart comparing
      subgroups; one row per item, colored by subgroup.
  }
}
\usage{
\method{plot}{maxdiff_hbmnl}(x, style = c("distribution", "bar"),
     scale = c("raw", "probability"), K_alts = NULL,
     sort = TRUE, label_wrap = 40, n_bins = 30,
     include_footnote = TRUE, ...)

\method{plot}{maxdiff_simple}(x, label_wrap = 40,
     include_footnote = TRUE, ...)

\method{plot}{maxdiff_subgroup}(x, metric = c("hbmnl", "simple"),
     sort_by = NULL, label_wrap = 40, include_footnote = TRUE, ...)
}
\arguments{
  \item{x}{A \code{maxdiff_hbmnl}, \code{maxdiff_simple}, or
    \code{maxdiff_subgroup} object.}
  \item{style}{\code{"distribution"} (default) or \code{"bar"}.}
  \item{scale}{\code{"raw"} (default) or \code{"probability"} (Sawtooth
    probability scale, anchor = 1/K).}
  \item{K_alts}{Integer. Alternatives per task for probability transform;
    defaults to \code{x$K}.}
  \item{sort}{Logical. Sort items by mean score, descending. Default
    \code{TRUE}.}
  \item{label_wrap}{Integer. Character width for wrapping item labels.
    Default 40.}
  \item{n_bins}{Integer. Number of histogram bins for the distribution style.
    Default 30.}
  \item{include_footnote}{Logical. Attach an interpretation caption. Default
    \code{TRUE}.}
  \item{metric}{\code{"hbmnl"} (default) or \code{"simple"} for subgroup
    plots.}
  \item{sort_by}{Name of a subgroup to sort by; defaults to mean across
    groups.}
  \item{...}{Additional arguments (ignored).}
}
\value{
  A \pkg{ggplot2} object.
}
\seealso{
  \code{\link{fit_hbmnl}}, \code{\link{fit_simple}},
  \code{\link{subgroup_summary}}
}
\examples{
\dontrun{
  # Distribution style (default): histogram + spine + mean point
  plot(fit)
  plot(fit, scale = "probability")

  # Bar style: colored by anchor status when anchored = TRUE
  plot(fit, style = "bar")
  plot(fit, style = "bar", scale = "probability")

  plot(simple_fit)
  plot(subgroup_fit)
}
}
