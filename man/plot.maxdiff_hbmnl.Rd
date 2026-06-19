\name{plot.maxdiff_hbmnl}
\alias{plot.maxdiff_hbmnl}
\alias{plot.maxdiff_simple}
\alias{plot.maxdiff_subgroup}
\title{Plot MaxDiff results}
\description{
  \itemize{
    \item \code{plot.maxdiff_hbmnl} — dot-and-range plot of posterior-mean
      partworths with 2.5\%–97.5\% cross-respondent heterogeneity bands.
    \item \code{plot.maxdiff_simple} — horizontal bar chart of net B-W scores.
    \item \code{plot.maxdiff_subgroup} — grouped dot plot comparing subgroups.
  }
}
\usage{
\method{plot}{maxdiff_hbmnl}(x, scale = c("raw", "probability"),
     K_alts = NULL, title = NULL, ...)

\method{plot}{maxdiff_simple}(x, title = NULL, ...)

\method{plot}{maxdiff_subgroup}(x, title = NULL, ...)
}
\arguments{
  \item{x}{A \code{maxdiff_hbmnl}, \code{maxdiff_simple}, or
    \code{maxdiff_subgroup} object.}
  \item{scale}{\code{"raw"} (default) or \code{"probability"} (Sawtooth
    probability scale).}
  \item{K_alts}{Integer. Number of alternatives per task (for probability
    scale). Defaults to \code{x$K}.}
  \item{title}{Optional plot title string.}
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
  plot(fit)
  plot(fit, scale = "probability")
  plot(simple_fit)
  plot(subgroup_fit)
}
}
