\name{to_probability_scale}
\alias{to_probability_scale}
\title{Convert HB partworths to the Sawtooth probability scale}
\description{
  Applies the Sawtooth-compatible probability transform:
  \deqn{p_j = \frac{e^{\beta_j}}{K - 1 + e^{\beta_j}}}
  where K is the number of alternatives per task. On this scale, an average
  item scores 1/K and a hypothetical item with infinite utility scores 1.

  Operates on numeric vectors, matrices, and 3D arrays (as returned by
  \code{fit_hbmnl$betadraw}).
}
\usage{
to_probability_scale(beta, K_alts)
}
\arguments{
  \item{beta}{Numeric vector, matrix (N × M), or 3D array (N × M × draws)
    of HB partworths.}
  \item{K_alts}{Integer. Number of alternatives per task.}
}
\value{
  Numeric object of the same shape as \code{beta}, with values in (0, 1).
}
\examples{
\dontrun{
  # Overall mean on probability scale
  to_probability_scale(colMeans(fit$beta_mean), K_alts = fit$K)

  # Per-respondent means on probability scale
  to_probability_scale(fit$beta_mean, K_alts = fit$K)
}
}
