to_probability_scale <- function(beta, K_alts) {
  if (is.null(K_alts) || !is.finite(K_alts) || K_alts < 2L)
    stop("K_alts must be an integer >= 2.")
  exp(beta) / (K_alts - 1 + exp(beta))
}
