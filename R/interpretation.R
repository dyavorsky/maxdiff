# Score interpretation footnotes for plots and tables. Not exported.

.interpretation_footnote <- function(scale = c("raw", "probability"),
                                     anchored = FALSE,
                                     K_alts = NULL) {
  scale <- match.arg(scale)
  if (scale == "raw") {
    if (anchored) {
      paste0(
        "Raw partworths on the anchored scale. The anchor sits at 0; ",
        "positive values are 'above the buy line', negative values 'below'. ",
        "A 1.0 gap means one feature is ~2.7x more likely to be chosen than the other (exp(1) ≈ 2.7)."
      )
    } else {
      paste0(
        "Raw partworths, sum-to-zero centered across items per respondent. ",
        "Positive values are above the per-respondent average; negative values below. ",
        "A 1.0 gap means one feature is ~2.7x more likely to be chosen than the other (exp(1) ≈ 2.7)."
      )
    }
  } else {
    K <- if (is.null(K_alts)) "K" else as.character(as.integer(K_alts))
    if (anchored) {
      paste0(
        "Probability-scaled scores. Each value is the probability of being chosen as MOST APPEALING ",
        "in a ", K, "-item comparison against alternatives at the buy-line anchor (utility 0). ",
        "The anchor's own score is 1/", K, "; items above the buy line score higher."
      )
    } else {
      paste0(
        "Probability-scaled scores. Each value is the probability of being chosen as MOST APPEALING ",
        "in a ", K, "-item comparison against alternatives at the per-respondent average (utility 0). ",
        "Items at the average score 1/", K, "."
      )
    }
  }
}
