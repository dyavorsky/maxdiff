summary_table <- function(fit_hbmnl,
                          subgroup_results = NULL,
                          scale            = c("raw", "probability"),
                          K_alts           = NULL,
                          include_footnote = TRUE,
                          digits           = 2L) {

  if (!inherits(fit_hbmnl, "maxdiff_hbmnl"))
    stop("fit_hbmnl must be a maxdiff_hbmnl object.")
  if (!requireNamespace("gt", quietly = TRUE))
    stop("Install the 'gt' package to use summary_table().")

  anchored <- isTRUE(fit_hbmnl$anchored)

  if (is.null(subgroup_results)) {
    scale  <- match.arg(scale)
    K_alts <- K_alts %||% fit_hbmnl$K
    s      <- summary(fit_hbmnl, scale = scale, K_alts = K_alts)

    tbl <- gt::gt(s) |>
      gt::cols_label(item = "Item ID", item_label = "Item",
                     mean = "Mean", sd = "SD", q025 = "2.5%", q975 = "97.5%") |>
      gt::fmt_number(columns = c("mean", "sd", "q025", "q975"), decimals = digits) |>
      gt::tab_header(
        title    = "MaxDiff partworth summary",
        subtitle = if (scale == "probability") "Probability-scaled scores" else "Raw partworths"
      ) |>
      gt::tab_spanner(label = "Heterogeneity range", columns = c("q025", "q975"))
  } else {
    sa_scale  <- attr(subgroup_results, "scale")    %||% "raw"
    sa_K_alts <- attr(subgroup_results, "K_alts")   %||% fit_hbmnl$K
    sa_anch   <- attr(subgroup_results, "anchored") %||% anchored
    scale     <- sa_scale; K_alts <- sa_K_alts; anchored <- sa_anch

    if (!requireNamespace("tidyr", quietly = TRUE))
      stop("Install the 'tidyr' package for by-subgroup summary_table().")

    wide <- tidyr::pivot_wider(
      subgroup_results[, c("group", "item", "item_label", "hb_mean")],
      names_from  = "group", values_from = "hb_mean"
    )
    group_cols <- setdiff(names(wide), c("item", "item_label"))

    tbl <- gt::gt(wide) |>
      gt::cols_label(item = "Item ID", item_label = "Item") |>
      gt::fmt_number(columns = group_cols, decimals = digits) |>
      gt::tab_header(
        title    = "MaxDiff partworths by subgroup",
        subtitle = if (scale == "probability") "Probability-scaled scores" else "Raw partworths"
      )
  }

  if (include_footnote)
    tbl <- gt::tab_footnote(tbl,
                            footnote  = .interpretation_footnote(scale, anchored, K_alts),
                            locations = gt::cells_title(groups = "subtitle"))
  tbl
}
