# Neutral color palette for maxdiff plots.
.md_dark  <- "#2C4A6E"
.md_blue  <- "#5B89B7"
.md_gray  <- "#909090"
.md_lgray <- "#D8D8D8"
.md_red   <- "#C0392B"

.md_palette <- c(.md_blue, .md_red, .md_dark, .md_gray, .md_lgray)

.theme_maxdiff <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(color = .md_dark, face = "bold",
                                            size  = base_size * 1.3, hjust = 0),
      plot.subtitle = ggplot2::element_text(color = .md_blue, size = base_size * 1.0, hjust = 0),
      axis.text     = ggplot2::element_text(color = .md_dark, size = base_size * 0.85),
      axis.title    = ggplot2::element_text(color = .md_dark, size = base_size * 0.9),
      panel.grid.major = ggplot2::element_line(color = .md_lgray, linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title  = ggplot2::element_text(face = "bold"),
      strip.text    = ggplot2::element_text(face = "bold", color = "white"),
      strip.background = ggplot2::element_rect(fill = .md_dark)
    )
}

plot.maxdiff_hbmnl <- function(x, fit_simple = NULL,
                               sort             = TRUE,
                               label_wrap       = 40,
                               zero_line        = TRUE,
                               scale            = c("raw", "probability"),
                               K_alts           = NULL,
                               include_footnote = TRUE,
                               ...) {
  if (!inherits(x, "maxdiff_hbmnl")) stop("x must be a maxdiff_hbmnl object.")
  scale    <- match.arg(scale)
  K_alts   <- K_alts %||% x$K
  anchored <- isTRUE(x$anchored)

  s <- summary(x, scale = scale, K_alts = K_alts)
  s$item_label_wrapped <- stringr::str_wrap(s$item_label, width = label_wrap)
  if (sort) {
    s$item_label_wrapped <- factor(s$item_label_wrapped,
                                   levels = s$item_label_wrapped[order(s$mean)])
  } else {
    s$item_label_wrapped <- factor(s$item_label_wrapped,
                                   levels = rev(s$item_label_wrapped))
  }

  x_label <- if (scale == "probability") {
    glue::glue("Probability of being chosen as MOST APPEALING (vs {K_alts - 1} avg alternatives)")
  } else {
    "Partworth utility (posterior mean across respondents)"
  }
  zero_x <- if (scale == "probability") 1 / K_alts else 0

  p <- ggplot2::ggplot(s, ggplot2::aes(x = .data$mean, y = .data$item_label_wrapped)) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = .data$q025, xmax = .data$q975),
                             color = .md_blue, linewidth = 0.4, size = 1.2) +
    ggplot2::labs(
      x        = x_label,
      y        = NULL,
      title    = "MAXDIFF HB RESULTS",
      subtitle = "Point: posterior mean.  Range: 2.5%–97.5% across respondent posterior means."
    ) +
    .theme_maxdiff(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  if (zero_line)
    p <- p + ggplot2::geom_vline(xintercept = zero_x, linetype = "dashed",
                                 color = .md_dark, linewidth = 0.4)

  if (!is.null(fit_simple)) {
    simple_df <- data.frame(
      item_label_wrapped = factor(stringr::str_wrap(fit_simple$item_label, width = label_wrap),
                                  levels = levels(s$item_label_wrapped)),
      simple = fit_simple$score
    )
    rng_pw   <- range(c(s$q025, s$q975), na.rm = TRUE)
    rng_pw   <- rng_pw + c(-1, 1) * 0.05 * diff(rng_pw)
    rng_simp <- c(-1, 1)
    simple_df$simple_scaled <- (simple_df$simple - rng_simp[1]) /
      diff(rng_simp) * diff(rng_pw) + rng_pw[1]
    p <- p +
      ggplot2::geom_point(data = simple_df,
                          ggplot2::aes(x = .data$simple_scaled, y = .data$item_label_wrapped),
                          shape = 4, color = .md_red, size = 2.2, inherit.aes = FALSE)
  }

  if (include_footnote)
    p <- p + ggplot2::labs(caption = .interpretation_footnote(scale, anchored, K_alts))
  p
}

plot.maxdiff_simple <- function(x, label_wrap = 40, include_footnote = TRUE, ...) {
  d           <- x
  d$pct_best  <-  d$n_best  / d$n_shown
  d$pct_worst <- -d$n_worst / d$n_shown
  d$item_label_wrapped <- stringr::str_wrap(d$item_label, width = label_wrap)
  d$item_label_wrapped <- factor(d$item_label_wrapped,
                                 levels = d$item_label_wrapped[order(d$score)])

  p <- ggplot2::ggplot(d, ggplot2::aes(y = .data$item_label_wrapped)) +
    ggplot2::geom_col(ggplot2::aes(x = .data$pct_worst), fill = .md_gray, width = 0.65) +
    ggplot2::geom_col(ggplot2::aes(x = .data$pct_best),  fill = .md_blue, width = 0.65) +
    ggplot2::geom_vline(xintercept = 0, color = .md_dark, linewidth = 0.4) +
    ggplot2::scale_x_continuous(labels = abs) +
    ggplot2::labs(x = NULL, y = NULL, title = "PERCENT BEST VS PERCENT WORST") +
    .theme_maxdiff(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  if (include_footnote)
    p <- p + ggplot2::labs(
      caption = "Blue bars: fraction of times item was picked BEST out of times shown. Gray bars: fraction picked WORST."
    )
  p
}

plot.maxdiff_subgroup <- function(x, metric = c("hbmnl", "simple"),
                                  sort_by = NULL, label_wrap = 40,
                                  include_footnote = TRUE, ...) {
  metric   <- match.arg(metric)
  d        <- x
  scale    <- attr(x, "scale")    %||% "raw"
  K_alts   <- attr(x, "K_alts")   %||% NA_integer_
  anchored <- attr(x, "anchored") %||% FALSE

  if (!is.factor(d$group)) d$group <- factor(d$group, levels = unique(d$group))
  d$item_label_wrapped <- stringr::str_wrap(d$item_label, width = label_wrap)

  if (!is.null(sort_by)) {
    if (!sort_by %in% levels(d$group))
      stop(glue::glue("sort_by='{sort_by}' not found among groups."))
    sort_vals <- d[d$group == sort_by, , drop = FALSE]
  } else {
    sort_vals <- dplyr::summarise(
      dplyr::group_by(d, .data$item, .data$item_label_wrapped),
      sort_metric = if (metric == "hbmnl") mean(.data$hb_mean, na.rm = TRUE)
                    else                   mean(.data$simple_score, na.rm = TRUE),
      .groups = "drop"
    )
    sort_vals$hb_mean      <- sort_vals$sort_metric
    sort_vals$simple_score <- sort_vals$sort_metric
  }
  metric_col <- if (metric == "hbmnl") "hb_mean" else "simple_score"
  ord        <- sort_vals$item[order(sort_vals[[metric_col]])]
  level_lookup <- setNames(d$item_label_wrapped[match(ord, d$item)], ord)
  d$item_label_wrapped <- factor(d$item_label_wrapped, levels = unname(level_lookup))

  if (metric == "hbmnl") {
    x_label <- if (scale == "probability") {
      glue::glue("Probability score (vs {K_alts - 1} avg alternatives)")
    } else {
      "HB partworth (posterior mean across respondents in group)"
    }
    zero_x <- if (scale == "probability") 1 / K_alts else 0
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$hb_mean,
                                         y = .data$item_label_wrapped,
                                         color = .data$group)) +
      ggplot2::geom_pointrange(
        ggplot2::aes(xmin = .data$hb_q025, xmax = .data$hb_q975),
        position  = ggplot2::position_dodge(width = 0.55),
        linewidth = 0.3, size = 0.75
      ) +
      ggplot2::labs(
        x = x_label, y = NULL, color = "Subgroup",
        title    = "MAXDIFF SUBGROUP COMPARISON (HB MNL)",
        subtitle = "Point: subgroup posterior mean.  Range: 2.5%–97.5% across respondents."
      )
  } else {
    zero_x <- 0
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$simple_score,
                                         y = .data$item_label_wrapped,
                                         color = .data$group)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55), size = 2.4) +
      ggplot2::labs(
        x = "Simple best-minus-worst score (-1 to 1)", y = NULL, color = "Subgroup",
        title = "MAXDIFF SUBGROUP COMPARISON (SIMPLE B-W)"
      )
  }

  p <- p +
    ggplot2::geom_vline(xintercept = zero_x, linetype = "dashed",
                        color = .md_dark, linewidth = 0.4) +
    ggplot2::scale_color_manual(values = .md_palette) +
    .theme_maxdiff(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  if (include_footnote && metric == "hbmnl")
    p <- p + ggplot2::labs(caption = .interpretation_footnote(scale, anchored, K_alts))
  p
}
