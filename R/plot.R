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

# ---- plot.maxdiff_hbmnl ---------------------------------------------------

plot.maxdiff_hbmnl <- function(x,
                               style            = c("distribution", "bar"),
                               scale            = c("raw", "probability"),
                               K_alts           = NULL,
                               sort             = TRUE,
                               label_wrap       = 40,
                               n_bins           = 30,
                               include_footnote = TRUE,
                               ...) {
  if (!inherits(x, "maxdiff_hbmnl")) stop("x must be a maxdiff_hbmnl object.")
  style  <- match.arg(style)
  scale  <- match.arg(scale)
  K_alts <- K_alts %||% x$K

  if (style == "distribution")
    .plot_hbmnl_dist(x, scale, K_alts, sort, label_wrap, n_bins, include_footnote)
  else
    .plot_hbmnl_bar(x, scale, K_alts, sort, label_wrap, include_footnote)
}

# Distribution style: point + CI spine + respondent histogram per item.
.plot_hbmnl_dist <- function(x, scale, K_alts, sort, label_wrap, n_bins, include_footnote) {
  anchored <- isTRUE(x$anchored)
  bm       <- x$beta_mean                         # N x M matrix
  if (scale == "probability") bm <- to_probability_scale(bm, K_alts)

  M     <- ncol(bm)
  means <- colMeans(bm)
  q025  <- apply(bm, 2, stats::quantile, probs = 0.025)
  q975  <- apply(bm, 2, stats::quantile, probs = 0.975)

  lab <- stringr::str_wrap(x$item_labels, width = label_wrap)

  # Assign y positions 1..M; highest mean → top (y = M)
  ord_asc        <- if (sort) order(means) else seq_len(M)
  y_pos          <- integer(M)
  y_pos[ord_asc] <- seq_len(M)

  summary_df <- data.frame(
    y    = y_pos,
    mean = means,
    q025 = q025,
    q975 = q975
  )

  # Histogram bins: shared breaks across items so x-positions are comparable
  x_pad  <- diff(range(bm, na.rm = TRUE)) * 0.02
  breaks <- seq(min(bm, na.rm = TRUE) - x_pad,
                max(bm, na.rm = TRUE) + x_pad,
                length.out = n_bins + 1L)

  max_bar_ht <- 0.40   # histogram bars extend upward from item's y baseline

  hist_rows <- vector("list", M)
  for (j in seq_len(M)) {
    bin_idx <- cut(bm[, j], breaks = breaks, labels = FALSE, include.lowest = TRUE)
    counts  <- tabulate(bin_idx, nbins = n_bins)
    mc      <- max(counts)
    if (mc > 0L) {
      ht <- counts / mc * max_bar_ht
      hist_rows[[j]] <- data.frame(
        xmin = breaks[-length(breaks)],
        xmax = breaks[-1L],
        ymin = y_pos[j],
        ymax = y_pos[j] + ht,
        show = counts > 0L
      )
    }
  }
  hist_df <- dplyr::bind_rows(hist_rows)
  hist_df <- hist_df[hist_df$show, ]

  zero_x  <- if (scale == "probability") 1 / K_alts else 0
  x_label <- if (scale == "probability") {
    glue::glue("Probability score (vs {K_alts - 1L} avg alternatives)")
  } else {
    "Partworth utility"
  }

  p <- ggplot2::ggplot() +
    # Reference line
    ggplot2::geom_vline(xintercept = zero_x, linetype = "dashed",
                        color = .md_gray, linewidth = 0.4) +
    # CI spine: 2.5th–97.5th percentile of respondent posterior means
    ggplot2::geom_segment(
      data = summary_df,
      ggplot2::aes(x = .data$q025, xend = .data$q975,
                   y = .data$y,    yend = .data$y),
      color = .md_dark, linewidth = 1.0, alpha = 0.35
    ) +
    # Respondent histogram (symmetric about item y center)
    ggplot2::geom_rect(
      data = hist_df,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = .md_blue, color = NA, alpha = 0.60
    ) +
    # Mean point (hollow circle on top)
    ggplot2::geom_point(
      data = summary_df,
      ggplot2::aes(x = .data$mean, y = .data$y),
      shape = 21, fill = "white", color = .md_dark, size = 2.5, stroke = 1.3
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(M),
      labels = lab[ord_asc],
      expand = ggplot2::expansion(add = 0.7)
    ) +
    ggplot2::labs(
      x        = x_label,
      y        = NULL,
      title    = "MAXDIFF HB RESULTS",
      subtitle = paste0(
        "Circle: mean across respondents.  ",
        "Spine: 2.5– 97.5th pct range.  ",
        "Histogram: respondent distribution."
      )
    ) +
    .theme_maxdiff(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )

  if (include_footnote)
    p <- p + ggplot2::labs(caption = .interpretation_footnote(scale, anchored, K_alts))
  p
}

# Bar style: simple horizontal bars, colored by anchor status when anchored.
.plot_hbmnl_bar <- function(x, scale, K_alts, sort, label_wrap, include_footnote) {
  anchored <- isTRUE(x$anchored)
  bm       <- x$beta_mean
  if (scale == "probability") bm <- to_probability_scale(bm, K_alts)

  means  <- colMeans(bm)
  zero_x <- if (scale == "probability") 1 / K_alts else 0

  df <- data.frame(
    label = stringr::str_wrap(x$item_labels, width = label_wrap),
    score = means
  )
  if (sort)
    df$label <- factor(df$label, levels = df$label[order(df$score)])

  if (anchored) {
    df$grp     <- ifelse(df$score >= zero_x, "above", "below")
    fill_vals  <- c(above = .md_blue, below = .md_red)
  } else {
    df$grp     <- "all"
    fill_vals  <- c(all = .md_blue)
  }

  x_label <- if (scale == "probability") "Probability score" else "Partworth utility"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$score, y = .data$label,
                                        fill = .data$grp)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_vline(xintercept = zero_x, color = .md_dark, linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
    ggplot2::labs(
      x     = x_label,
      y     = NULL,
      title = "MAXDIFF HB RESULTS"
    ) +
    .theme_maxdiff(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  if (include_footnote)
    p <- p + ggplot2::labs(caption = .interpretation_footnote(scale, anchored, K_alts))
  p
}

# ---- plot.maxdiff_simple --------------------------------------------------

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
      caption = paste(strwrap(
        "Blue bars: fraction of times item was picked BEST out of times shown. Gray bars: fraction picked WORST.",
        width = 100
      ), collapse = "\n")
    )
  p
}

# ---- plot.maxdiff_subgroup ------------------------------------------------

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
        subtitle = "Point: subgroup posterior mean.  Range: 2.5–97.5% across respondents."
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
