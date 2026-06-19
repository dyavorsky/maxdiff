subgroup_summary <- function(fit_hbmnl, subgroups,
                             fit_simple  = NULL,
                             responses   = NULL, design      = NULL,
                             version_col = NULL, best_cols   = NULL,
                             worst_cols  = NULL, item_labels = NULL,
                             scale       = c("raw", "probability"),
                             K_alts      = NULL) {

  if (!inherits(fit_hbmnl, "maxdiff_hbmnl"))
    stop("fit_hbmnl must be a maxdiff_hbmnl object.")
  if (!is.list(subgroups) || is.null(names(subgroups)) || any(!nzchar(names(subgroups))))
    stop("subgroups must be a fully named list of logical vectors.")

  scale  <- match.arg(scale)
  K_alts <- K_alts %||% fit_hbmnl$K

  beta_mean <- fit_hbmnl$beta_mean
  if (scale == "probability") beta_mean <- to_probability_scale(beta_mean, K_alts)
  M        <- ncol(beta_mean)
  resp_idx <- fit_hbmnl$respondent_index

  labels <- if (!is.null(item_labels)) coerce_item_labels(item_labels, M)
            else setNames(fit_hbmnl$item_labels, seq_len(M))

  rows <- list()
  for (gname in names(subgroups)) {
    flag <- subgroups[[gname]]
    if (!is.logical(flag)) flag <- as.logical(flag)
    flag[is.na(flag)] <- FALSE

    if (length(flag) == max(resp_idx, 0L) || length(flag) >= max(resp_idx, 0L)) {
      fit_mask <- flag[resp_idx]
    } else if (length(flag) == nrow(beta_mean)) {
      fit_mask <- flag
    } else {
      stop(glue::glue(
        "Subgroup '{gname}' length ({length(flag)}) doesn't match responses or fitted count ({nrow(beta_mean)})."
      ))
    }
    fit_mask[is.na(fit_mask)] <- FALSE
    if (!any(fit_mask)) {
      warning(glue::glue("Subgroup '{gname}' has 0 respondents in the fitted set; skipping."))
      next
    }
    bm      <- beta_mean[fit_mask, , drop = FALSE]
    hb_mean <- colMeans(bm)
    hb_sd   <- apply(bm, 2, stats::sd)
    hb_q025 <- apply(bm, 2, stats::quantile, probs = 0.025, names = FALSE)
    hb_q975 <- apply(bm, 2, stats::quantile, probs = 0.975, names = FALSE)

    simple_score <- rep(NA_real_,    M)
    n_shown      <- rep(NA_integer_, M)
    if (!is.null(fit_simple) && !is.null(responses) && !is.null(design) &&
        !is.null(version_col) && !is.null(best_cols) && !is.null(worst_cols)) {
      sub_resp <- responses[flag, , drop = FALSE]
      if (nrow(sub_resp) > 0L) {
        sb <- fit_simple(
          responses   = sub_resp, design = design,
          version_col = version_col, best_cols = best_cols, worst_cols = worst_cols,
          item_labels = item_labels %||% fit_hbmnl$item_labels
        )
        simple_score[sb$item] <- sb$score
        n_shown[sb$item]      <- sb$n_shown
      }
    } else if (!is.null(fit_simple)) {
      simple_score[fit_simple$item] <- fit_simple$score
      n_shown[fit_simple$item]      <- fit_simple$n_shown
    }

    rows[[gname]] <- tibble::tibble(
      group        = gname,
      item         = seq_len(M),
      item_label   = unname(labels[as.character(seq_len(M))]),
      n_in_group   = sum(fit_mask),
      simple_score = simple_score,
      n_shown      = n_shown,
      hb_mean      = hb_mean,
      hb_sd        = hb_sd,
      hb_q025      = hb_q025,
      hb_q975      = hb_q975
    )
  }

  out        <- dplyr::bind_rows(rows)
  out$group  <- factor(out$group, levels = names(subgroups))
  class(out) <- c("maxdiff_subgroup", class(out))
  attr(out, "scale")    <- scale
  attr(out, "K_alts")   <- K_alts
  attr(out, "anchored") <- fit_hbmnl$anchored
  out
}

print.maxdiff_subgroup <- function(x, ...) {
  cat("<maxdiff_subgroup>  scale:", attr(x, "scale"), "\n")
  cat("  Groups:", paste(levels(x$group), collapse = ", "), "\n\n")
  print(tibble::as_tibble(x), ...)
  invisible(x)
}
