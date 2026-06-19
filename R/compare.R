compare_groups <- function(fit_hbmnl,
                           group_a, group_b,
                           method  = c("welch", "bayesian"),
                           level   = 0.95,
                           scale   = c("raw", "probability"),
                           K_alts  = NULL) {

  if (!inherits(fit_hbmnl, "maxdiff_hbmnl"))
    stop("fit_hbmnl must be a maxdiff_hbmnl object.")
  method <- match.arg(method)
  scale  <- match.arg(scale)
  K_alts <- K_alts %||% fit_hbmnl$K
  if (!is.numeric(level) || level <= 0 || level >= 1)
    stop("level must be a number in (0, 1).")

  mask_a <- .to_fit_mask(group_a, fit_hbmnl)
  mask_b <- .to_fit_mask(group_b, fit_hbmnl)
  if (!any(mask_a) || !any(mask_b))
    stop("Both group_a and group_b must have at least one respondent in the fitted set.")
  if (any(mask_a & mask_b))
    warning(glue::glue("{sum(mask_a & mask_b)} respondent(s) are in both groups; results may be biased."))

  M      <- fit_hbmnl$M
  labels <- fit_hbmnl$item_labels
  alpha  <- 1 - level
  q_lo   <- alpha / 2
  q_hi   <- 1 - alpha / 2

  if (method == "welch") {
    bm <- fit_hbmnl$beta_mean
    if (scale == "probability") bm <- to_probability_scale(bm, K_alts)
    rows <- lapply(seq_len(M), function(m) .welch_one_diff(bm[mask_a, m], bm[mask_b, m], level))
    out  <- tibble::as_tibble(do.call(rbind, lapply(rows, as.data.frame)))
    out$item       <- seq_len(M)
    out$item_label <- labels
    out <- out[, c("item", "item_label", "n_a", "n_b", "mean_a", "mean_b", "mean_diff",
                   "t", "df", "p", "ci_lo", "ci_hi", "cohens_d", "sig")]
    names(out)[names(out) == "sig"] <- glue::glue("sig_{round(level * 100)}")
    return(out)
  }

  bd <- fit_hbmnl$betadraw
  if (scale == "probability") bd <- to_probability_scale(bd, K_alts)
  D            <- dim(bd)[3]
  diff_draws   <- matrix(NA_real_, nrow = D, ncol = M)
  mean_a_per   <- matrix(NA_real_, nrow = D, ncol = M)
  mean_b_per   <- matrix(NA_real_, nrow = D, ncol = M)
  for (d in seq_len(D)) {
    slab              <- bd[, , d]
    mean_a_per[d, ]   <- colMeans(slab[mask_a, , drop = FALSE])
    mean_b_per[d, ]   <- colMeans(slab[mask_b, , drop = FALSE])
    diff_draws[d, ]   <- mean_a_per[d, ] - mean_b_per[d, ]
  }
  mean_diff <- colMeans(diff_draws)
  ci_lo     <- apply(diff_draws, 2, stats::quantile, probs = q_lo, names = FALSE)
  ci_hi     <- apply(diff_draws, 2, stats::quantile, probs = q_hi, names = FALSE)
  prob_gt   <- colMeans(diff_draws > 0)
  sig_flag  <- (ci_lo > 0) | (ci_hi < 0)
  out <- tibble::tibble(
    item        = seq_len(M),
    item_label  = labels,
    n_a         = sum(mask_a),
    n_b         = sum(mask_b),
    mean_a      = colMeans(mean_a_per),
    mean_b      = colMeans(mean_b_per),
    mean_diff   = mean_diff,
    ci_lo       = ci_lo,
    ci_hi       = ci_hi,
    prob_a_gt_b = prob_gt
  )
  out[[glue::glue("sig_{round(level * 100)}")]] <- sig_flag
  out
}

.to_fit_mask <- function(flag, fit_hbmnl) {
  if (!is.logical(flag)) flag <- as.logical(flag)
  flag[is.na(flag)] <- FALSE
  resp_idx <- fit_hbmnl$respondent_index
  N_fit    <- nrow(fit_hbmnl$beta_mean)
  if (length(flag) == max(resp_idx, 0L) || length(flag) >= max(resp_idx, 0L)) {
    out <- flag[resp_idx]
  } else if (length(flag) == N_fit) {
    out <- flag
  } else {
    stop(glue::glue("Subgroup vector length ({length(flag)}) doesn't match responses or fitted count ({N_fit})."))
  }
  out[is.na(out)] <- FALSE
  out
}

.welch_one_diff <- function(x_a, x_b, level = 0.95) {
  m_a <- mean(x_a); m_b <- mean(x_b)
  n_a <- length(x_a); n_b <- length(x_b)
  v_a <- stats::var(x_a); v_b <- stats::var(x_b)
  se   <- sqrt(v_a / n_a + v_b / n_b)
  diff <- m_a - m_b
  t    <- diff / se
  df   <- (v_a / n_a + v_b / n_b)^2 /
          ((v_a / n_a)^2 / (n_a - 1) + (v_b / n_b)^2 / (n_b - 1))
  p    <- 2 * stats::pt(-abs(t), df)
  pooled_sd <- sqrt(((n_a - 1) * v_a + (n_b - 1) * v_b) / (n_a + n_b - 2))
  d    <- diff / pooled_sd
  q    <- stats::qt(1 - (1 - level) / 2, df)
  list(n_a = n_a, n_b = n_b, mean_a = m_a, mean_b = m_b, mean_diff = diff,
       t = t, df = df, p = p, ci_lo = diff - q * se, ci_hi = diff + q * se,
       cohens_d = d, sig = p < (1 - level))
}
