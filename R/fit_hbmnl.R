fit_hbmnl <- function(responses, design,
                      version_col, best_cols, worst_cols,
                      anchor_cols       = NULL,
                      anchor_codes      = c(all = 1, some = 2, none = 99),
                      response_encoding = c("item", "position"),
                      item_labels       = NULL,
                      R                 = 20000L,
                      burnin            = 10000L,
                      keep              = 10L,
                      nprint            = 1000L,
                      ncomp             = 1L,
                      seed              = 42L) {

  response_encoding <- match.arg(response_encoding)

  if (burnin %% keep != 0L)
    stop(glue::glue("burnin ({burnin}) must be a multiple of keep ({keep})."))
  if (burnin >= R)
    stop(glue::glue("burnin ({burnin}) must be less than R ({R})."))

  bd <- assemble_data(
    responses, design,
    version_col       = version_col,
    best_cols         = best_cols,
    worst_cols        = worst_cols,
    anchor_cols       = anchor_cols,
    anchor_codes      = anchor_codes,
    response_encoding = response_encoding
  )

  M        <- bd$M
  anchored <- bd$anchored
  labels   <- coerce_item_labels(item_labels, M)

  Data  <- list(lgtdata = bd$lgtdata)
  Prior <- list(ncomp = ncomp)
  Mcmc  <- list(R = R, keep = keep, nprint = nprint)

  message(glue::glue(
    "Fitting HB MNL: N={bd$N}, T={bd$T}, K={bd$K}, M={bd$M}, ",
    "anchored={anchored}, R={R}, keep={keep}, burnin={burnin}, ncomp={ncomp}."
  ))
  flush.console()

  set.seed(seed)
  t0            <- proc.time()[["elapsed"]]
  fit           <- .rhierMnlRwMixture(Data = Data, Prior = Prior, Mcmc = Mcmc)
  runtime_secs  <- proc.time()[["elapsed"]] - t0
  message(glue::glue("  ... sampler completed in {round(runtime_secs, 1)} seconds."))

  n_saved <- dim(fit$betadraw)[3]
  n_drop  <- burnin %/% keep
  if (n_drop >= n_saved)
    stop(glue::glue("burnin discards all saved draws (n_drop={n_drop}, n_saved={n_saved})."))
  betadraw_raw <- fit$betadraw[, , (n_drop + 1L):n_saved, drop = FALSE]
  n_kept       <- dim(betadraw_raw)[3]

  if (anchored) {
    betadraw <- betadraw_raw
  } else {
    betadraw <- array(0, dim = c(bd$N, M, n_kept))
    betadraw[, 1:(M - 1L), ] <- betadraw_raw
    for (d in seq_len(n_kept)) {
      slab          <- betadraw[, , d]
      betadraw[, , d] <- slab - rowMeans(slab)
    }
  }

  beta_mean <- apply(betadraw, c(1, 2), mean)

  Deltadraw_post <- if (!is.null(fit$Deltadraw)) {
    fit$Deltadraw[(n_drop + 1L):n_saved, , drop = FALSE]
  } else NULL

  out <- list(
    betadraw         = betadraw,
    beta_mean        = beta_mean,
    item_labels      = unname(labels[as.character(seq_len(M))]),
    Deltadraw        = Deltadraw_post,
    nmix             = fit$nmix,
    n_kept           = n_kept,
    iter_R           = R,
    iter_burnin      = burnin,
    iter_keep        = keep,
    runtime_secs     = runtime_secs,
    respondent_index = bd$respondent_index,
    anchored         = anchored,
    N                = bd$N,
    T                = bd$T,
    K                = bd$K,
    M                = bd$M,
    call             = match.call()
  )
  class(out) <- "maxdiff_hbmnl"
  out
}

print.maxdiff_hbmnl <- function(x, ...) {
  cat("<maxdiff_hbmnl>  (HB MNL, Rossi/Allenby/McCulloch Gibbs sampler)\n")
  cat("  Respondents (N) :", x$N, "\n")
  cat("  Items (M)       :", x$M, "\n")
  cat("  Tasks (T)       :", x$T, "\n")
  cat("  Alts/task (K)   :", x$K, "\n")
  cat("  Anchored        :", x$anchored, "\n")
  cat(glue::glue("  Iterations      : R={x$iter_R}, burnin={x$iter_burnin}, keep={x$iter_keep}"),
      "\n", sep = "")
  cat("  Post-burnin draws:", x$n_kept, "\n")
  cat("  Runtime         :", round(x$runtime_secs, 1), "seconds\n")
  invisible(x)
}

summary.maxdiff_hbmnl <- function(object,
                                  scale  = c("raw", "probability"),
                                  K_alts = NULL,
                                  ...) {
  scale  <- match.arg(scale)
  K_alts <- K_alts %||% object$K
  bm     <- object$beta_mean
  if (scale == "probability") bm <- to_probability_scale(bm, K_alts)
  tibble::tibble(
    item       = seq_len(ncol(bm)),
    item_label = object$item_labels,
    mean       = colMeans(bm),
    sd         = apply(bm, 2, stats::sd),
    q025       = apply(bm, 2, stats::quantile, probs = 0.025),
    q975       = apply(bm, 2, stats::quantile, probs = 0.975)
  )
}

coef.maxdiff_hbmnl <- function(object,
                               scale  = c("raw", "probability"),
                               K_alts = NULL,
                               ...) {
  scale  <- match.arg(scale)
  K_alts <- K_alts %||% object$K
  if (scale == "probability") to_probability_scale(object$beta_mean, K_alts)
  else                        object$beta_mean
}
