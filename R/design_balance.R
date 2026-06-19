# Internal helpers for make_design(). Not exported.

.balance_matrices <- function(blocks, n_items) {
  K    <- ncol(blocks)
  B    <- nrow(blocks)
  freq <- tabulate(as.integer(blocks), nbins = n_items)
  pair <- matrix(0L, nrow = n_items, ncol = n_items)
  for (b in seq_len(B)) {
    items <- as.integer(blocks[b, ])
    for (i in seq_len(K - 1L)) {
      for (j in (i + 1L):K) {
        a <- items[i]; bb <- items[j]
        if (a == bb) next
        pair[a, bb] <- pair[a, bb] + 1L
        pair[bb, a] <- pair[bb, a] + 1L
      }
    }
  }
  pos <- matrix(0L, nrow = n_items, ncol = K)
  for (k in seq_len(K)) pos[, k] <- tabulate(as.integer(blocks[, k]), nbins = n_items)
  list(freq = freq, pair = pair, pos = pos)
}

.block_has_forbidden <- function(items, forbidden_pairs) {
  if (is.null(forbidden_pairs) || !length(forbidden_pairs)) return(FALSE)
  for (fp in forbidden_pairs) if (all(fp %in% items)) return(TRUE)
  FALSE
}

.repair_forbidden <- function(blocks, forbidden_pairs, max_iter = 5000) {
  bad_idx <- function() which(apply(blocks, 1, function(r) .block_has_forbidden(r, forbidden_pairs)))
  for (iter in seq_len(max_iter)) {
    bad <- bad_idx()
    if (!length(bad)) return(blocks)
    i       <- bad[[1]]
    bad_row <- blocks[i, ]
    fp_match <- vapply(forbidden_pairs, function(p) all(p %in% bad_row), logical(1))
    fp       <- forbidden_pairs[[which(fp_match)[1]]]
    out_item <- sample(fp, 1)
    candidates <- sample(setdiff(seq_len(nrow(blocks)), i))
    for (j in candidates) {
      donor_row    <- blocks[j, ]
      if (out_item %in% donor_row) next
      donor_options <- setdiff(donor_row, bad_row)
      for (in_item in sample(donor_options)) {
        new_bad   <- replace(bad_row, which(bad_row == out_item), in_item)
        new_donor <- replace(donor_row, which(donor_row == in_item), out_item)
        if (!.block_has_forbidden(new_bad, forbidden_pairs) &&
            !.block_has_forbidden(new_donor, forbidden_pairs)) {
          blocks[i, ] <- new_bad
          blocks[j, ] <- new_donor
          break
        }
      }
      if (!.block_has_forbidden(blocks[i, ], forbidden_pairs)) break
    }
  }
  if (length(bad_idx()))
    warning(glue::glue("repair_forbidden: {length(bad_idx())} violations remain after {max_iter} iterations."))
  blocks
}

.optimize_position <- function(blocks, n_items, max_iter = 5000) {
  K       <- ncol(blocks)
  B       <- nrow(blocks)
  crit_fn <- function(b) {
    freq <- as.vector(vapply(seq_len(K), function(k) tabulate(b[, k], nbins = n_items), integer(n_items)))
    if (length(freq) != n_items * K) return(Inf)
    abs(sqrt(n_items * K) - sum(svd(freq)$u))
  }
  best_crit <- crit_fn(blocks)
  for (i in seq_len(max_iter)) {
    ii         <- sample.int(B, 1)
    try_blocks <- blocks
    try_blocks[ii, ] <- sample(try_blocks[ii, ])
    try_crit   <- crit_fn(try_blocks)
    if (try_crit < best_crit) { blocks <- try_blocks; best_crit <- try_crit }
  }
  blocks
}
