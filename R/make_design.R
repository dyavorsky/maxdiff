make_design <- function(n_items,
                        n_versions       = 100,
                        n_tasks          = 12,
                        n_items_per_task = 3,
                        forbidden_pairs  = NULL,
                        seed             = 42) {

  if (is.null(n_items)) stop("n_items is required.")
  if (n_items_per_task >= n_items)
    stop(glue::glue("n_items_per_task ({n_items_per_task}) must be less than n_items ({n_items})."))

  set.seed(seed)
  M            <- n_items
  K            <- n_items_per_task
  total_blocks <- n_versions * n_tasks

  message(glue::glue(
    "Generating design via choiceDes::tradeoff.des() for M={M}, B={total_blocks}, K={K} ..."
  ))
  if (total_blocks * M > 10000) message("  (Large study — this may take several minutes.)")
  flush.console()

  t0 <- proc.time()[["elapsed"]]
  td <- choiceDes::tradeoff.des(items = M, shown = K,
                                vers  = n_versions, tasks = n_tasks,
                                print = FALSE)
  elapsed <- round(proc.time()[["elapsed"]] - t0, 1)
  message(glue::glue("  ... tradeoff.des completed in {elapsed} seconds."))

  blocks <- as.matrix(td$design[, paste0("item", 1:K)])
  storage.mode(blocks) <- "integer"

  if (!is.null(forbidden_pairs) && length(forbidden_pairs)) {
    blocks <- .repair_forbidden(blocks, forbidden_pairs)
    blocks <- .optimize_position(blocks, M, max_iter = max(5000, 5 * total_blocks))
  }

  versions <- rep(seq_len(n_versions), each = n_tasks)
  tasks    <- rep(seq_len(n_tasks),    times = n_versions)
  out      <- tibble::tibble(Version = as.integer(versions), Task = as.integer(tasks))
  for (k in seq_len(K)) out[[paste0("Item", k)]] <- as.integer(blocks[, k])

  bm <- .balance_matrices(blocks, M)
  attr(out, "balance") <- list(
    freq_sd = stats::sd(bm$freq),
    pair_sd = stats::sd(bm$pair[lower.tri(bm$pair)]),
    pos_sd  = mean(apply(bm$pos, 2, stats::sd))
  )
  attr(out, "forbidden_pairs") <- forbidden_pairs
  out
}
