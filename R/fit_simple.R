fit_simple <- function(responses, design,
                       version_col, best_cols, worst_cols,
                       item_labels       = NULL,
                       response_encoding = c("item", "position")) {

  response_encoding <- match.arg(response_encoding)
  stop_if_not_dataframe(responses)
  stop_if_not_dataframe(design)
  require_columns(responses, c(version_col, unname(best_cols), unname(worst_cols)))

  task_col  <- if ("Task" %in% names(design)) "Task"
               else if ("Set" %in% names(design)) "Set"
               else stop("design must have a Task or Set column.")
  item_cols <- grep("^Item[0-9]+$", names(design), value = TRUE)
  if (length(item_cols) == 0L) stop("design must have columns named Item1, Item2, ...")
  if (length(best_cols) != length(worst_cols))
    stop("best_cols and worst_cols must have the same length.")
  if (!all(names(best_cols) == names(worst_cols)))
    stop("best_cols and worst_cols must be named by the same task numbers.")

  n_items      <- max(unlist(design[item_cols], use.names = FALSE), na.rm = TRUE)
  K            <- length(item_cols)
  labels       <- coerce_item_labels(item_labels, n_items)
  N            <- nrow(responses)
  task_numbers <- as.integer(names(best_cols))
  T_n          <- length(task_numbers)
  versions     <- as.integer(responses[[version_col]])

  long_version <- rep(versions,     times = T_n)
  long_task    <- rep(task_numbers, each  = N)
  long_best    <- unlist(lapply(best_cols,  function(c) as.integer(responses[[c]])),
                         use.names = FALSE)
  long_worst   <- unlist(lapply(worst_cols, function(c) as.integer(responses[[c]])),
                         use.names = FALSE)

  design_key   <- paste(design$Version, design[[task_col]], sep = "_")
  items_by_row <- as.matrix(design[item_cols])
  storage.mode(items_by_row) <- "integer"
  drow <- match(paste(long_version, long_task, sep = "_"), design_key)
  keep <- !is.na(drow)
  drow       <- drow[keep]
  long_best  <- long_best[keep]
  long_worst <- long_worst[keep]
  items_shown <- items_by_row[drow, , drop = FALSE]

  if (response_encoding == "position") {
    to_item <- function(pos) {
      out <- rep(NA_integer_, length(pos))
      ok  <- !is.na(pos) & pos >= 1L & pos <= K
      out[ok] <- items_shown[cbind(which(ok), pos[ok])]
      out
    }
    long_best  <- to_item(long_best)
    long_worst <- to_item(long_worst)
  }

  shown_count <- tabulate(as.vector(items_shown), nbins = n_items)
  best_count  <- tabulate(long_best,              nbins = n_items)
  worst_count <- tabulate(long_worst,             nbins = n_items)

  score <- ifelse(shown_count > 0L,
                  (best_count - worst_count) / shown_count,
                  NA_real_)

  out <- tibble::tibble(
    item       = seq_len(n_items),
    item_label = unname(labels[as.character(seq_len(n_items))]),
    n_shown    = shown_count,
    n_best     = best_count,
    n_worst    = worst_count,
    score      = score
  )
  class(out) <- c("maxdiff_simple", class(out))
  out
}

print.maxdiff_simple <- function(x, ...) {
  cat("<maxdiff_simple>  (simple best-minus-worst scores)\n")
  cat("  Items:", nrow(x), "\n\n")
  print(tibble::as_tibble(x), ...)
  invisible(x)
}
