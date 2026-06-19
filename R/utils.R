# Internal helpers for maxdiff. Not exported.

`%||%` <- function(x, y) if (!is.null(x)) x else y

# Auto-detect MaxDiff response columns from a wide data frame.
# Returns list(best_cols, worst_cols, anchor_cols), each a named character
# vector ordered by task number. Capture group 1 isolates the task number.
detect_md_columns <- function(df,
                              best_pattern   = "^FinalMD1_([0-9]+)most$",
                              worst_pattern  = "^FinalMD1_([0-9]+)least$",
                              anchor_pattern = "^MD2_([0-9]+)$") {
  nms <- names(df)
  pull <- function(pat) {
    if (is.null(pat)) return(setNames(character(0), character(0)))
    hits <- regmatches(nms, regexec(pat, nms))
    keep <- vapply(hits, length, integer(1)) > 0
    if (!any(keep)) return(setNames(character(0), character(0)))
    cols  <- nms[keep]
    tasks <- vapply(hits[keep], `[`, character(1), 2)
    ord   <- order(as.integer(tasks))
    setNames(cols[ord], tasks[ord])
  }
  list(
    best_cols   = pull(best_pattern),
    worst_cols  = pull(worst_pattern),
    anchor_cols = pull(anchor_pattern)
  )
}

stop_if_not_dataframe <- function(x, arg = deparse(substitute(x))) {
  if (!is.data.frame(x))
    stop(glue::glue("{arg} must be a data frame, not {class(x)[1]}."))
}

require_columns <- function(df, cols, arg = deparse(substitute(df))) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols))
    stop(glue::glue("{arg} is missing required columns: {paste(missing_cols, collapse = ', ')}."))
}

coerce_item_labels <- function(item_labels, n_items) {
  if (is.null(item_labels))
    return(setNames(paste0("Item ", seq_len(n_items)), seq_len(n_items)))
  if (is.character(item_labels) && length(item_labels) == n_items)
    return(setNames(item_labels, seq_len(n_items)))
  if (!is.null(names(item_labels))) {
    ids <- suppressWarnings(as.integer(names(item_labels)))
    if (any(is.na(ids))) stop("Names of item_labels must be integer item IDs.")
    out <- setNames(rep(NA_character_, n_items), seq_len(n_items))
    out[as.character(ids)] <- as.character(item_labels)
    missing_ids <- which(is.na(out))
    if (length(missing_ids)) out[missing_ids] <- paste0("Item ", missing_ids)
    return(out)
  }
  stop(glue::glue(
    "item_labels must be NULL, a length-{n_items} character vector, or a named character vector keyed by item ID."
  ))
}

# pandterm: stop with a message, no call stack (ported from bayesm).
pandterm <- function(message) stop(message, call. = FALSE)

# fsh: flush console (ported from bayesm).
fsh <- function() invisible(utils::flush.console())
