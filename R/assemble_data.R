assemble_data <- function(responses, design,
                          version_col, best_cols, worst_cols,
                          anchor_cols       = NULL,
                          anchor_codes      = c(all = 1, some = 2, none = 99),
                          response_encoding = c("item", "position")) {

  response_encoding <- match.arg(response_encoding)
  anchored <- !is.null(anchor_cols) && length(anchor_cols) > 0L

  stop_if_not_dataframe(responses)
  stop_if_not_dataframe(design)
  require_columns(responses, c(version_col, unname(best_cols), unname(worst_cols)))
  if (anchored) require_columns(responses, unname(anchor_cols))

  task_col  <- if ("Task" %in% names(design)) "Task"
               else if ("Set" %in% names(design)) "Set"
               else stop("design must have a Task or Set column.")
  item_cols <- grep("^Item[0-9]+$", names(design), value = TRUE)
  if (length(item_cols) == 0L) stop("design must have columns named Item1, Item2, ...")

  T_n <- length(best_cols)
  K   <- length(item_cols)
  if (length(worst_cols) != T_n)
    stop("best_cols and worst_cols must have the same length.")
  if (!all(names(best_cols) == names(worst_cols)))
    stop("best_cols and worst_cols must be named by the same task numbers.")
  if (anchored && length(anchor_cols) != T_n)
    stop("anchor_cols must have the same length as best_cols.")
  if (anchored && !all(c("all", "some", "none") %in% names(anchor_codes)))
    stop("anchor_codes must be a named vector with elements 'all', 'some', 'none'.")

  M <- max(unlist(design[item_cols], use.names = FALSE), na.rm = TRUE)
  if (M < 2L) stop("Design must have at least 2 distinct items.")

  n_var <- if (anchored) M else M - 1L

  item_row <- function(m, s = 1) {
    row <- numeric(n_var)
    if (anchored) {
      row[m] <- s
    } else if (m < M) {
      row[m] <- s
    }
    row
  }

  design_key   <- paste(design$Version, design[[task_col]], sep = "_")
  items_by_row <- as.matrix(design[item_cols])
  storage.mode(items_by_row) <- "integer"

  versions   <- as.integer(responses[[version_col]])
  best_mat   <- vapply(best_cols,  function(cn) as.integer(responses[[cn]]),
                       integer(nrow(responses)))
  worst_mat  <- vapply(worst_cols, function(cn) as.integer(responses[[cn]]),
                       integer(nrow(responses)))
  anchor_mat <- if (anchored) {
    vapply(anchor_cols, function(cn) as.integer(responses[[cn]]),
           integer(nrow(responses)))
  } else NULL

  in_design   <- !is.na(versions) & as.character(versions) %in% design$Version
  complete_bw <- !apply(is.na(best_mat),  1, any) &
                 !apply(is.na(worst_mat), 1, any)
  keep_idx    <- which(in_design & complete_bw)
  N <- length(keep_idx)
  if (N == 0L) stop("No respondents passed validation (version match + complete best/worst).")

  task_numbers <- as.integer(names(best_cols))

  lgtdata <- vector("list", N)
  for (j in seq_len(N)) {
    r <- keep_idx[[j]]
    v <- versions[[r]]

    max_occ_per_task <- if (anchored) 2L + K else 2L
    y_buf <- integer(T_n * max_occ_per_task)
    X_buf <- vector("list", T_n * max_occ_per_task)
    n_occ <- 0L

    for (i in seq_len(T_n)) {
      t    <- task_numbers[[i]]
      drow <- match(paste(v, t, sep = "_"), design_key)
      if (is.na(drow))
        stop(glue::glue("Design lookup missing for version={v}, task={t} (respondent row {r})."))
      items <- as.integer(items_by_row[drow, ])

      b_raw <- best_mat[r, i]
      w_raw <- worst_mat[r, i]

      if (response_encoding == "position") {
        bp <- b_raw; wp <- w_raw
        if (is.na(bp) || bp < 1L || bp > K || is.na(wp) || wp < 1L || wp > K)
          stop(glue::glue("Respondent row {r}, task {t}: best/worst position out of range 1..{K}."))
        b_item <- items[[bp]]; w_item <- items[[wp]]
      } else {
        bp <- match(b_raw, items); wp <- match(w_raw, items)
        if (is.na(bp) || is.na(wp)) {
          bad <- if (is.na(bp)) b_raw else w_raw
          stop(glue::glue("Respondent row {r}, task {t}: chosen item {bad} not in shown set."))
        }
        b_item <- b_raw; w_item <- w_raw
      }
      if (bp == wp)
        stop(glue::glue("Respondent row {r}, task {t}: best and worst pick the same item."))

      X_best <- do.call(rbind, lapply(items, item_row, s = 1))
      n_occ  <- n_occ + 1L
      y_buf[[n_occ]] <- bp
      X_buf[[n_occ]] <- X_best

      remaining   <- items[-bp]
      wp_in_rem   <- wp - (wp > bp)
      X_worst     <- do.call(rbind, lapply(remaining, item_row, s = -1))
      n_occ       <- n_occ + 1L
      y_buf[[n_occ]] <- wp_in_rem
      X_buf[[n_occ]] <- X_worst

      if (anchored) {
        a <- anchor_mat[r, i]
        if (!is.na(a)) {
          if (a == anchor_codes[["all"]]) {
            for (m in items) {
              X_anc <- rbind(item_row(m, s = 1), numeric(n_var))
              n_occ <- n_occ + 1L
              y_buf[[n_occ]] <- 1L
              X_buf[[n_occ]] <- X_anc
            }
          } else if (a == anchor_codes[["none"]]) {
            for (m in items) {
              X_anc <- rbind(item_row(m, s = 1), numeric(n_var))
              n_occ <- n_occ + 1L
              y_buf[[n_occ]] <- 2L
              X_buf[[n_occ]] <- X_anc
            }
          } else if (a == anchor_codes[["some"]]) {
            X_anc_b <- rbind(item_row(b_item, s = 1), numeric(n_var))
            n_occ   <- n_occ + 1L
            y_buf[[n_occ]] <- 1L
            X_buf[[n_occ]] <- X_anc_b

            X_anc_w <- rbind(item_row(w_item, s = 1), numeric(n_var))
            n_occ   <- n_occ + 1L
            y_buf[[n_occ]] <- 2L
            X_buf[[n_occ]] <- X_anc_w
          }
        }
      }
    }

    lgtdata[[j]] <- list(y = y_buf[seq_len(n_occ)], X = X_buf[seq_len(n_occ)])
  }

  list(
    lgtdata          = lgtdata,
    N                = N,
    T                = T_n,
    K                = K,
    M                = M,
    n_var            = n_var,
    anchored         = anchored,
    respondent_index = keep_idx
  )
}
