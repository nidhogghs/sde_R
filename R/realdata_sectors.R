# R/realdata_sectors.R
# Utilities to load sector CSVs (merged/) and build sim_data for estimation.
# Requirements: use base R only. Comments + console messages in English.

# ---- small helpers ----
.na_locf <- function(v) {
  last <- NA_real_
  for (i in seq_along(v)) {
    if (!is.na(v[i])) last <- v[i] else v[i] <- last
  }
  v
}
.na_bfill <- function(v) {
  nextv <- NA_real_
  for (i in length(v):1) {
    if (!is.na(v[i])) nextv <- v[i] else v[i] <- nextv
  }
  v
}
.clean_numeric_matrix <- function(M, na_drop_thresh = 0.50, do_fill = TRUE) {
  # Convert to numeric, drop columns with too many NAs, fill remaining NAs.
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  keep <- rep(TRUE, ncol(M))
  for (j in seq_len(ncol(M))) {
    frac_na <- mean(is.na(M[, j]))
    if (is.na(frac_na) || frac_na > na_drop_thresh) {
      keep[j] <- FALSE
    } else if (do_fill && anyNA(M[, j])) {
      M[, j] <- .na_bfill(.na_locf(M[, j]))
    }
  }
  M[, keep, drop = FALSE]
}

# ---- core loader ----
# data_dir: path to 'merged' folder
# date_col: name of the first column (e.g., "date")
# returns: list(ts, delta, n_vec, X_Delta)
load_sectors_folder_as_sim_data <- function(data_dir,
                                            date_col = "date",
                                            assume_price_not_log = TRUE,
                                            max_rows = NULL) {
  stopifnot(dir.exists(data_dir))
  files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0L) stop("No CSV files under '", data_dir, "'.")
  cat("[Info] Found", length(files), "sector CSVs.\n")

  # 1) read all files, collect dates and matrices
  date_lists <- list()
  value_mats <- list()
  sector_names <- basename(files)
  kept_names <- character(0) 

  for (idx in seq_along(files)) {
    f <- files[idx]
    df <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
                   error = function(e) stop("Failed to read: ", f, " :: ", e$message))
    if (!(date_col %in% names(df))) {
      stop("File '", f, "' must contain a date column named: ", date_col)
    }

    # parse dates like "2023/4/21"
    dt_raw <- df[[date_col]]
    # robust parse
    dt <- tryCatch(as.Date(dt_raw),
                   error = function(e) NA)
    if (anyNA(dt)) {
      dt <- tryCatch(as.Date(dt_raw, format = "%Y/%m/%d"),
                     error = function(e) NA)
    }
    if (anyNA(dt)) {
      dt <- tryCatch(as.Date(dt_raw, format = "%Y-%m-%d"),
                     error = function(e) NA)
    }
    if (anyNA(dt)) stop("Unrecognized date format in file: ", f)

    # values: all columns except date
    vals <- df[, setdiff(names(df), date_col), drop = FALSE]
    vals <- .clean_numeric_matrix(vals, na_drop_thresh = 0.50, do_fill = TRUE)
    if (ncol(vals) == 0L) {
      warning("All columns dropped for file: ", f, " (too many NAs). Skipped.")
      next
    }

    date_lists[[length(date_lists) + 1]] <- dt
    value_mats[[length(value_mats) + 1]] <- vals
    kept_names <- c(kept_names, sector_names[idx])
    cat(sprintf("[OK] %s -> %d rows x %d stocks\n", sector_names[idx], length(dt), ncol(vals)))
  }
  if (length(value_mats) == 0L) stop("No usable sector files after cleaning.")

  # 2) align dates across all sectors (intersection)
  common_dates <- Reduce(intersect, date_lists)
  common_dates <- sort(common_dates)

    # ===== 新增：仅取前 max_rows 行（按日期从早到晚）=====
  if (!is.null(max_rows)) {
    if (length(common_dates) < max_rows) {
      warning("max_rows (", max_rows, ") > available rows (", length(common_dates), "). Use all.")
    } else {
      common_dates <- head(common_dates, max_rows)
      cat("[Info] Using first", length(common_dates), "rows (dates) for a light test run.\n")
    }
  }
  # ================================================
  
  if (length(common_dates) < 2L) stop("Too few common dates across sectors.")
  cat("[Info] Common trading dates:", length(common_dates), "\n")

  # 3) build (m+1) x sum(n_k) matrix in sector order
  X_list <- list()
  n_vec <- integer(length(value_mats))
  for (i in seq_along(value_mats)) {
    dt <- date_lists[[i]]
    M  <- value_mats[[i]]
    # reorder rows by common_dates
    pos <- match(common_dates, dt)
    if (anyNA(pos)) stop("Date alignment failed for a sector file.")
    Mi <- M[pos, , drop = FALSE]
    X_list[[i]] <- Mi
    n_vec[i] <- ncol(Mi)
  }
  Xmat_raw <- do.call(cbind, X_list)  # (m+1) x sum(n_k)
  K <- length(n_vec)
  cat(sprintf("[Info] K(sectors)=%d; n_k range: [%d, %d]; total series=%d\n",
              K, min(n_vec), max(n_vec), sum(n_vec)))

  # 4) to sim_data (consistent with your pipeline): log -> variation -> X_Delta
  m <- nrow(Xmat_raw) - 1L
  if (!exists("T")) T <<- 1  # fallback if config.R not sourced yet
  delta <- T / m
  ts <- seq(delta, T, length.out = m)

  Xlog <- if (assume_price_not_log) log(Xmat_raw) else Xmat_raw
  # variation(log(X)) : (m x sum(n_k))
  X_Delta <- variation(Xlog, m)

  list(
    ts = ts,
    delta = delta,
    n_vec = n_vec,
    V_true = NULL,
    sigma2_true = NULL,
    mu_true = NULL,
    X_Delta = X_Delta,
    sector_names = kept_names  
  )
}
