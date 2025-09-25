# R/realdata_sectors.R
# Utilities to load sector CSVs (merged/) and build sim_data for estimation.
# Assumptions:
#   - CSV has a date column (date_col) and stock PRICE columns (not log).
#   - Data are pre-cleaned; we do NOT forward/backward fill or impute.
#   - We only align by common trading dates, take logs, and build X_Delta.
# Output fields:
#   ts, delta, n_vec, X_Delta, sector_names (no aggregation here)

# ---- helpers (lean, no imputation) ----
.as_numeric_matrix <- function(df_like) {
  M <- as.matrix(df_like)
  storage.mode(M) <- "double"
  M
}

# ---- core loader ----
# data_dir: path to 'merged' folder
# date_col: name of the date column (e.g., "Date" or "date")
# assume_price_not_log: if TRUE take log() before variation()
# max_rows / max_cols: optional throttling for light test runs
load_sectors_folder_as_sim_data <- function(data_dir,
                                            date_col = "date",
                                            assume_price_not_log = TRUE,
                                            max_rows = NULL,
                                            max_cols = NULL) {
  stopifnot(dir.exists(data_dir))
  files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0L) stop("No CSV files under '", data_dir, "'.")
  cat("[Info] Found", length(files), "sector CSVs.\n")

  date_lists   <- list()
  value_mats   <- list()
  sector_names <- basename(files)
  kept_names   <- character(0)

  for (idx in seq_along(files)) {
    f <- files[idx]
    df <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
                   error = function(e) stop("Failed to read: ", f, " :: ", e$message))

    if (!(date_col %in% names(df))) {
      stop("File '", f, "' must contain a date column named: ", date_col)
    }

    # robust parse of dates
    dt_raw <- df[[date_col]]
    dt <- suppressWarnings(try(as.Date(dt_raw), silent = TRUE))
    if (inherits(dt, "try-error") || anyNA(dt)) {
      dt <- suppressWarnings(try(as.Date(dt_raw, format = "%Y/%m/%d"), silent = TRUE))
    }
    if (inherits(dt, "try-error") || anyNA(dt)) {
      dt <- suppressWarnings(try(as.Date(dt_raw, format = "%Y-%m-%d"), silent = TRUE))
    }
    if (inherits(dt, "try-error") || anyNA(dt)) {
      stop("Unrecognized date format in file: ", f)
    }

    # values: all columns except date
    val_names_all <- setdiff(names(df), date_col)
    if (length(val_names_all) == 0L) {
      warning("No stock columns in file: ", f, ". Skipped.")
      next
    }

    # limit to first max_cols (if provided)
    val_names <- if (!is.null(max_cols)) head(val_names_all, max_cols) else val_names_all
    vals_raw  <- df[, val_names, drop = FALSE]
    vals      <- .as_numeric_matrix(vals_raw)

    # 强约束：真实数据阶段不做填充。若还有 NA，直接报错提示清洗数据源。
    if (anyNA(vals)) {
      bad_cols <- val_names[colSums(is.na(vals)) > 0L]
      stop(sprintf("NA detected in %s. Please pre-clean your data. Columns with NA: %s",
                   basename(f), paste(head(bad_cols, 10L), collapse = ", ")))
    }

    date_lists[[length(date_lists) + 1L]] <- dt
    value_mats[[length(value_mats) + 1L]] <- vals
    kept_names <- c(kept_names, sector_names[idx])
    cat(sprintf("[OK] %s -> %d rows x %d stocks (as-is)\n",
                sector_names[idx], length(dt), ncol(vals)))
  }
  if (length(value_mats) == 0L) stop("No usable sector files.")

  # align by intersection of dates
  common_dates <- Reduce(intersect, date_lists)
  common_dates <- sort(common_dates)

  # optional throttling by rows
  if (!is.null(max_rows)) {
    if (length(common_dates) < max_rows) {
      warning("max_rows (", max_rows, ") > available rows (", length(common_dates), "). Use all.")
    } else {
      common_dates <- head(common_dates, max_rows)
      cat("[Info] Using first", length(common_dates), "rows (dates) for a light test run.\n")
    }
  }
  if (length(common_dates) < 2L) stop("Too few common dates across sectors.")
  cat("[Info] Common trading dates:", length(common_dates), "\n")

  # build matrix by sector order
  X_list <- list()
  n_vec  <- integer(length(value_mats))
  for (i in seq_along(value_mats)) {
    dt <- date_lists[[i]]
    M  <- value_mats[[i]]
    pos <- match(common_dates, dt)
    if (anyNA(pos)) stop("Date alignment failed for a sector file.")
    Mi <- M[pos, , drop = FALSE]
    X_list[[i]] <- Mi
    n_vec[i]    <- ncol(Mi)
  }
  Xmat_raw <- do.call(cbind, X_list)  # (m+1) x sum(n_k)
  K <- length(n_vec)
  cat(sprintf("[Info] K(sectors)=%d; n_k range: [%d, %d]; total series=%d\n",
              K, min(n_vec), max(n_vec), sum(n_vec)))

  # time grid & Δ
  m <- nrow(Xmat_raw) - 1L
  if (!exists("T")) T <<- 1  # fallback if config.R not sourced yet
  delta <- T / m
  ts    <- seq(delta, T, length.out = m)

  # log then variation -> per-stock ΔlogP (no aggregation, no √Δ scaling here)
  Xlog    <- if (assume_price_not_log) log(Xmat_raw) else Xmat_raw
  X_Delta <- variation(Xlog, m)  # (m x sum(n_k))

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
