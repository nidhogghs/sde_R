#!/usr/bin/env Rscript
# plot_sector_means.R — visualize averaged ΔlogP trajectories by sector.
#Rscript plot_sector_means.R data_dir=Asm out_dir=output/test_sector_means row_start=1 row_end=217  
# Usage example:
#   Rscript plot_sector_means.R data_dir=merged out_dir=output/sector_means \
#       row_start=10 row_end=120 col_start=5 col_end=80
#
# Arguments:
#   data_dir   : folder that contains sector CSVs (same format as merged/)
#   out_dir    : where to write CSV + PNG outputs
#   row_start  : (optional) starting row index in the raw price table (1-based)
#   row_end    : (optional) ending   row index in the raw price table (1-based)
#   col_start  : (optional) starting global stock column index (1-based)
#   col_end    : (optional) ending   global stock column index (1-based)
#
# Notes:
#   - row_start / row_end are converted to the ΔlogP index space (ts indices)
#     after taking first differences; i.e. index i corresponds to the Δ between
#     rows (i) and (i+1) of the raw price table.
#   - col_start / col_end operate on the concatenated stock columns across all
#     sectors. Column slicing respects sector block boundaries and updates the
#     bookkeeping (n_vec, sector_names) so downstream summaries are coherent.

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

source("config.R")
source("R/simulation.R")
source("R/utils.R")
source("R/realdata_sectors.R")

# ----- parse CLI arguments ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
par <- list(
  data_dir  = "merged",
  out_dir   = "output/sector_means",
  row_start = NA_integer_,
  row_end   = NA_integer_,
  col_start = NA_integer_,
  col_end   = NA_integer_
)

if (length(args) > 0) {
  for (a in args) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    key <- kv[1]
    val <- kv[2]
    if (!nzchar(key) || !(key %in% names(par))) next

    if (key %in% c("row_start", "row_end", "col_start", "col_end")) {
      par[[key]] <- if (nzchar(val)) as.integer(val) else NA_integer_
    } else {
      par[[key]] <- val
    }
  }
}

dir.create(par$out_dir, showWarnings = FALSE, recursive = TRUE)

# ----- load sector data -------------------------------------------------------
args_loader <- list(
  data_dir = par$data_dir,
  date_col = "Date",
  assume_price_not_log = TRUE
)

sim_data <- do.call(load_sectors_folder_as_sim_data, args_loader)

# ----- restrict rows (ΔlogP time index) --------------------------------------
m_total <- length(sim_data$ts)
row_idx_default <- c(1L, m_total)
row_idx_used <- row_idx_default

convert_row_idx <- function(idx_raw, default_val, upper_bound) {
  if (is.na(idx_raw)) {
    return(default_val)
  }
  idx_ts <- idx_raw - 1L
  if (idx_ts < 1L) idx_ts <- 1L
  if (idx_ts > upper_bound) idx_ts <- upper_bound
  idx_ts
}

row_idx_used[1] <- convert_row_idx(par$row_start, row_idx_default[1], m_total)
row_idx_used[2] <- convert_row_idx(par$row_end,   row_idx_default[2], m_total)
if (row_idx_used[1] > row_idx_used[2]) {
  stop(sprintf("Invalid row restriction: start (%d) > end (%d) after conversion.",
               row_idx_used[1], row_idx_used[2]))
}

if (!all(row_idx_used == row_idx_default)) {
  sim_data$ts      <- sim_data$ts[row_idx_used[1]:row_idx_used[2]]
  sim_data$X_Delta <- sim_data$X_Delta[row_idx_used[1]:row_idx_used[2], , drop = FALSE]
  cat(sprintf("[Info] Using ΔlogP rows %d-%d (ts indices) out of %d available.\n",
              row_idx_used[1], row_idx_used[2], m_total))
} else {
  cat(sprintf("[Info] Using all ΔlogP rows (1-%d).\n", m_total))
}

# ----- restrict columns (stock blocks) ---------------------------------------
orig_n_vec    <- sim_data$n_vec
orig_names    <- sim_data$sector_names
if (is.null(orig_names) || length(orig_names) != length(orig_n_vec)) {
  orig_names <- paste0("sector_", seq_along(orig_n_vec))
}

total_cols_before <- sum(orig_n_vec)
col_idx_default   <- c(1L, total_cols_before)
col_idx_used      <- col_idx_default

normalize_col_idx <- function(idx_raw, default_val, lower_bound, upper_bound) {
  if (is.na(idx_raw)) {
    return(default_val)
  }
  idx_val <- idx_raw
  if (idx_val < lower_bound) idx_val <- lower_bound
  if (idx_val > upper_bound) idx_val <- upper_bound
  idx_val
}

col_idx_used[1] <- normalize_col_idx(par$col_start, col_idx_default[1], 1L, total_cols_before)
col_idx_used[2] <- normalize_col_idx(par$col_end,   col_idx_default[2], 1L, total_cols_before)
if (col_idx_used[1] > col_idx_used[2]) {
  stop(sprintf("Invalid column restriction: start (%d) > end (%d) after clamping.",
               col_idx_used[1], col_idx_used[2]))
}

if (!all(col_idx_used == col_idx_default)) {
  keep_blocks   <- list()
  new_n_vec     <- integer(0)
  new_names     <- character(0)
  kept_indices  <- integer(0)
  col_cursor    <- 0L

  for (k in seq_along(orig_n_vec)) {
    block_indices <- seq.int(col_cursor + 1L, col_cursor + orig_n_vec[k])
    col_cursor    <- col_cursor + orig_n_vec[k]
    keep_global   <- block_indices[block_indices >= col_idx_used[1] &
                                   block_indices <= col_idx_used[2]]
    if (length(keep_global) == 0L) next
    keep_local <- keep_global - block_indices[1] + 1L
    keep_blocks[[length(keep_blocks) + 1L]] <-
      sim_data$X_Delta[, block_indices[keep_local], drop = FALSE]
    new_n_vec <- c(new_n_vec, length(keep_local))
    new_names <- c(new_names, orig_names[k])
    kept_indices <- c(kept_indices, keep_global)
  }

  if (length(keep_blocks) == 0L) {
    stop(sprintf("Column restriction [%d, %d] removed all stock series.",
                 col_idx_used[1], col_idx_used[2]))
  }

  sim_data$X_Delta     <- do.call(cbind, keep_blocks)
  sim_data$n_vec       <- new_n_vec
  sim_data$sector_names <- new_names
  col_idx_used <- c(min(kept_indices), max(kept_indices))

  cat(sprintf("[Info] Using stock columns %d-%d (global indices). Retained %d sectors / %d stocks.\n",
              col_idx_used[1], col_idx_used[2], length(new_n_vec), sum(new_n_vec)))
} else {
  cat(sprintf("[Info] Using all stock columns (1-%d) across %d sectors.\n",
              total_cols_before, length(orig_n_vec)))
  sim_data$sector_names <- orig_names
}

# ----- compute sector means ---------------------------------------------------
K_eff <- length(sim_data$n_vec)
if (K_eff == 0L) {
  stop("No sectors available after applying row/column restrictions.")
}

Z_Delta <- cluster_mean(sim_data$X_Delta, K_eff, sim_data$n_vec)
colnames(Z_Delta) <- sim_data$sector_names

out_df <- data.frame(t = sim_data$ts, Z_Delta, check.names = FALSE)

row_tag <- if (all(row_idx_used == row_idx_default)) {
  sprintf("1-%d", row_idx_default[2])
} else {
  sprintf("%d-%d", row_idx_used[1], row_idx_used[2])
}
col_tag <- if (all(col_idx_used == col_idx_default)) {
  sprintf("1-%d", col_idx_default[2])
} else {
  sprintf("%d-%d", col_idx_used[1], col_idx_used[2])
}

csv_path <- file.path(par$out_dir,
                      sprintf("sector_means_rows_%s_cols_%s.csv", row_tag, col_tag))
write.csv(out_df, csv_path, row.names = FALSE)
cat(sprintf("[Write] Sector means saved: %s\n", csv_path))

df_long <- melt(out_df, id.vars = "t", variable.name = "Sector",
                value.name = "MeanDelta")

plot_title <- sprintf("Sector mean ΔlogP (rows %s, cols %s)", row_tag, col_tag)

p <- ggplot(df_long, aes(x = t, y = MeanDelta, colour = Sector)) +
  geom_line(linewidth = 0.6) +
  labs(title = plot_title, x = "t", y = expression(bar(Delta*log(P)))) +
  theme_minimal()

png_path <- file.path(par$out_dir,
                      sprintf("sector_means_rows_%s_cols_%s.png", row_tag, col_tag))

ggsave(filename = png_path, plot = p, width = 10, height = 6, dpi = 150)
cat(sprintf("[Plot] Figure saved: %s\n", png_path))

cat("[Done] Sector mean computation complete.\n")Rscript plot_sector_means.R data_dir=merged out_dir=output/test_sector_means row_start=10 row_end=50 col_start=1 col_end=40
