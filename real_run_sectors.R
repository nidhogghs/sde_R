# real_run_sectors.R — End-to-end run on real sector data (merged/).
# Usage:
#   Rscript real_run_sectors.R data_dir=merged out_dir=output/real_merged k_ci=1 alpha=0.05
# Notes:
#   - Assumes each CSV in data_dir has a 'date' column and stock close-price columns.
#   - Each file is a sector (cluster). n_k = number of stock columns in that file.
#   - Comments and console messages are in English.

# ---- load your modules ----
source("config.R")              # must define T, and optionally h_min
source("R/simulation.R")        # uses variation()
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/mu_sigma_est.R")
source("R/realdata_sectors.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- parse CLI args ----
args <- commandArgs(trailingOnly = TRUE)
par <- list(
  data_dir = "merged",
  out_dir  = "output/real_merged",
  k_ci     = 1,
  alpha    = 0.05,
  max_rows = NA_integer_         # <== 新增
)
for (a in args) {
  kv <- strsplit(a, "=")[[1]]
  if (length(kv) == 2 && kv[1] %in% names(par)) {
    key <- kv[1]; val <- kv[2]
    if (key %in% c("k_ci"))        par[[key]] <- as.integer(val)
    else if (key %in% c("alpha"))  par[[key]] <- as.numeric(val)
    else if (key %in% c("max_rows")) par[[key]] <- as.integer(val)  # <== 新增
    else par[[key]] <- val
  }
}
dir.create(par$out_dir, showWarnings = FALSE, recursive = TRUE)
if (!exists("h_min")) h_min <- 0.021  # fallback if not set in config.R

cat("[Start] Loading sector CSVs from: ", par$data_dir, "\n")
sim_data <- load_sectors_folder_as_sim_data(
  data_dir = par$data_dir,
  date_col = "date",
  assume_price_not_log = TRUE,
  max_rows = if (is.na(par$max_rows)) NULL else par$max_rows   # <== 新增
)
# quick sanity check
K <- length(sim_data$n_vec)
if (K != 18) {
  warning(sprintf("K(sectors)=%d (expected 18). Proceeding anyway.", K))
}

# ---- estimation & inference ----
cat("[Step] Estimating mu ...\n")
mu_res <- estimate_mu_from_data(sim_data)

cat("[Step] Estimating sigma^2 ...\n")
sigma_res <- inference_sigma_from_data(sim_data)

cat("[Step] Building pointwise CI for mu_k(t) ...\n")
ci <- compute_mu_ci(
  mu_hat     = mu_res$mu_hat,
  PCs        = mu_res$PCs_hat,
  sigma2_hat = sigma_res$sigma2_hat,
  ts         = sim_data$ts,
  n_vec      = sim_data$n_vec,
  alpha      = par$alpha,
  L          = mu_res$L
)

# ---- save outputs ----
write.csv(mu_res$mu_hat,    file.path(par$out_dir, "mu_hat.csv"),    row.names = FALSE)
write.csv(mu_res$m_mu_hat,  file.path(par$out_dir, "m_mu_hat.csv"),  row.names = FALSE)
write.csv(mu_res$G_mu_hat,  file.path(par$out_dir, "G_mu_hat.csv"),  row.names = FALSE)
write.csv(mu_res$PCs_hat,   file.path(par$out_dir, "PCs_mu.csv"),    row.names = FALSE)
write.csv(mu_res$lams_hat,  file.path(par$out_dir, "lams_mu.csv"),   row.names = FALSE)

write.csv(sigma_res$V_hat,       file.path(par$out_dir, "V_hat.csv"),       row.names = FALSE)
write.csv(sigma_res$sigma2_hat,  file.path(par$out_dir, "sigma2_hat.csv"),  row.names = FALSE)
write.csv(sigma_res$m_V_hat,     file.path(par$out_dir, "m_V_hat.csv"),     row.names = FALSE)
write.csv(sigma_res$G_V_hat,     file.path(par$out_dir, "G_V_hat.csv"),     row.names = FALSE)
write.csv(sigma_res$PCs_hat,     file.path(par$out_dir, "PCs_V.csv"),       row.names = FALSE)
write.csv(sigma_res$lams_hat,    file.path(par$out_dir, "lams_V.csv"),      row.names = FALSE)

write.csv(ci$lower, file.path(par$out_dir, "mu_ci_lower.csv"), row.names = FALSE)
write.csv(ci$upper, file.path(par$out_dir, "mu_ci_upper.csv"), row.names = FALSE)

# ---- visualization for one sector k ----
p <- plot_mu_with_ci(
  ts        = sim_data$ts,
  mu_true   = mu_res$mu_hat,  # no ground truth; plot estimate + CI
  mu_hat    = mu_res$mu_hat,
  ci_lower  = ci$lower,
  ci_upper  = ci$upper,
  k         = par$k_ci
)
ggsave(filename = file.path(par$out_dir, sprintf("mu_with_CI_k%d.png", par$k_ci)),
       plot = p, width = 8, height = 5, dpi = 150)

cat("[Done] Outputs saved to: ", par$out_dir, "\n")

# ---- per-sector plots: one PNG each + a multi-page PDF + a facet overview ----
sanitize <- function(x) {
  # make safe filenames on Windows/macOS/Linux
  x <- gsub('[/:*?"<>|]', "_", x)
  x <- gsub("\\s+", "_", x)
  x
}

K <- ncol(mu_res$mu_hat)
sec_names <- sim_data$sector_names
if (length(sec_names) != K) sec_names <- paste0("sector_", seq_len(K))

# (A) 每个行业一张：mu + 95%CI
for (k in seq_len(K)) {
  p_k <- plot_mu_with_ci(
    ts        = sim_data$ts,
    mu_true   = mu_res$mu_hat,  # 真实数据无真值，这里只画估计+CI
    mu_hat    = mu_res$mu_hat,
    ci_lower  = ci$lower,
    ci_upper  = ci$upper,
    k         = k
  )
  ggsave(
    filename = file.path(par$out_dir,
      sprintf("mu_with_CI_%02d_%s.png", k, sanitize(sec_names[k]))),
    plot = p_k, width = 8, height = 5, dpi = 150
  )
}
cat("[Plot] Saved per-sector μ plots with CI.\n")

# (B) 多页 PDF（每页一个行业）
pdf(file.path(par$out_dir, "mu_with_CI_all.pdf"), width = 8, height = 5)
for (k in seq_len(K)) {
  print(plot_mu_with_ci(
    ts        = sim_data$ts,
    mu_true   = mu_res$mu_hat,
    mu_hat    = mu_res$mu_hat,
    ci_lower  = ci$lower,
    ci_upper  = ci$upper,
    k         = k
  ))
}
dev.off()
cat("[Plot] Saved multi-page PDF: mu_with_CI_all.pdf\n")

# (C) 九宫格总览（facet）
suppressPackageStartupMessages(library(reshape2))
df_all <- do.call(rbind, lapply(seq_len(K), function(k) {
  data.frame(
    t       = sim_data$ts,
    Estimate= mu_res$mu_hat[, k],
    lower   = ci$lower[, k],
    upper   = ci$upper[, k],
    Sector  = sec_names[k]
  )
}))
p_facet <- ggplot(df_all, aes(x = t)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(aes(y = Estimate)) +
  facet_wrap(~ Sector, scales = "free_y") +
  labs(title = "Estimated μ_k(t) with 95% CI by Sector",
       x = "t", y = expression(mu[k](t))) +
  theme_minimal()
ggsave(file.path(par$out_dir, "mu_facets.png"),
       plot = p_facet, width = 14, height = max(6, ceiling(K/3) * 2.5), dpi = 150)
cat("[Plot] Saved facet overview: mu_facets.png\n")

