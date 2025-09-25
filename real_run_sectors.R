# real_run_sectors.R — End-to-end run on real sector data (merged/)
# Usage:
#   Rscript real_run_sectors.R data_dir=merged out_dir=output/real_merged k_ci=1 alpha=0.05 [max_rows=..] [max_cols=..]
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
  library(reshape2)
  library(gridExtra)
})

# ---- parse CLI args ----
args <- commandArgs(trailingOnly = TRUE)
par <- list(
  data_dir = "Asm",
  out_dir  = "output/real_Asm",
  k_ci     = 1,
  alpha    = 0.05,
  max_rows = NA_integer_,
  max_cols = NA_integer_      # 限制股票列数（不含日期列）；仅当 loader 支持时才传
)
for (a in args) {
  kv <- strsplit(a, "=")[[1]]
  if (length(kv) == 2 && kv[1] %in% names(par)) {
    key <- kv[1]; val <- kv[2]
    if (key %in% c("k_ci"))          par[[key]] <- as.integer(val)
    else if (key %in% c("alpha"))    par[[key]] <- as.numeric(val)
    else if (key %in% c("max_rows","max_cols")) par[[key]] <- as.integer(val)
    else par[[key]] <- val
  }
}
dir.create(par$out_dir, showWarnings = FALSE, recursive = TRUE)
if (!exists("h_min")) h_min <- 0.021  # fallback if not set in config.R

cat("[Start] Loading sector CSVs from: ", par$data_dir, "\n")
if (!is.na(par$max_cols)) {
  cat(sprintf("[Info] Limiting each sector to the first %d stock columns (if supported by loader).\n", par$max_cols))
}

# ---- build args safely depending on the loader signature ----
args_loader <- list(
  data_dir = par$data_dir,
  date_col = "Date",
  assume_price_not_log = TRUE
)
if ("max_rows" %in% names(formals(load_sectors_folder_as_sim_data))) {
  args_loader$max_rows <- if (is.na(par$max_rows)) NULL else par$max_rows
}
if ("max_cols" %in% names(formals(load_sectors_folder_as_sim_data))) {
  args_loader$max_cols <- if (is.na(par$max_cols)) NULL else par$max_cols
} else if (!is.na(par$max_cols)) {
  warning("The loader does not support 'max_cols'; ignoring this CLI argument.")
}

# 可选调试：确认生效的函数签名（需要时取消注释）
# cat("[DEBUG] loader formals: ", paste(names(formals(load_sectors_folder_as_sim_data)), collapse=", "), "\n")

sim_data <- do.call(load_sectors_folder_as_sim_data, args_loader)

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

# ---- save core outputs ----
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

# ---- one-sector visualization with CI ----
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

# ---- helpers ----
sanitize <- function(x) {
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

# ====== Eigenfunction & eigenvalue visualizations (μ part) ======
L_eval <- min(max(3L, mu_res$L), 8L, ncol(mu_res$PCs_hat))  # 至少3个，最多8个

df_phi_mu <- do.call(rbind, lapply(1:L_eval, function(j) {
  data.frame(t = sim_data$ts,
             value = mu_res$PCs_hat[, j],
             PC = paste0("PC", j))
}))
p_phi_mu <- ggplot(df_phi_mu, aes(t, value)) +
  geom_line(linewidth = 1, color = "red") +
  facet_wrap(~ PC, scales = "free_y") +
  labs(title = sprintf("Eigenfunctions (μ), first %d components", L_eval),
       x = "t", y = expression(phi[j](t))) +
  theme_minimal()
ggsave(file.path(par$out_dir, "eigenfunctions_mu_topL.png"),
       plot = p_phi_mu, width = 12, height = 6, dpi = 150)

lam_mu <- as.numeric(mu_res$lams_hat)
df_scree_mu <- data.frame(
  j = seq_along(lam_mu),
  lambda = lam_mu,
  FVE = cumsum(lam_mu) / sum(lam_mu)
)
p_mu1 <- ggplot(df_scree_mu, aes(j, lambda)) +
  geom_point() + geom_line() +
  labs(title = "Scree (μ): eigenvalues", x = "Component j", y = expression(lambda[j])) +
  theme_minimal()
p_mu2 <- ggplot(df_scree_mu, aes(j, FVE)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  geom_vline(xintercept = mu_res$L, color = "red", linetype = "dotted") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = sprintf("Cumulative FVE (μ), L selected = %d", mu_res$L),
       x = "Component j", y = "FVE") +
  theme_minimal()
ggsave(file.path(par$out_dir, "scree_mu.png"),
       plot = grid.arrange(p_mu1, p_mu2, ncol = 2), width = 12, height = 5, dpi = 150)

# 另存一份带 FVE 的 CSV
write.csv(df_scree_mu, file.path(par$out_dir, "eigenvalues_mu_with_FVE.csv"), row.names = FALSE)
write.csv(mu_res$lams_hat[1:L_eval], file.path(par$out_dir, "eigenvalues_mu_topL.csv"), row.names = FALSE)

# (C) scores（投影系数）
get_Z_for_mu <- function(sim_data) {
  if (max(sim_data$n_vec) == 1) sim_data$X_Delta / sqrt(sim_data$delta)
  else                          cluster_mean(sim_data$X_Delta, length(sim_data$n_vec), sim_data$n_vec) / sqrt(sim_data$delta)
}
Z_mu <- get_Z_for_mu(sim_data)
Xi_mu <- matrix(NA_real_, nrow = ncol(Z_mu), ncol = mu_res$L)
for (i in 1:ncol(Z_mu)) {
  for (j in 1:mu_res$L) {
    Xi_mu[i, j] <- mean((Z_mu[, i] - mu_res$m_mu_hat) * mu_res$PCs_hat[, j])
  }
}
colnames(Xi_mu) <- paste0("PC", seq_len(mu_res$L))
rownames(Xi_mu) <- if (!is.null(sim_data$sector_names) && length(sim_data$sector_names)==ncol(Z_mu))
  sim_data$sector_names else paste0("sector_", seq_len(ncol(Z_mu)))
write.csv(Xi_mu, file.path(par$out_dir, "scores_mu_topL.csv"))
cat("[Plot] Saved eigenfunctions & scree (μ). Exported scores_mu_topL.csv\n")

# ====== Centered-by-mean curves: μ_k(t) - m_μ(t) ======
cat("[Step] Plotting centered curves: mu_k(t) - m_mu(t) ...\n")
res_mu <- sweep(mu_res$mu_hat, 1, mu_res$m_mu_hat, FUN = "-")
write.csv(res_mu, file.path(par$out_dir, "mu_centered_by_mmu.csv"), row.names = FALSE)

for (k in seq_len(K)) {
  df_k <- data.frame(t = sim_data$ts, Residual = res_mu[, k])
  p_k <- ggplot(df_k, aes(x = t, y = Residual)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_line() +
    labs(
      title = sprintf("Centered μ_k(t): %s", sanitize(sec_names[k])),
      x = "t",
      y = expression(mu[k](t) - m[mu](t))
    ) +
    theme_minimal()
  ggsave(
    filename = file.path(par$out_dir, sprintf("mu_centered_%02d_%s.png", k, sanitize(sec_names[k]))),
    plot = p_k, width = 8, height = 5, dpi = 150
  )
}
pdf(file.path(par$out_dir, "mu_centered_all.pdf"), width = 8, height = 5)
for (k in seq_len(K)) {
  df_k <- data.frame(t = sim_data$ts, Residual = res_mu[, k])
  print(
    ggplot(df_k, aes(x = t, y = Residual)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_line() +
      labs(
        title = sprintf("Centered μ_k(t): %s", sanitize(sec_names[k])),
        x = "t",
        y = expression(mu[k](t) - m[mu](t))
      ) +
      theme_minimal()
  )
}
dev.off()
df_center <- do.call(rbind, lapply(seq_len(K), function(k) {
  data.frame(t = sim_data$ts, Residual = res_mu[, k], Sector = sec_names[k])
}))
p_center <- ggplot(df_center, aes(x = t, y = Residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_line() +
  facet_wrap(~ Sector, scales = "free_y") +
  labs(
    title = "Centered curves by sector: μ_k(t) − m_μ(t)",
    x = "t",
    y = expression(mu[k](t) - m[mu](t))
  ) +
  theme_minimal()
ggsave(file.path(par$out_dir, "mu_centered_facets.png"),
       plot = p_center, width = 14, height = max(6, ceiling(K/3) * 2.5), dpi = 150)
cat("[Plot] Saved facet overview: mu_centered_facets.png\n")

# ====== Visualize m_mu_hat & Raw vs Centered ======
cat("[Step] Visualizing m_mu_hat and raw vs centered comparisons ...\n")
df_mmu <- data.frame(t = sim_data$ts, m_mu = as.numeric(mu_res$m_mu_hat))
p_mmu <- ggplot(df_mmu, aes(x = t, y = m_mu)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_line(linewidth = 1) +
  labs(title = expression(paste("Estimated overall mean ", m[mu], "(t)")),
       x = "t", y = expression(m[mu](t))) +
  theme_minimal()
ggsave(file.path(par$out_dir, "m_mu_hat.png"),
       plot = p_mmu, width = 8, height = 4.5, dpi = 150)

k_show <- par$k_ci
if (k_show < 1 || k_show > ncol(mu_res$mu_hat)) k_show <- 1
df_one <- data.frame(
  t = sim_data$ts,
  Raw = mu_res$mu_hat[, k_show],
  Centered = res_mu[, k_show]
)
p_one <- ggplot(df_one, aes(x = t)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_line(aes(y = Raw), alpha = 0.9) +
  geom_line(aes(y = Centered), linetype = "dotted") +
  labs(title = sprintf("Sector %s: raw μ_k(t) vs centered by m_μ(t)",
                       sanitize(sec_names[k_show])),
       x = "t", y = "value") +
  theme_minimal()
ggsave(file.path(par$out_dir, sprintf("mu_raw_vs_centered_k%d_%s.png",
                                      k_show, sanitize(sec_names[k_show]))),
       plot = p_one, width = 8, height = 4.5, dpi = 150)

df_raw  <- melt(data.frame(t = sim_data$ts, mu_res$mu_hat), id.vars = "t",
                variable.name = "Sector", value.name = "Value")
df_cent <- melt(data.frame(t = sim_data$ts, res_mu), id.vars = "t",
                variable.name = "Sector", value.name = "Value")
df_raw$Type  <- "Raw"
df_cent$Type <- "Centered"
df_cmp <- rbind(df_raw, df_cent)
if (!is.null(sim_data$sector_names) && length(sim_data$sector_names) == K) {
  levels(df_cmp$Sector) <- sim_data$sector_names
}
p_cmp <- ggplot(df_cmp, aes(x = t, y = Value)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_line() +
  facet_grid(Type ~ Sector, scales = "free_y") +
  labs(title = "Raw vs Centered by m_μ(t) across sectors",
       x = "t", y = "") +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 7))
ggsave(file.path(par$out_dir, "mu_raw_vs_centered_facets.png"),
       plot = p_cmp, width = max(12, 3 + 1.2 * K), height = 6, dpi = 150)

mmu_mean <- mean(mu_res$m_mu_hat)
mmu_sd   <- sd(mu_res$m_mu_hat)
mu_sd_by_sector <- apply(mu_res$mu_hat, 2, sd)
scale_ratio <- mmu_sd / median(mu_sd_by_sector)
diag_tbl <- data.frame(
  metric = c("m_mu_mean", "m_mu_sd", "median_sd_of_mu_k", "sd_ratio_mmu_vs_mu"),
  value  = c(mmu_mean, mmu_sd, median(mu_sd_by_sector), scale_ratio)
)
write.csv(diag_tbl, file.path(par$out_dir, "m_mu_diagnostics.csv"), row.names = FALSE)
cat(sprintf("[Diag] m_mu_hat: mean=%.4g, sd=%.4g ; median sd of mu_k=%.4g ; sd ratio=%.4g\n",
            mmu_mean, mmu_sd, median(mu_sd_by_sector), scale_ratio))

