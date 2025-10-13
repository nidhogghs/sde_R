# =========================================
# convergence_analysis_bw.R
# -----------------------------------------
# Run simulations for K = 500 with varying n_ave to study convergence.
# Generates black & white plots (journal-friendly): no titles, linetype/shape only.
# Usage:
#   Rscript convergence_analysis_bw.R [num_seeds=100]
# Example (quick):
#   Rscript convergence_analysis_bw.R num_seeds=2
# =========================================

# ---- Load modules ----
source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/mu_sigma_est.R")
source("R/traditional.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

# ---- Parse args ----
args <- commandArgs(trailingOnly = TRUE)
num_seeds <- 1
for (arg in args) {
  if (grepl("^num_seeds=", arg)) {
    num_seeds <- as.integer(sub("num_seeds=", "", arg))
  }
}

# ---- Experiment setup ----
K <- 500
n_values <- c(1, 30, 70, 100)

results <- data.frame()
mu_store <- list()        # i = 1 curve across n_ave (seed == 1)
G_mu_store <- list()      # G_mu estimates (seed == 1)
eigenfun_store <- list()  # leading eigenfunctions (seed == 1)
eigenval_store <- list()  # leading eigenvalues (seed == 1)

# ---- Main loop over n_ave ----
for (n in n_values) {
  rmse_fda_list  <- numeric(num_seeds)
  rmse_trad_list <- numeric(num_seeds)
  rmse_G_list    <- numeric(num_seeds)

  for (seed in seq_len(num_seeds)) {
    set.seed(seed)

    sim_data <- prepare_simulation_data(
      K = K,
      n_ave = n,
      m = 50,
      k = 4,
      l = 4,
      scale_sigma = 5,
      seed_sigma = 311 + seed,
      seed_mu = 311 + seed,
      seed_X = 69 + seed
    )

    mu_res    <- estimate_mu_from_data(sim_data)
    mu_trad   <- estimate_mu_traditional_from_data(sim_data)
    sigma_res <- inference_sigma_from_data(sim_data)

    # --- RMSE for mu (FDA & Traditional) ---
    rmse_fda  <- sqrt(mean((mu_res$mu_hat - sim_data$mu_true)^2))
    rmse_trad <- sqrt(mean((mu_trad       - sim_data$mu_true)^2))
    rmse_fda_list[seed]  <- rmse_fda
    rmse_trad_list[seed] <- rmse_trad

    # --- RMSE for G_mu ---
    G_mu_true <- (sim_data$mu_true - rowMeans(sim_data$mu_true)) %*%
                 t(sim_data$mu_true - rowMeans(sim_data$mu_true)) / K
    rmse_G    <- sqrt(mean((mu_res$G_mu_hat - G_mu_true)^2))
    rmse_G_list[seed] <- rmse_G

    # --- Store first-seed detailed outputs for plotting ---
    if (seed == 1) {
      ts_vals <- sim_data$ts

      # i = 1 curve across n_ave
      mu_store[[as.character(n)]] <- data.frame(
        index = seq_len(nrow(mu_res$mu_hat)),
        time = ts_vals,
        True = sim_data$mu_true[, 1],
        FDA = mu_res$mu_hat[, 1],
        Traditional = mu_trad[, 1],
        n_ave = n
      )

      # G_mu entries
      m_dim <- nrow(mu_res$G_mu_hat)
      G_mu_store[[as.character(n)]] <- data.frame(
        i = rep(seq_len(m_dim), times = m_dim),
        j = rep(seq_len(m_dim), each = m_dim),
        G_mu_hat = as.vector(mu_res$G_mu_hat),
        G_mu_true = as.vector(G_mu_true),
        n_ave = n
      )

      # Leading eigenfunctions/eigenvalues (align sign to True)
      eig_true   <- eigen(G_mu_true)
      PCs_true   <- Re(eig_true$vectors) * sqrt(m_dim)
      lams_true  <- Re(eig_true$values) / m_dim
      PCs_hat    <- mu_res$PCs_hat
      lams_hat   <- mu_res$lams_hat

      num_components <- min(3, ncol(PCs_hat))
      if (num_components > 0) {
        PCs_true_sel <- PCs_true[, seq_len(num_components), drop = FALSE]
        PCs_hat_sel  <- PCs_hat[,  seq_len(num_components), drop = FALSE]

        for (comp in seq_len(num_components)) {
          align_sign <- sign(sum(PCs_hat_sel[, comp] * PCs_true_sel[, comp]))
          if (align_sign == 0) align_sign <- 1
          PCs_hat_sel[, comp] <- PCs_hat_sel[, comp] * align_sign
        }

        eigenfun_store[[as.character(n)]] <- rbind(
          data.frame(
            time = rep(ts_vals, times = num_components),
            component = factor(rep(seq_len(num_components), each = length(ts_vals)),
                               levels = seq_len(num_components),
                               labels = paste0("phi", seq_len(num_components))),
            value = as.vector(PCs_true_sel),
            Method = "True",
            n_ave = n
          ),
          data.frame(
            time = rep(ts_vals, times = num_components),
            component = factor(rep(seq_len(num_components), each = length(ts_vals)),
                               levels = seq_len(num_components),
                               labels = paste0("phi", seq_len(num_components))),
            value = as.vector(PCs_hat_sel),
            Method = "Estimated",
            n_ave = n
          )
        )

        eigenval_store[[as.character(n)]] <- data.frame(
          component = seq_len(num_components),
          value = c(lams_true[seq_len(num_components)], lams_hat[seq_len(num_components)]),
          Method = rep(c("True", "Estimated"), each = num_components),
          n_ave = n
        )
      }
    }
  }

  results <- rbind(results, data.frame(
    n_ave = n,
    avg_rmse_mu_fda         = mean(rmse_fda_list),
    avg_rmse_mu_traditional = mean(rmse_trad_list),
    avg_rmse_G_mu           = mean(rmse_G_list)
  ))
}

# ---- Save numeric tables ----
write.csv(results, "convergence_summary.csv", row.names = FALSE)
G_mu_df <- do.call(rbind, G_mu_store)
write.csv(G_mu_df, "G_mu_estimates_seed1.csv", row.names = FALSE)

# ---- Visualization (Black & White, no titles) ----

# 1) RMSE vs n_ave (linetype + shape; no color)
results_long <- melt(results, id.vars = "n_ave",
                     measure.vars = c("avg_rmse_mu_fda", "avg_rmse_mu_traditional"),
                     variable.name = "Method", value.name = "avg_rmse_mu")
results_long$Method <- factor(results_long$Method,
  levels = c("avg_rmse_mu_fda", "avg_rmse_mu_traditional"),
  labels = c("FDA", "Traditional")
)

p1_bw <- ggplot(results_long, aes(x = n_ave, y = avg_rmse_mu,
                                  linetype = Method, shape = Method)) +
  geom_line(color = "black") +
  geom_point(color = "black", size = 2) +
  scale_linetype_manual(name = "Legend",
                        values = c("FDA" = "longdash", "Traditional" = "solid")) +
  scale_shape_manual(name = "Legend",
                     values = c("FDA" = 17, "Traditional" = 16)) +  # triangle / solid circle
  labs(x = "n_ave", y = expression("Average RMSE(" * mu * ")")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

ggsave("rmse_convergence_bw.png", p1_bw, width = 5, height = 4, dpi = 300)


# 2) μ_i(t) by n_ave: facet black & white (linetype only)
mu_plot_df   <- do.call(rbind, mu_store)  # seed==1, i == 1
mu_plot_long <- melt(mu_plot_df, id.vars = c("index", "time", "n_ave"),
                     variable.name = "Method", value.name = "mu")
mu_plot_long$Method <- factor(mu_plot_long$Method,
                              levels = c("True", "FDA", "Traditional"))

p2_bw <- ggplot(mu_plot_long, aes(x = time, y = mu, linetype = Method)) +
  geom_line(color = "black") +
  facet_wrap(~ n_ave, ncol = 2, labeller = label_both) +
  labs(x = "Time", y = bquote(mu[i == .(1)](t))) +  # i=1
  scale_linetype_manual(name = "Legend",
                        values = c("True" = "solid", "FDA" = "longdash", "Traditional" = "dotdash")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

ggsave("mu_compare_traditional_bw.png", p2_bw, width = 8, height = 6, dpi = 300)


# 2b) NEW: Overlay plot of the same mu_i(t) with different n_ave (True fixed, FDA varies)
mu_overlay_long <- subset(mu_plot_long, Method %in% c("True", "FDA"))

# True should be identical across n_ave for the same i when seed==1; take one representative
true_once <- subset(mu_overlay_long, Method == "True")[, c("time", "mu")]
true_once <- true_once[!duplicated(true_once$time), ]

fda_df <- subset(mu_overlay_long, Method == "FDA")

# Stable linetypes by n_ave:
linetypes_by_n <- c("1" = "dotted", "30" = "dotdash", "70" = "longdash", "100" = "twodash")
fda_df$n_ave <- factor(fda_df$n_ave, levels = names(linetypes_by_n))

p2_overlay_bw <- ggplot() +
  geom_line(data = true_once, aes(x = time, y = mu),
            color = "black", linewidth = 1.05) +
  geom_line(data = fda_df, aes(x = time, y = mu, linetype = n_ave),
            color = "black") +
  scale_linetype_manual(name = "n_ave (FDA)", values = linetypes_by_n,
                        breaks = names(linetypes_by_n)) +
  labs(x = "Time", y = bquote(mu[i == .(1)](t))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

ggsave("mu_overlay_by_n_bw.png", p2_overlay_bw, width = 7, height = 5, dpi = 300)


# 3) Leading eigenfunctions of G_mu (linetype for True vs Estimated)
if (length(eigenfun_store) > 0) {
  eigenfun_df <- do.call(rbind, eigenfun_store)
  eigenfun_df$Method <- factor(eigenfun_df$Method, levels = c("True", "Estimated"))

  p3_bw <- ggplot(eigenfun_df, aes(x = time, y = value, linetype = Method)) +
    geom_line(color = "black") +
    facet_grid(component ~ n_ave, labeller = label_both) +
    labs(x = "Time", y = "Eigenfunction value") +
    scale_linetype_manual(name = "Legend",
                          values = c("True" = "solid", "Estimated" = "longdash")) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank())

  ggsave("G_mu_eigenfunctions_bw.png", p3_bw, width = 10, height = 6, dpi = 300)
}

# 4) Leading eigenvalues of G_mu (linetype + shape for True vs Estimated)
if (length(eigenval_store) > 0) {
  eigenval_df <- do.call(rbind, eigenval_store)
  eigenval_df$Method <- factor(eigenval_df$Method, levels = c("True", "Estimated"))

  p4_bw <- ggplot(eigenval_df, aes(x = component, y = value,
                                   linetype = Method, shape = Method)) +
    geom_line(color = "black") +
    geom_point(color = "black") +
    facet_wrap(~ n_ave, labeller = label_both) +
    scale_x_continuous(breaks = unique(eigenval_df$component)) +
    labs(x = "Component", y = "Eigenvalue") +
    scale_linetype_manual(name = "Legend",
                          values = c("True" = "solid", "Estimated" = "longdash")) +
    scale_shape_manual(name = "Legend",
                       values = c("True" = 16, "Estimated" = 17)) +  # solid circle / triangle
    theme_bw(base_size = 13) +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank())

  ggsave("G_mu_eigenvalues_bw.png", p4_bw, width = 8, height = 5, dpi = 300)
}

cat("BW plots saved:\n",
    "- rmse_convergence_bw.png\n",
    "- mu_compare_traditional_bw.png\n",
    "- mu_overlay_by_n_bw.png\n",
    "- G_mu_eigenfunctions_bw.png (if available)\n",
    "- G_mu_eigenvalues_bw.png (if available)\n", sep = "")
cat("\nData tables saved: convergence_summary.csv, G_mu_estimates_seed1.csv\n")
