# =========================================
# convergence_analysis.R
# -----------------------------------------
# Run simulations for K = 500 with varying n_ave to study convergence.
# Generates RMSE plot and mu comparison plots.
# Usage:
#   Rscript convergence_analysis.R [num_seeds=100]
# Example (quick test with fewer seeds):
#   Rscript convergence_analysis.R num_seeds=2
# =========================================

# ---- Load modules ----
source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/mu_sigma_est.R")
source("R/traditional.R")

# ---- Parse args ----
args <- commandArgs(trailingOnly = TRUE)
num_seeds <- 1
for (arg in args) {
  if (grepl("^num_seeds=", arg)) {
    num_seeds <- as.integer(sub("num_seeds=", "", arg))
  }
}

K <- 500
n_values <- c(1, 30, 70, 100)

results <- data.frame()
mu_store <- list()
G_mu_store <- list()
eigenfun_store <- list()
eigenval_store <- list()

for (n in n_values) {
  rmse_fda_list <- numeric(num_seeds)
  rmse_trad_list <- numeric(num_seeds)
  rmse_G_list <- numeric(num_seeds)
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

    mu_res <- estimate_mu_from_data(sim_data)
    mu_trad <- estimate_mu_traditional_from_data(sim_data)
    sigma_res <- inference_sigma_from_data(sim_data)

    rmse_fda <- sqrt(mean((mu_res$mu_hat - sim_data$mu_true)^2))
    rmse_trad <- sqrt(mean((mu_trad - sim_data$mu_true)^2))
    rmse_fda_list[seed] <- rmse_fda
    rmse_trad_list[seed] <- rmse_trad

    G_mu_true <- (sim_data$mu_true - rowMeans(sim_data$mu_true)) %*%
      t(sim_data$mu_true - rowMeans(sim_data$mu_true)) / K
    rmse_G <- sqrt(mean((mu_res$G_mu_hat - G_mu_true)^2))
    rmse_G_list[seed] <- rmse_G

    if (seed == 1) {
      ts_vals <- sim_data$ts
      mu_store[[as.character(n)]] <- data.frame(
        index = seq_len(nrow(mu_res$mu_hat)),
        time = ts_vals,
        True = sim_data$mu_true[, 1],
        FDA = mu_res$mu_hat[, 1],
        Traditional = mu_trad[, 1],
        n_ave = n
      )

      m_dim <- nrow(mu_res$G_mu_hat)
      G_mu_store[[as.character(n)]] <- data.frame(
        i = rep(seq_len(m_dim), times = m_dim),
        j = rep(seq_len(m_dim), each = m_dim),
        G_mu_hat = as.vector(mu_res$G_mu_hat),
        G_mu_true = as.vector(G_mu_true),
        n_ave = n
      )

      eig_true <- eigen(G_mu_true)
      PCs_true <- Re(eig_true$vectors) * sqrt(m_dim)
      lams_true <- Re(eig_true$values) / m_dim

      PCs_hat <- mu_res$PCs_hat
      lams_hat <- mu_res$lams_hat

      num_components <- min(3, ncol(PCs_hat))
      if (num_components > 0) {
        PCs_true_sel <- PCs_true[, seq_len(num_components), drop = FALSE]
        PCs_hat_sel <- PCs_hat[, seq_len(num_components), drop = FALSE]

        for (comp in seq_len(num_components)) {
          align_sign <- sign(sum(PCs_hat_sel[, comp] * PCs_true_sel[, comp]))
          if (align_sign == 0) {
            align_sign <- 1
          }
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
    avg_rmse_mu_fda = mean(rmse_fda_list),
    avg_rmse_mu_traditional = mean(rmse_trad_list),
    avg_rmse_G_mu = mean(rmse_G_list)
  ))
}

# ---- Save numeric tables ----
write.csv(results, "convergence_summary.csv", row.names = FALSE)
G_mu_df <- do.call(rbind, G_mu_store)
write.csv(G_mu_df, "G_mu_estimates_seed1.csv", row.names = FALSE)

# ---- Visualization ----
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}
library(ggplot2)
library(reshape2)

results_long <- melt(results, id.vars = "n_ave",
                     measure.vars = c("avg_rmse_mu_fda", "avg_rmse_mu_traditional"),
                     variable.name = "Method", value.name = "avg_rmse_mu")
results_long$Method <- factor(results_long$Method,
                             levels = c("avg_rmse_mu_fda", "avg_rmse_mu_traditional"),
                             labels = c("FDA", "Traditional"))

p1 <- ggplot(results_long, aes(x = n_ave, y = avg_rmse_mu, color = Method)) +
  geom_line() + geom_point(size = 2) +
  labs(title = sprintf("RMSE of mu vs n_ave (N=%d)", K),
       x = "n_ave", y = "Average RMSE(mu)", color = "Method") +
  theme_minimal()

ggsave("rmse_convergence.png", p1, width = 5, height = 4)

mu_plot_df <- do.call(rbind, mu_store)
mu_plot_long <- melt(mu_plot_df, id.vars = c("index", "time", "n_ave"),
                     variable.name = "Method", value.name = "mu")

p2 <- ggplot(mu_plot_long, aes(x = time, y = mu, color = Method)) +
  geom_line() +
  facet_wrap(~ n_ave, ncol = 2) +
  labs(title = sprintf("True vs Estimates (N=%d, seed=1, n=1)", K),
       x = "Time", y = "mu") +
  theme_minimal()

ggsave("mu_compare_traditional.png", p2, width = 8, height = 6)

if (length(eigenfun_store) > 0) {
  eigenfun_df <- do.call(rbind, eigenfun_store)
  p3 <- ggplot(eigenfun_df, aes(x = time, y = value, color = Method, linetype = Method)) +
    geom_line() +
    facet_grid(component ~ n_ave, labeller = label_both) +
    labs(title = sprintf("Leading eigenfunctions of G_mu (N=%d, seed=1)", K),
         x = "Time", y = "Eigenfunction value") +
    theme_minimal()

  ggsave("G_mu_eigenfunctions.png", p3, width = 10, height = 6)
}

if (length(eigenval_store) > 0) {
  eigenval_df <- do.call(rbind, eigenval_store)
  p4 <- ggplot(eigenval_df, aes(x = component, y = value, color = Method)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ n_ave, labeller = label_both) +
    scale_x_continuous(breaks = unique(eigenval_df$component)) +
    labs(title = sprintf("Leading eigenvalues of G_mu (N=%d, seed=1)", K),
         x = "Component", y = "Eigenvalue") +
    theme_minimal()

  ggsave("G_mu_eigenvalues.png", p4, width = 8, height = 5)
}

cat("Simulation complete. Plots saved to rmse_convergence.png, mu_compare_traditional.png, G_mu_eigenfunctions.png, and G_mu_eigenvalues.png\n")
cat("Data tables saved to convergence_summary.csv and G_mu_estimates_seed1.csv\n")