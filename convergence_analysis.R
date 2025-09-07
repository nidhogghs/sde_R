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

# ---- Parse args ----
args <- commandArgs(trailingOnly = TRUE)
num_seeds <- 100
for (arg in args) {
  if (grepl("^num_seeds=", arg)) {
    num_seeds <- as.integer(sub("num_seeds=", "", arg))
  }
}

K <- 500
n_values <- c(1, 30, 70, 100)

results <- data.frame()
mu_store <- list()

for (n in n_values) {
  rmse_list <- numeric(num_seeds)
  for (seed in seq_len(num_seeds)) {
    set.seed(seed)
    sim_data <- prepare_simulation_data(
      K = K,
      n_ave = n,
      m = 50,
      k = 2,
      l = 2,
      scale_sigma = 5,
      seed_sigma = 311 + seed,
      seed_mu = 311 + seed,
      seed_X = 69 + seed
    )

    mu_res <- estimate_mu_from_data(sim_data)
    sigma_res <- inference_sigma_from_data(sim_data)
    rmse_mu <- sqrt(mean((mu_res$mu_hat - sim_data$mu_true)^2))
    rmse_list[seed] <- rmse_mu

    if (seed == 1) {
      mu_store[[as.character(n)]] <- data.frame(
        index = seq_len(nrow(mu_res$mu_hat)),
        mu_hat = mu_res$mu_hat[, 1],
        mu_true = sim_data$mu_true[, 1],
        n_ave = n
      )
    }
  }
  results <- rbind(results, data.frame(n_ave = n, avg_rmse_mu = mean(rmse_list)))
}

# ---- Visualization ----
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

p1 <- ggplot(results, aes(x = n_ave, y = avg_rmse_mu)) +
  geom_line() + geom_point(size = 2) +
  labs(title = sprintf("RMSE of mu vs n_ave (K=%d)", K),
       x = "n_ave", y = "Average RMSE(mu)") +
  theme_minimal()

ggsave("rmse_convergence.png", p1, width = 5, height = 4)

mu_plot_df <- do.call(rbind, mu_store)

p2 <- ggplot(mu_plot_df, aes(x = index)) +
  geom_line(aes(y = mu_true), colour = "black", linetype = "dashed") +
  geom_line(aes(y = mu_hat, colour = factor(n_ave))) +
  facet_wrap(~ n_ave, ncol = 2) +
  labs(title = sprintf("Estimated mu vs True mu (K=%d, seed=1)", K),
       x = "Time index", y = "mu", colour = "n_ave") +
  theme_minimal()

ggsave("mu_approach_true.png", p2, width = 8, height = 6)

cat("Simulation complete. Plots saved to rmse_convergence.png and mu_approach_true.png\n")
