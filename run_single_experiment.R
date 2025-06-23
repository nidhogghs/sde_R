args <- commandArgs(trailingOnly = TRUE)
K <- as.integer(args[1])
n_ave <- as.integer(args[2])
m <- as.integer(args[3])
seed_id <- as.integer(args[4])
Nsim <- as.integer(args[5])
l <- as.integer(args[6])
k <- as.integer(args[7])

cat(sprintf("===== BEGIN: K=%d, n_ave=%d, m=%d, l=%d, k=%d, seed_id=%d =====\n",
            K, n_ave, m, l, k, seed_id))

# 加载模块
source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/mu_sigma_est.R")

res_group <- data.frame()

for (s in 1:Nsim) {
  cat(sprintf("  > Simulation %d of %d (seed base: %d)\n", s, Nsim, seed_id))

  sim_data <- prepare_simulation_data(
    K = K, n_ave = n_ave, m = m,
    l = l, k = k, scale_sigma = 10,
    seed_sigma = seed_id * 100 + s,
    seed_mu = seed_id * 200 + s,
    seed_X = seed_id * 300 + s
  )

  cat("    - Data generated. Estimating mu...\n")
  mu_res <- estimate_mu_from_data(sim_data)

  cat("    - Estimating sigma²...\n")
  sigma_res <- inference_sigma_from_data(sim_data)

  rmse_mu <- sqrt(mean((mu_res$mu_hat - sim_data$mu_true)^2))
  rmse_sigma <- sqrt(mean((sigma_res$sigma2_hat - sim_data$sigma2_true)^2))

  cat(sprintf("    - RMSE(mu): %.6f, RMSE(sigma²): %.6f\n", rmse_mu, rmse_sigma))

  res_group <- rbind(res_group, data.frame(
    K = K, n_ave = n_ave, m = m,
    l = l, k = k,
    seed = seed_id,
    sim = s,
    rmse_mu = rmse_mu,
    rmse_sigma = rmse_sigma
  ))
}

outfile <- sprintf("results/K%d_n%d_m%d_l%d_k%d_seed%d.rds", K, n_ave, m, l, k, seed_id)
saveRDS(res_group, file = outfile)

cat(sprintf("===== DONE: Results saved to %s =====\n", outfile))
