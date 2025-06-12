source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
library(np)

cat("===== TEST SCRIPT START: Estimating sigma²_k(t) =====\n")

# Step 1: Read parameters (from CLI or default) --------------------------
args <- commandArgs(trailingOnly = TRUE)

# Default values
K <- 1
n_ave <- 300
m <- 50

# Parse CLI arguments like "K=50"
for (arg in args) {
  eval(parse(text = arg))
}
cat(sprintf("[Step 1] Parameter setup: K = %d, n_ave = %d, m = %d\n", K, n_ave, m))

delta <- T / m
ts <- seq(delta, T, length.out = m)
n_vec <- rep(n_ave, K)

# Step 2: Simulate sigma and mu ------------------------------------------
cat("[Step 2] Simulating sigma_k(t) and mu_k(t)...\n")
set.seed(456)
sigma <- generate_K_trajectory(K, a, delta, alpha, k, m)
sigma <- sigma / max(sigma)  # Normalize

set.seed(789)
mu <- generate_K_trajectory(K, b, delta, beta, l, m)
cat(sprintf("[Step 2] Done: dim(sigma) = [%d x %d], dim(mu) = [%d x %d]\n",
            nrow(sigma), ncol(sigma), nrow(mu), ncol(mu)))

# Step 3: Simulate X and compute log-differences --------------------------
cat("[Step 3] Generating X and computing Δlog X...\n")
set.seed(1)
X <- generate_n_X(n_vec, sigma, mu, X0)
logX <- log(X)
Z_Delta_raw <- variation(logX, m) * sqrt(delta)
cat(sprintf("[Step 3] Done: dim(X) = [%d x %d], dim(Z_Delta) = [%d x %d]\n",
            nrow(X), ncol(X), nrow(Z_Delta_raw), ncol(Z_Delta_raw)))

# Step 4: Extract data for class k = 1 ------------------------------------
cat("[Step 4] Extracting data for class k = 1...\n")
Z_k <- Z_Delta_raw[, 1:n_ave]
cat(sprintf("[Step 4] Done: dim(Z_k) = [%d x %d]\n", nrow(Z_k), ncol(Z_k)))

# Step 5: Estimate sigma²_k(t) --------------------------------------------
cat("[Step 5] Estimating sigma²_k(t)...\n")
sigma2_hat_k <- estimate_sigma2_k(Z_k, ts, L = 3)
cat(sprintf("[Step 5] Done: range = [%.4f, %.4f]\n", 
            min(sigma2_hat_k), max(sigma2_hat_k)))

# Step 6: Get true sigma²_k(t) --------------------------------------------
sigma2_true_k <- sigma[2:(m+1), 1]^2
cat("[Step 6] True sigma²_k(t) computed.\n")

# Step 7: Plotting --------------------------------------------------------
cat("[Step 7] Plotting estimated vs. true sigma²_k(t)...\n")
plot(ts, sigma2_hat_k, type = "l", col = "blue", lwd = 2,
     ylim = range(c(sigma2_hat_k, sigma2_true_k)),
     main = expression(hat(sigma)^2(t)~vs.~sigma[true]^2(t)),
     xlab = "t", ylab = expression(sigma^2(t)))
lines(ts, sigma2_true_k, col = "black", lty = 3, lwd = 2)
legend("topright", legend = c("Estimated", "True"),
       col = c("blue", "black"), lty = c(1, 3), lwd = 2)
cat("[Step 7] Plot finished.\n")

# Step 8: Compute RMSE ----------------------------------------------------
rmse <- sqrt(mean((sigma2_hat_k - sigma2_true_k)^2))
cat(sprintf("[Step 8] RMSE = %.6f\n", rmse))
cat("===== TEST SCRIPT END =====\n")

