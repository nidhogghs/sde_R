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
K <- 300
n_ave <- 300
m <- 50

# Parse CLI arguments like "K=50"
for (arg in args) {
  eval(parse(text = arg))
}


delta <- T / m
ts <- seq(delta, T, length.out = m)
n_vec <- rep(n_ave, K)

# Step 2: Simulate sigma and mu ------------------------------------------



set.seed(456)
sigma_res <- generate_V_and_sigma2(K, a, delta, alpha, k, m)
sigma  <- sigma_res$sigma        # 已归一化
sigma2 <- sigma_res$sigma2
V_true <- sigma_res$V

set.seed(789)
mu <- generate_K_trajectory(K, b, delta, beta, l, m)



# Step 3: Simulate X and compute log-differences --------------------------

set.seed(1)
X <- generate_n_X(n_vec, sigma, mu, X0)
logX <- log(X)
X_Delta <- variation(logX, m)




# Step 4: Extract data for class k = 1 ------------------------------------

Z_Delta <- cluster_mean(X_Delta, K, n_vec) / sqrt(delta)  # dim = [m × K]
epsilon <- 1e-8
Z_Delta[abs(Z_Delta) < epsilon | is.na(Z_Delta) | is.infinite(Z_Delta)] <- epsilon

q0 <- -1.270
Y_all <- log(Z_Delta^2) - q0
Y_k <- Y_all[, 1, drop = FALSE]  # 提取第 1 个 process，确保为矩阵形式

# Step 5: Estimate sigma²_k(t) --------------------------------------------


sigma2_hat_k <- estimate_sigma2_Muller(Y_all, Y_k, ts, L = 3)

# Step 6: Get true sigma²_k(t) from V_true
V_true_k <- V_true[2:(m+1), 1]           # log σ²_k(t)
sigma2_true_k <- exp(V_true_k)

# Step 7: Plot and save to file
png("output/sigma2_est_vs_true.png", width = 800, height = 600)
plot(ts, sigma2_hat_k, type = "l", col = "blue", lwd = 2,
     ylim = range(c(sigma2_hat_k, sigma2_true_k)),
     main = expression(hat(sigma)^2(t)~vs.~sigma[true]^2(t)),
     xlab = "t", ylab = expression(sigma^2(t)))
lines(ts, sigma2_true_k, col = "black", lty = 3, lwd = 2)
legend("topright", legend = c("Estimate", "True"),
       col = c("blue", "black"), lty = c(1, 3), lwd = 2)
dev.off()


# Step 8: Compute RMSE ----------------------------------------------------
rmse <- sqrt(mean((sigma2_hat_k - sigma2_true_k)^2))
cat(sprintf("[Step 8] RMSE = %.6f\n", rmse))
cat("===== TEST SCRIPT END =====\n")

