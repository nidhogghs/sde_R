# test/test_inference_sigma2_large.R
source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
library(np)

# Step 1: 参数设置
K <- 500
n_ave <- 300
m <- m_default
delta <- T / m
ts <- seq(delta, T, length.out = m)
n_vec <- rep(n_ave, K)

# Step 2: 模拟 sigma, mu 并规范化
set.seed(456)
sigma <- generate_K_trajectory(K, a, delta, alpha, k, m)
sigma <- sigma / max(sigma)  # 归一化防止数值不稳定

set.seed(789)
mu <- generate_K_trajectory(K, b, delta, beta, l, m)

# Step 3: 生成 X, 并构造 Z_Delta 与每类平均
set.seed(1)
X <- generate_n_X(n_vec, sigma, mu, X0)
logX <- log(X)
X_Delta <- variation(logX, m)                   # m x (K*n_ave)
Z_Delta <- cluster_mean(X_Delta, K, n_vec) / sqrt(delta)  # m x K

# Step 4: 提取第1类观测的“原始差分”
Z_Delta_raw <- X_Delta * sqrt(delta)            # ΔlogX
Z1 <- Z_Delta_raw[, 1:n_ave]                    # 第1类

# Step 5: 用 Müller 方法估计 σ²_k(t)
sigma2_hat_k <- estimate_sigma2_muller(Z1, ts)

# Step 6: 真值（注意对齐差分）
sigma2_true_k <- sigma[2:(m+1), 1]^2

# Step 7: 绘图 + RMSE
plot(ts, sigma2_hat_k, type = "l", col = "blue", lwd = 2,
     ylim = range(c(sigma2_hat_k, sigma2_true_k)),
     main = expression(hat(sigma)^2(t)~vs.~sigma[true]^2(t)),
     xlab = "t", ylab = expression(sigma^2(t)))
lines(ts, sigma2_true_k, col = "black", lty = 3, lwd = 2)
legend("topright", legend = c("估计值", "真实值"),
       col = c("blue", "black"), lty = c(1, 3), lwd = 2)

rmse <- sqrt(mean((sigma2_hat_k - sigma2_true_k)^2))
cat("RMSE for σ²_k(t) estimation (k=1):", rmse, "\n")
