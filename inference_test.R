source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
library(np)
library(ggplot2)
library(reshape2)

cat("===== TEST SCRIPT START: Estimating sigma²_k(t) =====\n")

# 1: Read parameters (from CLI or default) --------------------------
args <- commandArgs(trailingOnly = TRUE)

# Default values
K <- 300
n_ave <- 300
m <- 50
L <- 4
for (arg in args) {
  eval(parse(text = arg))
}


delta <- T / m
ts <- seq(delta, T, length.out = m)  #时间间隔网
n_vec <- rep(n_ave, K) #K个过程的观测数为n_ave

compute_qk <- function(nk) {
  - 1 / (2 * nk) - 1 / (12 * nk^2)
}


q0 <- compute_qk(n_ave)
# cat(sprintf("q0 = %f\n", q0))


# q0 <- -1.270  # 预设 q0 值，通常为 E[log(W²)] 的估计值
# 2. 真值生成 ------------------------------------------


set.seed(456)
sigma <- generate_K_trajectory(K, a, delta, alpha, k, m)#m+1* K 矩阵
sigma <- sigma / max(sigma)  # 确保最大值为 1
sigma2 <- sigma^2 # m × K 矩阵，所有curve的sigma²_k(t)
V_true <- log(sigma2) # m × K 矩阵，所有curve的V_k(t)


m_V_true_value <- rowMeans(V_true)[-1]        # m × 1 向量值为估计值
G_V_true_value <- (V_true[-1, ] - m_V_true_value) %*% t(V_true[-1, ] - m_V_true_value) / K # m × m 矩阵
PCs_true <- Re(eigen(G_V_true_value)$vectors) * sqrt(m) # m × m 矩阵.主成分离散点值
lams_true <- Re(eigen(G_V_true_value)$values) / m # m × 1 向量.主成分贡献


sigma2_true <- exp(V_true)[-1, ] # m × K 矩阵，所有curve的sigma²_k(t)
sigma_true <- sqrt(sigma2_true) # m × K 矩阵，所有curve的sigma_k(t)


# 生成漂移轨迹
set.seed(789)
mu <- generate_K_trajectory(K, b, delta, beta, l, m)

set.seed(1)
X <- generate_n_X(n_vec, sigma, mu, X0)

# 3. 标准化差分数据
X_Delta <- variation(log(X), m) #差分 得到m*sum(n_vec)矩阵
Y_Delta <-construct_Y(X_Delta, n_vec, q0)  # m × K 矩阵

# 4. 估计 m_mu 与 G_mu
m_V_hat <- estimate_m_mu(Y_Delta, ts) # m × 1 向量值为估计值
G_V_hat <- estimate_G_mu(Y_Delta, m_V_hat, m, K, h_min = h_min)  # m × m 矩阵

# 5. 主成分分析 (FPCA)
PCs_hat <- Re(eigen(G_V_hat)$vectors) * sqrt(m)
lams_hat <- Re(eigen(G_V_hat)$values) / m

# 6.自动选择主成分数 L（若未指定）
if (is.null(L)) {
  L <- select_L_by_AIC(Y_Delta, m_V_hat, PCs_hat, max_L = 10)
  if (is.na(L)) stop("AIC failed: selected L is NA.")
  cat(sprintf("Selected L by AIC: %d\n", L))
} else {
  cat(sprintf("Using manually specified L: %d\n", L))
}

V_hat <- recover_mu(Y_Delta, m_V_hat, PCs_hat, L, K)

sigma2_hat <- exp(V_hat) # m × K 矩阵，所有curve的sigma²_k(t)
sigma_hat <- sqrt(sigma2_hat) # m × K 矩阵，所有curve的sigma_k(t)




p1 <- plot_compare_single(ts, m_V_true_value, V_hat, k = 1, label = "V_k(t)")
ggsave("figures/V_k1.png", plot = p1, width = 6, height = 4)

p2 <- plot_compare_single(ts, sigma2_true, sigma2_hat, k = 1, label = "sigma²_k(t)")
ggsave("figures/sigma2_k1.png", plot = p2, width = 6, height = 4)


