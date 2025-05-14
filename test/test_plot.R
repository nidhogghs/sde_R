# test_plot.R
# ========= 可视化测试（ASCII安全） ==========

source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/one_run.R")

library(ggplot2)

# === 参数设置 ===
K <- 100
n_ave <- 100
m <- 50
sd <- 42

delta <- T / m
ts <- seq(delta, 1, delta)

# === 真值生成 ===
set.seed(456); sigma <- generate_K_trajectory(K, a, delta, alpha, k, m); sigma <- sigma / max(sigma)
set.seed(789); mu <- generate_K_trajectory(K, b, delta, beta, l, m)
m_mu_true <- rowMeans(mu)[-1]
G_mu_true <- (mu[-1,] - m_mu_true) %*% t(mu[-1,] - m_mu_true) / K
PCs_true <- Re(eigen(G_mu_true)$vectors) * sqrt(m)
lams_true <- Re(eigen(G_mu_true)$values) / m

# === 模拟观测 ===
set.seed(sd)
n_vec <- rep(n_ave, K)
X <- generate_n_X(n_vec, sigma, mu, X0)
X_Delta <- variation(log(X), m)
Z_Delta <- cluster_mean(X_Delta, K, n_vec) / sqrt(delta)

# === 估计 ===
m_mu_hat <- estimate_m_mu(Z_Delta, ts)
G_mu_hat <- estimate_G_mu(Z_Delta, m_mu_hat, m, K, h_min = h_min)
PCs_hat <- Re(eigen(G_mu_hat)$vectors) * sqrt(m)
lams_hat <- Re(eigen(G_mu_hat)$values) / m
L <- select_L_by_AIC(Z_Delta, m_mu_hat, PCs_hat, max_L = 10)
mu_hat <- recover_mu(Z_Delta, m_mu_hat, PCs_hat, L, K)

# === 图像保存 ===
dir.create("figures", showWarnings = FALSE)

# 1. m_mu(t) 对比
df_m <- data.frame(t = ts, True = m_mu_true, Estimated = m_mu_hat)
p1 <- ggplot(df_m, aes(x = t)) +
  geom_line(aes(y = True, color = "True")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "Mean Function m_mu(t)", y = "Value") +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  theme_minimal()
ggsave("figures/plot_m_mu.pdf", p1, width = 6, height = 4)

# 2. 第一主成分函数 phi_1(t)
df_pc1 <- data.frame(
  t = ts,
  True = PCs_true[,1] * sign(PCs_true[1,1]),
  Estimated = PCs_hat[,1] * sign(PCs_hat[1,1])
)
p2 <- ggplot(df_pc1, aes(x = t)) +
  geom_line(aes(y = True, color = "True")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "First Principal Component phi_1(t)", y = "Value") +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  theme_minimal()
ggsave("figures/plot_pc1.pdf", p2, width = 6, height = 4)

# 3. 第1个 mu_i(t) 重构效果
df_mu1 <- data.frame(t = ts, True = mu[-1,1], Estimated = mu_hat[,1])
p3 <- ggplot(df_mu1, aes(x = t)) +
  geom_line(aes(y = True, color = "True")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "First Sample mu_1(t)", y = "Value") +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  theme_minimal()
ggsave("figures/plot_mu1.pdf", p3, width = 6, height = 4)

# 4. 特征值谱 lambda_j
df_lam <- data.frame(Index = 1:length(lams_true), True = lams_true, Estimated = lams_hat)
p4 <- ggplot(df_lam, aes(x = Index)) +
  geom_point(aes(y = True, color = "True")) +
  geom_point(aes(y = Estimated, color = "Estimated")) +
  labs(title = "Eigenvalue Spectrum lambda_j", y = "Value") +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  theme_minimal()
ggsave("figures/plot_eigenvalues.pdf", p4, width = 6, height = 4)
