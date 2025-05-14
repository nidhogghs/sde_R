# R/one_run.R
# ========== 单次仿真流程主函数 ==========


one_run <- function(K, n_ave, m, sd, L = NULL) {
  # 设置基本时间参数
  delta <- T / m
  ts <- seq(delta, 1, delta)

  # 1. 生成漂移和扩散轨迹
  set.seed(456)
  sigma <- generate_K_trajectory(K, a, delta, alpha, k, m)
  sigma <- sigma / max(sigma)  # 归一化防止数值不稳定

  set.seed(789)
  mu <- generate_K_trajectory(K, b, delta, beta, l, m)

  # 真值（不带 t=0）
  m_mu_true_val <- rowMeans(mu)[-1]
  G_mu_true_val <- (mu[-1,] - m_mu_true_val) %*% t(mu[-1,] - m_mu_true_val) / K
  PCs_true <- Re(eigen(G_mu_true_val)$vectors) * sqrt(m)
  lams_true <- Re(eigen(G_mu_true_val)$values) / m

  # 2. 生成观测轨迹 X_i(t)
  set.seed(sd)
  n_vec <- rep(n_ave, K)
  X <- generate_n_X(n_vec, sigma, mu, X0)

  # 3. 标准化差分数据 Z_Δ
  X_Delta <- variation(log(X), m)
  Z_Delta <- if (max(n_vec) == 1) {
    X_Delta
  } else {
    cluster_mean(X_Delta, K, n_vec) / sqrt(delta)
  }

  # 4. 估计 m_mu 与 G_mu
  m_mu_hat <- estimate_m_mu(Z_Delta, ts)
  G_mu_hat <- estimate_G_mu(Z_Delta, m_mu_hat, m, K, h_min = h_min)

# 5. 主成分分析 (FPCA)
PCs_hat <- Re(eigen(G_mu_hat)$vectors) * sqrt(m)
lams_hat <- Re(eigen(G_mu_hat)$values) / m

# 自动选择主成分数 L（若未指定）
if (is.null(L)) {
  L <- select_L_by_AIC(Z_Delta, m_mu_hat, PCs_hat, max_L = 10)
  if (is.na(L)) stop("AIC failed: selected L is NA.")
  cat(sprintf("Selected L by AIC: %d\n", L))
} else {
  cat(sprintf("Using manually specified L: %d\n", L))
}

# 6. 估计每条 μ_i(t)
mu_hat <- recover_mu(Z_Delta, m_mu_hat, PCs_hat, L, K)
  # 7. 误差评估
  mu <- mu[-1,]  # 去掉 t=0 行

  err_m <- max(abs(m_mu_hat - m_mu_true_val))
  err_G <- max(abs(G_mu_hat - G_mu_true_val))
  err_lams <- max(abs(lams_hat[1:length(lams_true)] - lams_true))

  err_PCs <- max(sapply(1:l, function(i) {
    true <- PCs_true[, i] * sign(PCs_true[1, i])
    est <- PCs_hat[, i] * sign(PCs_hat[1, i])
    max(abs(true - est))
  }))

  err_mu <- max(abs(mu_hat[, 1] - mu[, 1])) / (max(mu[, 1]) - min(mu[, 1]))

  # 打印结果
  cat(sprintf("##### K = %d, n_ave = %d, m = %d #####\n", K, n_ave, m))
  print_var(err_m)
  print_var(err_G)
  print_var(err_lams)
  print_var(err_PCs)
  print_var(err_mu)

  return(c(err_m, err_G, err_lams, err_PCs, err_mu))

}
