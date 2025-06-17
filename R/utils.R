# R/utils.R
# ========== 工具函数集合 ==========

# 打印变量名和其值，用于调试
print_var <- function(var) {
  var_name <- deparse(substitute(var))  # 获取变量的名称作为字符串
  cat(var_name, "=", var, "\n")         # 打印“变量名 = 值”
}

# 按 cluster（每个个体）对观测进行平均，得到 K 条路径的平均值（用于 Z_Δ 聚合）
# 输入:
#   X: 所有观测轨迹组成的矩阵，列为多个观测点（总长度 = sum(n)）
#   K: 个体数量
#   n: 每个个体的观测次数向量，长度为 K
# 输出:
#   X_bar: 每个个体的平均轨迹，维度为 m × K
cluster_mean <- function(X, K, n) {
  n <- c(0, n)  # 插入初始 0，方便用 sum(n[1:k]) 定位区间起点
  X_bar <- matrix(NA, nrow = nrow(X), ncol = K)  # 初始化输出矩阵

  for (k in 1:K) {
    start_idx <- sum(n[1:k]) + 1
    end_idx <- sum(n[1:(k + 1)])
    X_bar[, k] <- rowMeans(X[, start_idx:end_idx])  # 对每个 cluster 求平均轨迹
  }

  return(X_bar)
}

# 基于 AIC 准则选择最优主成分数 L（截断维数）
# 输入:
#   Z_Delta: 标准化后的观测矩阵（维度 m × K）
#   m_mu: 已估计出的漂移均值函数（长度 m）
#   PCs: 主成分基向量（m × L_max 维，来自 G_mu 的特征向量）
#   max_L: 最大允许的主成分数
# 输出:
#   L_hat: 使 AIC 最小的主成分数
select_L_by_AIC <- function(Z_Delta, m_mu, PCs, max_L = 10) {
  K <- ncol(Z_Delta)  # 样本数（个体数）
  m <- nrow(Z_Delta)  # 时间点数
  aic_values <- numeric(max_L)  # 储存各个 L 值对应的 AIC

  for (L in 1:max_L) {
    Z_hat <- matrix(0, m, K)  # 重构的 Z 值
    for (i in 1:K) {
      Z_hat[, i] <- m_mu  # 初始化为均值函数
      for (j in 1:L) {
        # 主成分系数估计：投影
        xi_hat <- mean((Z_Delta[, i] - m_mu) * PCs[, j])
        # 逐个主成分叠加重构
        Z_hat[, i] <- Z_hat[, i] + xi_hat * PCs[, j]
      }
    }

    # 残差平方和
    rss <- sum((Z_Delta - Z_hat)^2)

    # AIC 计算：拟合误差 + 模型复杂度惩罚
    aic_values[L] <- K * log(rss / (K * m)) + 2 * L * m
  }

  # 返回使 AIC 最小的 L
  return(which.min(aic_values))
}



plot_compare_single <- function(ts, true_mat, est_mat, k, label = "V") {
  if (is.vector(true_mat)) true_mat <- matrix(true_mat, ncol = 1)
  if (is.vector(est_mat))  est_mat  <- matrix(est_mat,  ncol = 1)

  df <- data.frame(
    t = ts,
    True = true_mat[, k],
    Estimate = est_mat[, k]
  )
  df_melt <- reshape2::melt(df, id.vars = "t", variable.name = "Type", value.name = "Value")

  p <- ggplot(df_melt, aes(x = t, y = Value, color = Type)) +
    geom_line(linewidth = 1) +
    labs(title = sprintf("Comparison of %s: True vs Estimated (k = %d)", label, k),
         x = "t", y = label) +
    theme_minimal()

  print(p)  # 显示图像
}

# 计算 q0 的值
compute_qk_mc <- function(nk, B = 1e6, seed = 123) {
  # 设置随机种子以保证可复现
  set.seed(seed)
  
  # 生成 B 个自由度为 nk 的卡方随机变量
  samples <- rchisq(B, df = nk)
  
  # 计算期望 E[log(χ²_nk)] - log(nk)
  qk <- mean(log(samples)) - log(nk)
  
  return(qk)
}

evaluate_V_estimation <- function(V_true, V_hat, sigma2_true, sigma2_hat, k_eval = 1) {
  stopifnot(dim(V_true) == dim(V_hat))
  stopifnot(dim(sigma2_true) == dim(sigma2_hat))
  stopifnot(k_eval >= 1 && k_eval <= ncol(V_true))

  # 所有轨迹的 log σ² 误差
  err_V_col_max <- apply(abs(V_hat - V_true), 2, max)
  err_V_max <- max(err_V_col_max)
  err_V_mean <- mean(err_V_col_max)

  # 第 k_eval 条轨迹：V 和 σ² 的误差
  V_diff_k <- V_hat[, k_eval] - V_true[, k_eval]
  sigma2_diff_k <- sigma2_hat[, k_eval] - sigma2_true[, k_eval]

  # 相对误差（已归一化）
  err_V_k_rel <- max(abs(V_diff_k)) / (max(V_true[, k_eval]) - min(V_true[, k_eval]))
  err_sigma2_k_rel <- sqrt(mean(sigma2_diff_k^2)) / mean(sigma2_true[, k_eval])

  # 绝对误差（MAE / RMSE）
  err_sigma2_k_mae <- mean(abs(sigma2_diff_k))
  err_sigma2_k_rmse <- sqrt(mean(sigma2_diff_k^2))

  # 输出
  cat("===== Estimation Errors =====\n")
  cat(sprintf("Max abs error over all V_k(t):     %0.6f\n", err_V_max))
  cat(sprintf("Mean abs error over all V_k(t):    %0.6f\n", err_V_mean))
  cat(sprintf("Relative error of V_k(t), k=%d:     %0.6f\n", k_eval, err_V_k_rel))
  cat(sprintf("Relative error of sigma²_k(t), k=%d:    %0.6f\n", k_eval, err_sigma2_k_rel))
  cat(sprintf("Absolute MAE  of sigma²_k(t), k=%d:     %0.6f\n", k_eval, err_sigma2_k_mae))
  cat(sprintf("Absolute RMSE of sigma²_k(t), k=%d:     %0.6f\n", k_eval, err_sigma2_k_rmse))

  # 返回结果
  return(list(
    err_V_max = err_V_max,
    err_V_mean = err_V_mean,
    err_V_k_rel = err_V_k_rel,
    err_sigma2_k_rel = err_sigma2_k_rel,
    err_sigma2_k_mae = err_sigma2_k_mae,
    err_sigma2_k_rmse = err_sigma2_k_rmse
  ))
}

