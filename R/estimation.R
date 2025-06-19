# R/estimation.R
# ========== 核估计与FPCA ==========

library("np")
library("plotly")
estimate_m_mu <- function(Z_Delta, ts){
  data1 <- data.frame(Y = rowMeans(Z_Delta), X = ts)
  fit1 <- npreg(Y ~ X, data = data1, regtype = "ll", bwmethod = "cv.aic")
  predict(fit1, newdata = data.frame(X = ts))
}

estimate_G_mu <- function(Z_Delta, m_mu_hat, m, K, h_min = 0.021){
  coor1 <- rep(1:m, m) / m
  coor2 <- rep(1:m, each = m) / m

  sample_cov <- (Z_Delta - m_mu_hat) %*% t(Z_Delta - m_mu_hat) / K
  idx_diag <- which(coor1 == coor2)
  data2 <- data.frame(
    Y = as.vector(sample_cov)[-idx_diag],
    X1 = coor1[-idx_diag],
    X2 = coor2[-idx_diag]
  )

  h1 <- max((K * ncol(Z_Delta) / m)^(-1/3), h_min)
  fit2 <- npreg(Y ~ X1 + X2, data = data2, regtype = "ll", bws = c(h1, h1) * 0.2)
  G_hat <- predict(fit2, newdata = data.frame(X1 = coor1, X2 = coor2))
  return(matrix(G_hat, m, m))
}

# 构造 Y_Delta，支持每个过程单独的 q_k
construct_Y <- function(Z_Delta, n_vec, q_vec) {
  m <- nrow(Z_Delta)
  K <- length(n_vec)
  Y_Delta <- matrix(NA, nrow = m, ncol = K)

  col_idx <- 1
  for (k in 1:K) {
    idx_range <- col_idx:(col_idx + n_vec[k] - 1)
    Z_k <- Z_Delta[, idx_range, drop = FALSE]

    Z_sq_sum <- rowSums(Z_k^2)
    Y_k <- log(Z_sq_sum) - log(n_vec[k]) - q_vec[k]  # 用第 k 个 q_k

    Y_Delta[, k] <- Y_k
    col_idx <- col_idx + n_vec[k]
  }

  return(Y_Delta)  # 返回 m × K 的矩阵
}
