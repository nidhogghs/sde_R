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

construct_Y <- function(Z_Delta, n_vec, q0) {
  m <- nrow(Z_Delta)
  K <- length(n_vec)
  Y_all <- matrix(NA, nrow = m, ncol = K)

  col_idx <- 1
  for (k in 1:K) {
    idx_range <- col_idx:(col_idx + n_vec[k] - 1)
    Z_k <- Z_Delta[, idx_range, drop = FALSE]

    # 图中所需形式：log(sum(Z^2)) - log(n_k) - q0
    Z_sq_sum <- rowSums(Z_k^2)
    Y_k <- log(Z_sq_sum) - log(n_vec[k]) - q0

    Y_all[, k] <- Y_k
    col_idx <- col_idx + n_vec[k]
  }

  return(Y_all)  # 返回 m × K 的矩阵，每列是 Y_k^Δ(t_j)
}



# #' Estimate σ²_k(t) using Müller (2011) — cross-process smoothing, within-process recovery
# #'
# #' @param Y_all [m × N] matrix: all Y_{ij} = log(Z²) - q0
# #' @param Y_k [m × n_k] matrix: only kth process Y
# #' @param ts vector of time points
# #' @param L number of FPCA components
# #' @param q0 E[log(W²)] ~ -1.270
# #' @return vector of σ²_k(t)
# estimate_sigma2_Muller <- function(Y_all, Y_k, ts, L = 3, q0 = -1.270, h_min = 0.021) {
#   m <- nrow(Y_all) 
#   N <- ncol(Y_all)

#   ## Step 1: Estimate μ_V(t) from all Y_{ij}
#   Y_vec <- as.vector(Y_all)   # 所有 Y_{ij}
#   t_vec <- rep(ts, N)         # 相应的 t_j 重复 N 次
#   df_mu <- data.frame(t = t_vec, Y = Y_vec)
#   fit_mu <- npreg(Y ~ t, data = df_mu, regtype = "ll", bws = 0.05)
#   mu_V <- predict(fit_mu, newdata = data.frame(t = ts))

#   ## Step 2: Estimate G_V(s, t)
#   Y_centered <- Y_all - matrix(mu_V, m, N, byrow = FALSE)
#   raw_cov <- Y_centered %*% t(Y_centered) / N
#   coor1 <- rep(1:m, m) / m
#   coor2 <- rep(1:m, each = m) / m
#   idx_diag <- which(coor1 == coor2)
#   data_cov <- data.frame(Y = as.vector(raw_cov)[-idx_diag], X1 = coor1[-idx_diag], X2 = coor2[-idx_diag])
#   h <- max((N * m)^(-1/3), h_min)
#   fit_cov <- npreg(Y ~ X1 + X2, data = data_cov, regtype = "ll", bws = c(h, h) * 0.2)
#   G_V <- matrix(predict(fit_cov, newdata = data.frame(X1 = coor1, X2 = coor2)), m, m)

#   ## Step 3: FPCA
#   eig <- eigen(G_V, symmetric = TRUE)
#   phi_list <- lapply(1:L, function(i) Re(eig$vectors[, i]) * sqrt(m))

#   ## Step 4: Project kth process onto eigenbasis
#   Y_bar_k <- rowMeans(Y_k)
#   V_k <- mu_V
#   for (l in 1:L) {
#     phi_l <- phi_list[[l]]
#     xi_hat <- mean(Y_bar_k * phi_l)
#     V_k <- V_k + xi_hat * phi_l
#   }

#   ## Step 5: Return exp(V_k)
#   return(exp(V_k))
# }


# #' Perform FPCA on G_V(s, t) to obtain eigenfunctions and eigenvalues
# #'
# #' @param G_V matrix [m x m], estimated covariance function
# #' @param m number of grid points
# #' @param L_max maximum number of components to retain (default: 10)
# #'
# #' @return list with elements:
# #'   phi_list: list of length L, each element is vector φ_l(t) [length m]
# #'   lambda_vec: numeric vector [length L], eigenvalues
# #' @export
# fpca_from_G <- function(G_V, m, L_max = 10) {
#   eig <- eigen(G_V, symmetric = TRUE)

#   # 提取前 L 个实数主成分
#   phi_raw <- Re(eig$vectors[, 1:L_max]) * sqrt(m)  # 标准化
#   lambda_raw <- Re(eig$values[1:L_max]) / m        # 缩放回原比例

#   return(list(
#     phi_list = lapply(1:L_max, function(i) phi_raw[, i]),
#     lambda_vec = lambda_raw
#   ))
# }
