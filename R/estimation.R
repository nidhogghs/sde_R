# R/estimation.R
# ========== 核估计与FPCA ==========

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

#' Estimate σ²_k(t) via Müller (2011) within-process method
#'
#' @param Z_k matrix [m × n_k] of Δlog X for class k
#' @param ts vector of time points
#' @param L number of FPCA components
#' @param q0 constant E[log(W²)] ~ -1.270
#' @return vector of σ²_k(t) over ts
#' @export
estimate_sigma2_k <- function(Z_k, ts, L = 3, q0 = -1.270, h_min = 0.021) {
  m <- nrow(Z_k)
  n_k <- ncol(Z_k)
  
  # Step 1: log(Z²) - q0
  Y_k <- log(Z_k^2) - q0
  
  # Step 2: μ_V(t)
  Y_vec <- as.vector(Y_k)
  t_vec <- rep(ts, n_k)
  df_mu <- data.frame(t = t_vec, y = Y_vec)
  fit_mu <- npreg(y ~ t, data = df_mu, regtype = "ll", bwmethod = "cv.aic")
  mu_V <- predict(fit_mu, newdata = data.frame(t = ts))
  
  # # Step 3: G_V(s,t)
  # Y_centered <- Y_k - matrix(mu_V, m, n_k, byrow = FALSE)
  # raw_cov <- Y_centered %*% t(Y_centered) / n_k
  # coor1 <- rep(1:m, m) / m
  # coor2 <- rep(1:m, each = m) / m
  # idx_diag <- which(coor1 == coor2)
  # df_cov <- data.frame(
  #   Y = as.vector(raw_cov)[-idx_diag],
  #   X1 = coor1[-idx_diag],
  #   X2 = coor2[-idx_diag]
  # )
  # h <- max((n_k * m)^(-1/3), h_min)
  # fit_cov <- npreg(Y ~ X1 + X2, data = df_cov, regtype = "ll", bws = c(h, h) * 0.2)
  # G_V <- matrix(predict(fit_cov, newdata = data.frame(X1 = coor1, X2 = coor2)), m, m)
  
  # # Step 4: FPCA
  # eig <- eigen(G_V, symmetric = TRUE)
  # phi_list <- lapply(1:L, function(i) Re(eig$vectors[, i]) * sqrt(m))
  
  # Step 5: Reconstruct V_k(t)
  V_k <- mu_V
  # for (l in 1:L) {
  #   xi_hat <- mean(rowMeans(Y_k) * phi_list[[l]])
  #   V_k <- V_k + xi_hat * phi_list[[l]]
  # }
  
  # Step 6: Return σ²_k(t)
  return(exp(V_k))
}


#' Perform FPCA on G_V(s, t) to obtain eigenfunctions and eigenvalues
#'
#' @param G_V matrix [m x m], estimated covariance function
#' @param m number of grid points
#' @param L_max maximum number of components to retain (default: 10)
#'
#' @return list with elements:
#'   phi_list: list of length L, each element is vector φ_l(t) [length m]
#'   lambda_vec: numeric vector [length L], eigenvalues
#' @export
fpca_from_G <- function(G_V, m, L_max = 10) {
  eig <- eigen(G_V, symmetric = TRUE)

  # 提取前 L 个实数主成分
  phi_raw <- Re(eig$vectors[, 1:L_max]) * sqrt(m)  # 标准化
  lambda_raw <- Re(eig$values[1:L_max]) / m        # 缩放回原比例

  return(list(
    phi_list = lapply(1:L_max, function(i) phi_raw[, i]),
    lambda_vec = lambda_raw
  ))
}

#' Reconstruct V_k(t) from Y_k using FPCA scores and basis
#'
#' @param Y_k matrix [m × n_k] of log-squared Z_{ij} - q0 for class k
#' @param mu_V vector of length m, estimated μ_V(t)
#' @param phi_list list of φ_l(t), each a vector of length m
#' @param L number of components to use
#'
#' @return vector of length m, reconstructed V_k(t)
#' @export
reconstruct_V_k <- function(Y_k, mu_V, phi_list, L) {
  m <- nrow(Y_k)
  Y_bar <- rowMeans(Y_k)

  V_k_hat <- mu_V
  for (l in 1:L) {
    phi_l <- phi_list[[l]]
    xi_hat <- mean(Y_bar * phi_l)  # inner product approximation
    V_k_hat <- V_k_hat + xi_hat * phi_l
  }

  return(V_k_hat)
}
