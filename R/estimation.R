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


#' Estimate σ²_k(t) using Müller et al. (2011) log-squared method
#'
#' This method estimates σ²_k(t) by first computing log(Z^2) - q0,
#' then smoothing it via local linear regression, and exponentiating the result.
#'
#' @param Z_k matrix of shape [m x n_k], each column is Z_Δ of one sample in group k
#' @param ts vector of time points t_j (length m)
#' @param q0 optional adjustment constant, default -1.270 (E[log(W²)] for N(0,1))
#'
#' @return numeric vector of length m, the estimated σ²_k(t)
#' @export
estimate_sigma2_muller <- function(Z_k, ts, q0 = -1.270) {
  m <- nrow(Z_k)
  n_k <- ncol(Z_k)

  # Step 1: construct log(Z^2) - q0
  Y_mat <- log(Z_k^2) - q0  # matrix [m x n_k]

  # Step 2: take column-wise average over n_k samples
  Y_bar <- rowMeans(Y_mat)

  # Step 3: smooth Y_bar over t_j using local linear regression
  data <- data.frame(Y = Y_bar, t = ts)
  fit <- npreg(Y ~ t, data = data, regtype = "ll", bwmethod = "cv.aic")

  # Step 4: predict V(t) = smoothed log σ²(t), then exponentiate
  V_hat <- predict(fit, newdata = data.frame(t = ts))
  sigma2_hat <- exp(V_hat)

  return(sigma2_hat)
}

