# R/inference.R
one_run_mu_inference <- function(K, n_ave, m, sd, L = NULL) {
  delta <- T / m
  ts <- seq(delta, T, length.out = m)
  n_vec <- rep(n_ave, K)

  # Step 1: simulate true log sigma²
  set.seed(456)
  V_true <- generate_K_trajectory(K, a, delta, alpha, L, m)  # (m+1) × K
  sigma2 <- exp(V_true)
  sigma  <- sqrt(sigma2)

  m_V_true <- rowMeans(V_true)[-1]  # m × 1
  G_V_true <- (V_true[-1, ] - m_V_true) %*% t(V_true[-1, ] - m_V_true) / K
  PCs_true <- Re(eigen(G_V_true)$vectors) * sqrt(m)
  lams_true <- Re(eigen(G_V_true)$values) / m
  sigma2_true <- sigma2[-1, ]

  # Step 2: simulate mu(t) and X(t)
  set.seed(789)
  mu <- generate_K_trajectory(K, b, delta, beta, L, m)

  set.seed(sd)
  X <- generate_n_X(n_vec, sigma, mu, X0)

  # Step 3: construct log(Z²) responses
  X_Delta <- variation(log(X), m)
  q0 <- compute_qk_mc(n_ave)
  Y_Delta <- construct_Y(X_Delta, n_vec, q0)

  # Step 4: estimate mean + covariance
  m_V_hat <- estimate_m_mu(Y_Delta, ts)
  G_V_hat <- estimate_G_mu(Y_Delta, m_V_hat, m, K, h_min = h_min)

  # Step 5: FPCA
  PCs_hat <- Re(eigen(G_V_hat)$vectors) * sqrt(m)
  lams_hat <- Re(eigen(G_V_hat)$values) / m

  # Step 6: select L
  if (is.null(L)) {
    L <- select_L_by_AIC(Y_Delta, m_V_hat, PCs_hat, max_L = 10)
    if (is.na(L)) stop("AIC failed.")
  }

  # Step 7: estimate V_k(t)
  V_hat <- recover_mu(Y_Delta, m_V_hat, PCs_hat, L, K)
  sigma2_hat <- exp(V_hat)

  # Step 8: error evaluation
  err_m <- max(abs(m_V_hat - m_V_true))
  err_G <- max(abs(G_V_hat - G_V_true))
  err_lams <- max(abs(lams_hat[1:length(lams_true)] - lams_true))
  err_PCs <- max(sapply(1:l, function(i) {
    true <- PCs_true[, i] * sign(PCs_true[1, i])
    est <- PCs_hat[, i] * sign(PCs_hat[1, i])
    max(abs(true - est))
  }))
  err_sigma2 <- max(abs(sigma2_hat[, 1] - sigma2_true[, 1])) / (max(sigma2_true[, 1]) - min(sigma2_true[, 1]))

  return(c(err_m, err_G, err_lams, err_PCs, err_sigma2))
}
