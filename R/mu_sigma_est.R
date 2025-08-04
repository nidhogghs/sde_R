# R/mu_sigma_est.R

estimate_mu_from_data <- function(sim_data, L = NULL) {

  ts <- sim_data$ts
  m <- length(ts)
  K <- length(sim_data$n_vec)

  if (max(sim_data$n_vec) == 1) {
    Z_Delta <- sim_data$X_Delta / sqrt(sim_data$delta)
  } else {
    Z_Delta <- cluster_mean(sim_data$X_Delta, K, sim_data$n_vec) / sqrt(sim_data$delta)
  }

  m_mu_hat <- estimate_m_mu(Z_Delta, ts)
  G_mu_hat <- estimate_G_mu(Z_Delta, m_mu_hat, m, K, h_min = h_min)
  PCs_hat <- Re(eigen(G_mu_hat)$vectors) * sqrt(m)
  lams_hat <- Re(eigen(G_mu_hat)$values) / m

  if (is.null(L)) {
    L <- select_L_by_AIC(Z_Delta, m_mu_hat, PCs_hat, max_L = 10)
    cat(sprintf("Selected L (mu) = %d\n", L))
  }

  mu_hat <- recover_mu(Z_Delta, m_mu_hat, PCs_hat, L, K)

  return(list(
    mu_hat = mu_hat,
    m_mu_hat = m_mu_hat,
    G_mu_hat = G_mu_hat,
    PCs_hat = PCs_hat,
    lams_hat = lams_hat,
    L = L
  ))
}

inference_sigma_from_data <- function(sim_data, L = NULL) {
  ts <- sim_data$ts
  m <- length(ts)
  K <- length(sim_data$n_vec)

  q_vec <- compute_qk(sim_data$n_vec)
  Y_Delta <- construct_Y(sim_data$X_Delta, sim_data$n_vec, q_vec)

  m_V_hat <- estimate_m_mu(Y_Delta, ts)
  G_V_hat <- estimate_G_mu(Y_Delta, m_V_hat, m, K, h_min = h_min)
  PCs_hat <- Re(eigen(G_V_hat)$vectors) * sqrt(m)
  lams_hat <- Re(eigen(G_V_hat)$values) / m

  if (is.null(L)) {
    L <- select_L_by_AIC(Y_Delta, m_V_hat, PCs_hat, max_L = 10)
    cat(sprintf("Selected L (sigma²) = %d\n", L))
  }

  V_hat <- recover_mu(Y_Delta, m_V_hat, PCs_hat, L, K)
  sigma2_hat <- exp(V_hat)

  return(list(
    V_hat = V_hat,
    sigma2_hat = sigma2_hat,
    m_V_hat = m_V_hat,
    G_V_hat = G_V_hat,
    PCs_hat = PCs_hat,
    lams_hat = lams_hat,
    L = L
  ))
}

compute_mu_ci <- function(mu_hat, PCs, sigma2_hat, ts, n_vec, alpha = 0.05, L) {
  m <- length(ts)
  K <- ncol(mu_hat)
  z_alpha <- qnorm(1 - alpha / 2)

  # 初始化上下界矩阵
  ci_lower <- matrix(NA, m, K)
  ci_upper <- matrix(NA, m, K)

  for (k in 1:K) {
    sigma2_k <- sigma2_hat[, k]
    n_k <- n_vec[k]

    # 构造 V_μ(t)
    V_mu_t <- numeric(m)
    for (t in 1:m) {
      sum_term <- 0
      for (j in 1:m) {
        phi_sum <- sum(PCs[j, 1:L] * PCs[t, 1:L]) 
        sum_term <- sum_term + (phi_sum^2) * sigma2_k[j]
      }
      V_mu_t[t] <- sum_term / m
    }

    se <- sqrt(V_mu_t / n_k)
    ci_lower[, k] <- mu_hat[, k] - z_alpha * se
    ci_upper[, k] <- mu_hat[, k] + z_alpha * se
  }

  return(list(lower = ci_lower, upper = ci_upper))
}

