#先生成指数，再产生sigma
prepare_simulation_data <- function(K, n_ave, m, k = 2, l = 2,
                                    seed_sigma =789, seed_mu = 456, seed_X = 69) {
  delta <- T / m
  ts <- seq(delta, T, length.out = m)
  n_vec <- rep(n_ave, K)

  set.seed(seed_sigma)
  V_true <- generate_K_trajectory(K, a, delta, alpha, k, m)
  sigma2 <- exp(V_true)
  sigma <- sqrt(sigma2)

  set.seed(seed_mu)
  mu <- generate_K_trajectory(K, b, delta, beta, l, m)

  set.seed(seed_X)
  X <- generate_n_X(n_vec, sigma, mu, X0)

  X_Delta <- variation(log(X), m)

  list(
    ts = ts,
    delta = delta,
    n_vec = n_vec,
    V_true = V_true,
    sigma2_true = sigma2[-1, ],
    mu_true = mu[-1, ],
    X_Delta = X_Delta
  )
}

#直接生成sigma，此部分仅用于测试对μ的估计
prepare_simulation_data2 <- function(K, n_ave, m, k = 2, l = 2, 
                                    seed_sigma = 789, seed_mu = 456, seed_X = 69) {
  delta <- T / m
  ts <- seq(delta, T, length.out = m)
  n_vec <- rep(n_ave, K)

  set.seed(seed_sigma)
  V_true <- generate_K_trajectory(K, a, delta, alpha, k, m)#此处仅仅保持接口统一，未使用
  sigma2 <- exp(V_true)
  sigma <- generate_K_trajectory(K, a, delta, alpha, k, m)
  sigma <- sigma / max(sigma)  # 归一化防止数值不稳定

  set.seed(seed_mu)
  mu <- generate_K_trajectory(K, b, delta, beta, l, m)

  set.seed(seed_X)
  X <- generate_n_X(n_vec, sigma, mu, X0)

  X_Delta <- variation(log(X), m)

  list(
    ts = ts,
    delta = delta,
    n_vec = n_vec,
    V_true = V_true,
    sigma2_true = sigma2[-1, ],
    mu_true = mu[-1, ],
    X_Delta = X_Delta
  )
}

#先生成指数，再对μ估计，结果不好
test_mu_estimation <- function(K = 300, n_ave = 500, m = 50, L = 2) {
  sim_data <- prepare_simulation_data(K, n_ave, m, l = 2)
  ts <- sim_data$ts
  X_Delta <- sim_data$X_Delta

  if (max(sim_data$n_vec) == 1) {
    Z_Delta <- X_Delta
  } else {
    Z_Delta <- cluster_mean(X_Delta, K, sim_data$n_vec) / sqrt(sim_data$delta)
  }

  # Step 1: estimate m_mu
  m_mu_hat <- estimate_m_mu(Z_Delta, ts)

  # Step 2: estimate G_mu
  G_mu_hat <- estimate_G_mu(Z_Delta, m_mu_hat, m, K, h_min = h_min)
  PCs_hat <- Re(eigen(G_mu_hat)$vectors) * sqrt(m)
  lams_hat <- Re(eigen(G_mu_hat)$values) / m

  # Step 3: select L
  if (is.null(L)) {
    L <- select_L_by_AIC(Z_Delta, m_mu_hat, PCs_hat, max_L = 10)
    cat(sprintf("Selected L = %d\n", L))
  }

  # Step 4: recover mu_k(t)
  mu_hat <- recover_mu(Z_Delta, m_mu_hat, PCs_hat, L, K)

  # Step 5: compute errors
  mu_true <- sim_data$mu_true
  m_mu_true <- rowMeans(mu_true)
  G_mu_true <- (mu_true - m_mu_true) %*% t(mu_true - m_mu_true) / K

  PCs_true <- Re(eigen(G_mu_true)$vectors) * sqrt(m)
  lams_true <- Re(eigen(G_mu_true)$values) / m

  err_m <- max(abs(m_mu_hat - m_mu_true))
  err_G <- max(abs(G_mu_hat - G_mu_true))
  err_lams <- max(abs(lams_hat[1:length(lams_true)] - lams_true))

  err_PCs <- max(sapply(1:L, function(i) {
    true <- PCs_true[, i] * sign(PCs_true[1, i])
    est  <- PCs_hat[, i] * sign(PCs_hat[1, i])
    max(abs(true - est))
  }))

  err_mu <- max(abs(mu_hat[, 1] - mu_true[, 1])) / (max(mu_true[, 1]) - min(mu_true[, 1]))

  cat(sprintf("err_m_mu         = %0.6f\n", err_m))
  cat(sprintf("err_G_mu         = %0.6f\n", err_G))
  cat(sprintf("err_lams         = %0.6f\n", err_lams))
  cat(sprintf("err_PCs          = %0.6f\n", err_PCs))
  cat(sprintf("err_mu_k         = %0.6f\n", err_mu))

  # Step 6: plot comparison for 1st process
  plot_compare_single(ts, mu_true, mu_hat, k = 1, label = "mu_k(t)")

  # Return for comparison if needed
  invisible(list(
    mu_true = mu_true,
    mu_hat = mu_hat,
    err_m   = err_m,
    err_G   = err_G,
    err_lams= err_lams,
    err_PCs = err_PCs,
    err_mu  = err_mu
  ))
}


#直接生成sigma，估计μ，结果较好
test_mu_estimation2 <- function(K = 300, n_ave = 500, m = 50, L = 2) {
  sim_data <- prepare_simulation_data2(K, n_ave, m, l = 2)
  ts <- sim_data$ts
  X_Delta <- sim_data$X_Delta

  if (max(sim_data$n_vec) == 1) {
    Z_Delta <- X_Delta
  } else {
    Z_Delta <- cluster_mean(X_Delta, K, sim_data$n_vec) / sqrt(sim_data$delta)
  }

  # Step 1: estimate m_mu
  m_mu_hat <- estimate_m_mu(Z_Delta, ts)

  # Step 2: estimate G_mu
  G_mu_hat <- estimate_G_mu(Z_Delta, m_mu_hat, m, K, h_min = h_min)
  PCs_hat <- Re(eigen(G_mu_hat)$vectors) * sqrt(m)
  lams_hat <- Re(eigen(G_mu_hat)$values) / m

  # Step 3: select L
  if (is.null(L)) {
    L <- select_L_by_AIC(Z_Delta, m_mu_hat, PCs_hat, max_L = 10)
    cat(sprintf("Selected L = %d\n", L))
  }

  # Step 4: recover mu_k(t)
  mu_hat <- recover_mu(Z_Delta, m_mu_hat, PCs_hat, L, K)

  # Step 5: compute errors
  mu_true <- sim_data$mu_true
  m_mu_true <- rowMeans(mu_true)
  G_mu_true <- (mu_true - m_mu_true) %*% t(mu_true - m_mu_true) / K

  PCs_true <- Re(eigen(G_mu_true)$vectors) * sqrt(m)
  lams_true <- Re(eigen(G_mu_true)$values) / m

  err_m <- max(abs(m_mu_hat - m_mu_true))
  err_G <- max(abs(G_mu_hat - G_mu_true))
  err_lams <- max(abs(lams_hat[1:length(lams_true)] - lams_true))

  err_PCs <- max(sapply(1:L, function(i) {
    true <- PCs_true[, i] * sign(PCs_true[1, i])
    est  <- PCs_hat[, i] * sign(PCs_hat[1, i])
    max(abs(true - est))
  }))

  err_mu <- max(abs(mu_hat[, 1] - mu_true[, 1])) / (max(mu_true[, 1]) - min(mu_true[, 1]))

  cat(sprintf("err_m_mu         = %0.6f\n", err_m))
  cat(sprintf("err_G_mu         = %0.6f\n", err_G))
  cat(sprintf("err_lams         = %0.6f\n", err_lams))
  cat(sprintf("err_PCs          = %0.6f\n", err_PCs))
  cat(sprintf("err_mu_k         = %0.6f\n", err_mu))

  # Step 6: plot comparison for 1st process
  plot_compare_single(ts, mu_true, mu_hat, k = 1, label = "mu_k(t)")

  # Return for comparison if needed
  invisible(list(
    mu_true = mu_true,
    mu_hat = mu_hat,
    err_m   = err_m,
    err_G   = err_G,
    err_lams= err_lams,
    err_PCs = err_PCs,
    err_mu  = err_mu
  ))
}

#生成指数再估计sigma，结果较好
test_sigma_inference <- function(K = 100, n_ave = 300, m = 50, L = NULL) {
  sim_data <- prepare_simulation_data(K, n_ave, m,l=2)
  q_vec <- compute_qk_mc(sim_data$n_vec)
  Y_Delta <- construct_Y(sim_data$X_Delta, sim_data$n_vec, q_vec)

  m_V_hat <- estimate_m_mu(Y_Delta, sim_data$ts)
  G_V_hat <- estimate_G_mu(Y_Delta, m_V_hat, m, K, h_min = h_min)
  PCs_hat <- Re(eigen(G_V_hat)$vectors) * sqrt(m)
  lams_hat <- Re(eigen(G_V_hat)$values) / m

  if (is.null(L)) {
    L <- select_L_by_AIC(Y_Delta, m_V_hat, PCs_hat, max_L = 10)
    cat(sprintf("Selected L = %d\n", L))
  }

  V_hat <- recover_mu(Y_Delta, m_V_hat, PCs_hat, L, K)
  sigma2_hat <- exp(V_hat)

  plot_compare_single(sim_data$ts, sim_data$V_true[-1, ], V_hat, k = 1, label = "V_k(t)")
  plot_compare_single(sim_data$ts, sim_data$sigma2_true, sigma2_hat, k = 1, label = "sigma²_k(t)")
}


