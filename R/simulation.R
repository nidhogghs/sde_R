# R/simulation.R
# ========== 模拟相关函数 ==========

generate_one_trajectory <- function(a, delta, alpha, k, N){
  kappa <- function(x) 100 * x^3 * (1 - x)^3
  seq1 <- seq(0, 1, 1 / N)
  part1 <- sapply(seq1, kappa)

  phi1 <- function(x, j) sqrt(2) * sin(j * pi * x)
  phi2 <- function(x, j) sqrt(2) * cos(j * pi * x)

  part2 <- rep(0, length(seq1))
  for (j in 1:k) {
    xi <- rnorm(1, 0, sqrt(2^(1 - j)))
    j1 <- (j + 1) %/% 2
    if (j %% 2 == 1) part2 <- part2 + xi * sapply(seq1, function(x) phi1(x, j1))
    if (j %% 2 == 0) part2 <- part2 + xi * sapply(seq1, function(x) phi2(x, j1))
  }

  return(a + part1 + alpha * part2)
}



#返回值N+1*K的矩阵
generate_K_trajectory <- function(K, a, delta, alpha, k, N){
  sapply(1:K, function(i) generate_one_trajectory(a, delta, alpha, k, N))
}

#生成轨迹
generate_n_X <- function(n, sigma, mu, X0){#矩阵，维度为 (m+1) × ∑ n_k
  K <- ncol(sigma)
  m <- nrow(sigma) - 1
  Delta <- T / m
  X <- rnorm(sum(n), X0, 1)
  X <- 1 + X - min(X)

  d_log_X <- matrix(NA, nrow = m, ncol = sum(n))
  col_idx <- 1
  for (k in 1:K) {
    for (j in 1:n[k]) {
      d_log_X[, col_idx] <- mu[2:(m+1), k] * Delta + 
        sigma[2:(m+1), k] * rnorm(m, 0, sqrt(Delta))
      col_idx <- col_idx + 1
    }
  }

  d_log_X <- rbind(log(X), d_log_X)
  log_X <- apply(d_log_X, 2, cumsum)
  return(exp(log_X))
}

variation <- function(X, m){
  Delta <- T / m
  y <- apply(X, 2, diff) / sqrt(Delta)
  return(y)
}

prepare_simulation_data <- function(K, n_ave = NULL, n_vec = NULL, m, k = 2, l = 2,
                                    seed_sigma = 311, seed_mu = 311, seed_X = 69, scale_sigma = 5) {
  if (is.null(n_vec)) {
    if (is.null(n_ave)) stop("You must provide either n_vec or n_ave.")
    n_vec <- rep(n_ave, K)
  } else {
    if (length(n_vec) != K) stop("Length of n_vec must be equal to K.")
  }

  delta <- T / m
  ts <- seq(delta, T, length.out = m)

  set.seed(seed_sigma)
  V_true <- generate_K_trajectory(K, a, delta, alpha, k, m)
  V_scaled <- V_true  # 保留未缩放版本以供恢复 σ²

  sigma2 <- exp(V_scaled)
  sigma <- sqrt(sigma2) / scale_sigma  # 控制扩散强度
  sigma2 <- sigma^2

  set.seed(seed_mu)
  mu <- generate_K_trajectory(K, b, delta, beta, l, m)

  set.seed(seed_X)
  X <- generate_n_X(n_vec, sigma, mu, X0)

  X_Delta <- variation(log(X), m)

  list(
    ts = ts,
    delta = delta,
    n_vec = n_vec,
    V_true = V_scaled,
    sigma2_true = sigma2[-1, ],
    mu_true = mu[-1, ],
    X_Delta = X_Delta
  )
}

prepare_simulation_data2 <- function(K, n_ave = NULL, n_vec = NULL, m, k = 2, l = 2,
                                    seed_sigma = 311, seed_mu = 311, seed_X = 69, scale_sigma = 5) {
  if (is.null(n_vec)) {
    if (is.null(n_ave)) stop("You must provide either n_vec or n_ave.")
    n_vec <- rep(n_ave, K)
  } else {
    if (length(n_vec) != K) stop("Length of n_vec must be equal to K.")
  }

  delta <- T / m
  ts <- seq(delta, T, length.out = m)

  set.seed(seed_sigma)
  sigma <- generate_K_trajectory(K, a, delta, alpha, k, m)
  sigma <- sigma / max(sigma)  # 归一化防止数值不稳定
  sigma2 <- sigma^2

  V_true <- log(sigma2)
  V_scaled <- V_true


  set.seed(seed_mu)
  mu <- generate_K_trajectory(K, b, delta, beta, l, m)

  set.seed(seed_X)
  X <- generate_n_X(n_vec, sigma, mu, X0)

  X_Delta <- variation(log(X), m)

  list(
    ts = ts,
    delta = delta,
    n_vec = n_vec,
    V_true = V_scaled,
    sigma2_true = sigma2[-1, ],
    mu_true = mu[-1, ],
    X_Delta = X_Delta
  )
}