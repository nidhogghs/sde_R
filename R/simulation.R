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

generate_K_trajectory <- function(n, a, delta, alpha, k, N){
  sapply(1:n, function(i) generate_one_trajectory(a, delta, alpha, k, N))
}

generate_n_X <- function(n, sigma, mu, X0){
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
