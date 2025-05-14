# R/recovery.R
# ========== FPCA 恢复 μ_i(t) ==========

recover_mu <- function(Z_Delta, m_mu, PCs, L, K){
  m <- nrow(Z_Delta)
  mu_hat <- matrix(0, m, K)

  for (i in 1:K) {
    mu_hat[, i] <- m_mu
    for (j in 1:L) {
      xi_hat <- mean((Z_Delta[, i] - m_mu) * PCs[, j])
      mu_hat[, i] <- mu_hat[, i] + xi_hat * PCs[, j]
    }
  }

  return(mu_hat)
}
