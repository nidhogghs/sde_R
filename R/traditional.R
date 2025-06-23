# R/traditional.R

estimate_mu_traditional_from_data <- function(sim_data, bw_method = "cv.aic") {
  X_Delta <- sim_data$X_Delta
  ts <- sim_data$ts
  delta <- sim_data$delta
  n_vec <- sim_data$n_vec
  m <- length(ts)
  K <- length(n_vec)

  # Step 1: 取均值
  if (max(n_vec) == 1) {
    Z_Delta <- X_Delta
  } else {
    Z_Delta <- cluster_mean(X_Delta, K, n_vec) / sqrt(delta)
  }

  # Step 2: 对每条聚合轨迹核平滑
  mu_hat_mat <- matrix(NA, m, K)
  for (k in 1:K) {
    y_k <- Z_Delta[, k]
    df <- data.frame(t = ts, y = y_k)
    fit <- npreg(y ~ t, data = df, regtype = "ll", bwmethod = bw_method)
    mu_hat_mat[, k] <- predict(fit, newdata = data.frame(t = ts))
  }

  return(mu_hat_mat)
}
evaluate_traditional_mu_estimation <- function(sim_data, mu_hat_mat, show_plot_k = c(1)) {
  ts <- sim_data$ts
  mu_true <- sim_data$mu_true
  m <- nrow(mu_hat_mat)
  K <- ncol(mu_hat_mat)

  # 误差评估
  mse_vec <- colMeans((mu_hat_mat - mu_true)^2)
  rmse_vec <- sqrt(mse_vec)
  avg_rmse <- mean(rmse_vec)
  max_rmse <- max(rmse_vec)

  cat("=== Traditional Method RMSE Results ===\n")
  cat(sprintf("Average RMSE across K: %.6f\n", avg_rmse))
  cat(sprintf("Maximum RMSE among K: %.6f\n", max_rmse))

  # 可视化
  par(mfrow = c(length(show_plot_k), 1), mar = c(4, 4, 2, 1))
  for (k in show_plot_k) {
    plot(ts, mu_true[, k], type = "l", lwd = 2, col = "black",
         ylab = bquote(mu[.(k)](t)), xlab = "t",
         main = paste0("Process k = ", k, ": True vs Traditional"))
    lines(ts, mu_hat_mat[, k], col = "red", lwd = 2, lty = 2)
    legend("topright", legend = c("True", "Traditional"), col = c("black", "red"), lty = c(1, 2), lwd = 2)
  }

  invisible(list(
    rmse = rmse_vec,
    avg_rmse = avg_rmse,
    mu_hat = mu_hat_mat
  ))
}
