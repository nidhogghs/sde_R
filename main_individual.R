# ======= main_compare_n10_n100.R =======

# ---- 载入配置与函数模块 ----
source("config.R")               # 模型参数配置
source("R/simulation.R")         # 数据生成函数
source("R/estimation.R")         # μ 和 G_μ 的估计
source("R/recovery.R")           # 主成分重建函数
source("R/utils.R")              # 工具函数：如 cluster_mean, compute_qk_mc 等
source("R/mu_sigma_est.R")       # 封装的 estimate_mu_from_data / inference_sigma_from_data
source("R/traditional.R")        # 如果内部依赖了绘图函数，保留；本脚本不会调用传统估计

# ---- 绘图相关库 ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(reshape2)
})

# ---- 输出目录 ----
out_dir <- "output/compare_n10_n100"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 一个小工具：生成“真实 vs 估计”的 μ_k(t) 对比图 ----
plot_mu_compare <- function(ts, mu_true, mu_hat, k, title_suffix = "") {
  df <- data.frame(
    t = ts,
    True = mu_true[, k],
    Estimated = mu_hat[, k]
  )
  dfl <- reshape2::melt(df, id.vars = "t", variable.name = "Type", value.name = "Value")
  ggplot(dfl, aes(x = t, y = Value, color = Type)) +
    geom_line(linewidth = 1) +
    labs(
      title = sprintf("mu_k(t): True vs Estimated (k = %d)%s", k, title_suffix),
      x = "t", y = "mu_k(t)"
    ) +
    theme_minimal()
}

# ---- 一个小工具：绘制 μ_k(t) 的置信区间（使用 compute_mu_ci 的输出）----
plot_mu_ci_band <- function(ts, mu_true, mu_hat, ci_lower, ci_upper, k, title_suffix = "") {
  df <- data.frame(
    t = ts,
    mu_true = mu_true[, k],
    mu_hat  = mu_hat[, k],
    lower   = ci_lower[, k],
    upper   = ci_upper[, k]
  )
  ggplot(df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(aes(y = mu_hat), linewidth = 1) +
    geom_line(aes(y = mu_true), linewidth = 0.9, linetype = "dashed") +
    labs(
      title = sprintf("mu_k(t) with 95%% CI (k = %d)%s", k, title_suffix),
      x = "t", y = "mu_k(t)"
    ) +
    theme_minimal()
}

# ---- 主流程：给定 n_ave，生成数据并绘制 k=1,2 的对比与置信区间 ----
run_once <- function(n_ave, K = 500, m = 50, k_true = 4, l_true = 4, scale_sigma = 5, alpha = 0.05) {
  message(sprintf(">>> Running with n_ave = %d ...", n_ave))
  
  # 1) 模拟数据
  sim_data <- prepare_simulation_data(
    K = K,
    n_ave = n_ave,   # 每个过程的子样本数平均水平
    m = m,
    k = k_true,
    l = l_true,
    scale_sigma = scale_sigma
  )
  
  # 2) 估计 μ
  mu_res <- estimate_mu_from_data(sim_data)
  
  # 3) 推断 σ²（compute_mu_ci 需要）
  sigma_res <- inference_sigma_from_data(sim_data)
  
  # 4) 计算 μ 的置信区间
  ci_mu <- compute_mu_ci(
    mu_hat   = mu_res$mu_hat,
    PCs      = mu_res$PCs_hat,
    sigma2_hat = sigma_res$sigma2_hat,
    ts       = sim_data$ts,
    n_vec    = sim_data$n_vec,
    alpha    = alpha,
    L        = mu_res$L
  )
  
  # 5) 仅绘制 k = 1, 2 两条轨迹
  ks <- c(2, 5)
  grobs <- list()
  idx <- 1
  for (kk in ks) {
    # 5.1 真 vs 估计
    p_compare <- plot_mu_compare(
      ts = sim_data$ts,
      mu_true = sim_data$mu_true,
      mu_hat  = mu_res$mu_hat,
      k = kk,
      title_suffix = sprintf(" | n_ave=%d", n_ave)
    )
    # 5.2 置信区间图
    p_ci <- plot_mu_ci_band(
      ts = sim_data$ts,
      mu_true = sim_data$mu_true,
      mu_hat  = mu_res$mu_hat,
      ci_lower = ci_mu$lower,
      ci_upper = ci_mu$upper,
      k = kk,
      title_suffix = sprintf(" | n_ave=%d", n_ave)
    )
    grobs[[idx]]     <- p_compare; idx <- idx + 1
    grobs[[idx]]     <- p_ci;      idx <- idx + 1
  }
  
  # 6) 排版与导出
  grid_title <- sprintf("Comparison for n_ave = %d (k = 1, 2)", n_ave)
  g <- do.call(grid.arrange, c(grobs, ncol = 2, top = grid_title))
  
  out_file <- file.path(out_dir, sprintf("mu_compare_ci_n%d_k1k2.png", n_ave))
  ggsave(out_file, g, width = 14, height = 8, dpi = 150)
  message(sprintf("Saved: %s", out_file))
}

# ---- 实际运行：n_ave = 10 和 100 ----
run_once(n_ave = 30)
run_once(n_ave = 300)

# ===== End of main_compare_n10_n100.R =====
