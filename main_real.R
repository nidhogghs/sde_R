# =================== main_real_fast_rows.R ====================
# 快速试跑：仅用前 50 行 (交易日) × 628 列 (股票)
# --------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(np)
  library(ggplot2)
  library(gridExtra)
})

# ---------- 参数 ----------
logret_path <- "output/logret_2015_2024.csv"
keep_rows   <- 50            # 只取前 50 行
h_fix       <- 0.05          # 固定带宽
workers     <- max(1, parallel::detectCores() - 1)
options(np.par = workers)    # np 包并行核数
# --------------------------

# ---------- 载入你的函数模块 ----------
# 载入配置与函数模块
source("config.R")               # 模型参数配置
source("R/simulation.R")         # 数据生成函数
source("R/estimation.R")         # μ 和 G_μ 的估计
source("R/recovery.R")           # 主成分重建函数
source("R/utils.R")              # 工具函数：如 cluster_mean, compute_qk_mc 等
source("R/mu_sigma_est.R")       # 封装的 estimate_mu_from_data / inference_sigma_from_data
source("R/traditional.R")

# ---------- 覆写核估计函数 (带 h_min 占位) ----------
estimate_m_mu <- function(Z_Delta, ts){
  mu_bar <- rowMeans(Z_Delta)
  fit <- npreg(Y ~ X,
               data    = data.frame(Y = mu_bar, X = ts),
               regtype = "ll",
               bws     = h_fix)
  predict(fit, newdata = data.frame(X = ts))
}

estimate_G_mu <- function(Z_Delta, m_mu_hat, m, K, h_min = 0.021){
  ts <- seq_len(m) / m
  coor1 <- rep(ts, m);  coor2 <- rep(ts, each = m)

  sample_cov <- (Z_Delta - m_mu_hat) %*% t(Z_Delta - m_mu_hat) / K
  idx_diag   <- coor1 == coor2
  data2 <- data.frame(
    Y  = as.vector(sample_cov)[!idx_diag],
    X1 = coor1[!idx_diag],
    X2 = coor2[!idx_diag]
  )

  fit <- npreg(Y ~ X1 + X2,
               data    = data2,
               regtype = "ll",
               bws     = c(h_fix, h_fix))
  matrix(predict(fit, newdata = data.frame(X1 = coor1, X2 = coor2)), m, m)
}

# ---------- 读取 Δlog 收益并裁切行 ----------
logret_dt <- fread(logret_path)
Z_full    <- as.matrix(logret_dt[, 2:ncol(logret_dt), with = FALSE])
Z_Delta   <- Z_full[1:keep_rows, ]        # ★ 仅前 50 行
m         <- nrow(Z_Delta)                # 50
K         <- ncol(Z_Delta)                # 628
delta     <- 1 / m
X_Delta   <- Z_Delta * sqrt(delta)

cat("✔ 使用维度: ", m, "×", K, " (rows × cols)\n")

# ---------- 组装 sim_data ----------
sim_data <- list(
  X_Delta = X_Delta,
  ts      = seq_len(m) / m,
  n_vec   = rep(1, K),
  delta   = delta
)

# ---------- 计时 & 估计 μ ----------
t0 <- Sys.time()
mu_res <- estimate_mu_from_data(sim_data)
cat("⏱ 总耗时:", round(difftime(Sys.time(), t0, units = "secs"), 1), "秒\n")

# ---------- 可视化 ----------
k_show <- 10
p_mu <- ggplot(
  data.frame(t = sim_data$ts, mu_hat = mu_res$mu_hat[, k_show]),
  aes(t, mu_hat)) +
  geom_line(color = "#1f78b4") +
  labs(title = paste0("μ̂_k(t) (stock #", k_show, ")"),
       x = "Normalized time", y = "μ̂_k(t)") +
  theme_minimal()

p_mmu <- ggplot(
  data.frame(t = sim_data$ts, m_mu_hat = mu_res$m_mu_hat),
  aes(t, m_mu_hat)) +
  geom_line(color = "#33a02c") +
  labs(title = paste0("m_μ(t) over ", K, " stocks"),
       x = "Normalized time", y = "m_μ(t)") +
  theme_minimal()

grid.arrange(p_mu, p_mmu, ncol = 2)
