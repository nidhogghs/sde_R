# ===================== main_real_strict.R =====================
# 读取 output/logret_2015_2024.csv  →  构造 X_Delta  → 估计 μ
# -------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

# ---- 1. 载入你的函数模块 ----
source("config.R")
source("R/estimation.R")
source("R/mu_sigma_est.R")   # 如需要 σ² 推断

# ---- 2. 读取 Δlog 收益并还原 X_Delta ----
logret_path <- "output/logret_2015_2024.csv"
logret_dt   <- fread(logret_path)

# 去掉日期列，确保仍是 data.table，再转 matrix
Z_Delta <- as.matrix(logret_dt[, 2:ncol(logret_dt), with = FALSE])
m <- nrow(Z_Delta)             # 247
K <- ncol(Z_Delta)             # 628
delta <- 1 / m                 # Δt   (你的代码里要 sqrt(delta))

# ★ 恢复 X_Delta：Z_Delta = X_Delta / √Δt  ⇒   X_Delta = Z_Delta * √Δt
X_Delta <- Z_Delta * sqrt(delta)

cat("✔ X_Delta dim:", m, "×", K, "\n")

# ---- 3. 组装 sim_data ----
sim_data <- list(
  X_Delta = X_Delta,   # ★ 关键字段
  ts      = seq_len(m) / m,
  n_vec   = rep(1, K), # 每股票 n_k = 1
  delta   = delta
)

# ---- 4. 估计 μ_k(t) ----
cat("🚀 估计 μ_k(t)…\n")
mu_res <- estimate_mu_from_data(sim_data)   # 会自动选 L 并打印
cat("✔ 完成 μ 估计，L_hat =", mu_res$L, "\n")

# ---- 5. (可选) σ² 推断 ----
# sigma_res <- inference_sigma_from_data(sim_data)

# ---- 6. 绘图：第 1 列 μ̂₁(t) & 整体 m_μ(t) ----
k_show <- 1
df_mu  <- data.frame(t = sim_data$ts,
                     mu_hat = mu_res$mu_hat[, k_show])
p_mu <- ggplot(df_mu, aes(t, mu_hat)) +
  geom_line(color = "#1f78b4") +
  labs(title = paste0("μ̂_k(t) · stock #", k_show),
       x = "Normalized time", y = "μ̂_k(t)") +
  theme_minimal()

df_mmu <- data.frame(t = sim_data$ts,
                     m_mu_hat = mu_res$m_mu_hat)
p_mmu <- ggplot(df_mmu, aes(t, m_mu_hat)) +
  geom_line(color = "#33a02c") +
  labs(title = paste0("m_μ(t) over ", K, " stocks"),
       x = "Normalized time", y = "m_μ(t)") +
  theme_minimal()

grid.arrange(p_mu, p_mmu, ncol = 2)
