# ========= main_real_3.R =========
# 只替换数据来源，其余算法文件一概不动

# ------- 加载原有依赖 -------
source("config.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/mu_sigma_est.R")
source("R/traditional.R")

library(data.table)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(np)

# ------- 1. 读取三个板块 CSV（各 500 股票 × 120 天） -------
csv_paths <- c(
  "data/US_combined_close_prices.csv",
  "data/HK_combined_close_prices.csv",
  "data/CN_combined_close_prices.csv"
)

K        <- length(csv_paths)   # =3
n_each   <- 500                 # 每板块 500 只
price_ls <- list()

for (p in csv_paths) {
  # 跳过首行股票名，保留 120 行收盘价；跳日期列，取前 500 列
  dt  <- fread(p, skip = 1, header = FALSE, encoding = "UTF-8", nrows = 120)
  mat <- as.matrix(dt[, 2:(n_each + 1)])          # 120 × 500
  mat <- apply(mat, 2, as.numeric)                # 保证 numeric
  price_ls[[length(price_ls) + 1]] <- mat
}

# ------- 2. 计算 X_Delta (= diff(log(price))/sqrt(delta)) -------
m_raw <- 120            # 原始收盘价行数
m     <- m_raw - 1      # 差分后行数 (119)
Delta <- 1 / m          # 与模拟保持 T=1 的设定一致
scale_fac <- 1 / sqrt(Delta)

X_ls <- lapply(price_ls, function(mat){
  diff(log(mat)) * scale_fac      # 119 × 500
})

X_Delta <- do.call(cbind, X_ls)   # 119 × 1500
n_vec   <- rep(n_each, K)
ts      <- seq(Delta, 1, length.out = m)

# ------- 3. 构造 sim_data 按你原有字段命名 -------
sim_data <- list(
  ts     = ts,
  delta  = Delta,
  n_vec  = n_vec,
  X_Delta = X_Delta,
  mu_true = NULL,          # 无真值
  sigma2_true = NULL
)

cat(sprintf("dataloader completed:  m=%d  K=%d  n_vec=%s\n",
            m, K, paste(n_vec, collapse = "/")))

# ------- 4. 估计 μ_k(t) -------
mu_res <- estimate_mu_from_data(sim_data)

# ------- 5. 推断 σ²_k(t) -------
sigma_res <- inference_sigma_from_data(sim_data)

# ------- 6. 简单可视化 (任选一个 k) -------
k_show <- 10
p_mu <- plot_compare_single(ts, mu_res$mu_hat, mu_res$mu_hat,
                            k = k_show, label = "mu_k(t)  (est.)")
p_sigma <- plot_compare_single(ts, sigma_res$sigma2_hat, sigma_res$sigma2_hat,
                               k = k_show, label = "sigma²_k(t)  (est.)")

ci_mu <- compute_mu_ci(
  mu_hat     = mu_res$mu_hat,
  PCs        = mu_res$PCs_hat,
  sigma2_hat = sigma_res$sigma2_hat,
  ts         = ts,
  n_vec      = n_vec,
  alpha      = 0.05,
  L          = mu_res$L
)
p_ci <- plot_mu_with_ci(ts, mu_res$mu_hat, mu_res$mu_hat,
                        ci_mu$lower, ci_mu$upper, k = k_show)

m_mu_hat <- mu_res$m_mu_hat
p_m_mu <- ggplot(data.frame(t = ts, m_mu_hat),
                 aes(t, m_mu_hat)) +
  geom_line(colour = "steelblue") +
  labs(title = "Estimated Mean Drift  m_mu(t)", x = "t", y = "m_mu(t)") +
  theme_minimal()

grid.arrange(p_mu, p_sigma, p_ci, p_m_mu, ncol = 4)
# ========= END =========
