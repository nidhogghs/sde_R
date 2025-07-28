# ======= main_batch.R =======

# 载入函数模块
source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/mu_sigma_est.R")

# ---------- Step 1: 解析命令行参数 ----------
args <- commandArgs(trailingOnly = TRUE)
params <- list(seed = NA, K = 1000, n_ave = 1, k = 2, l = 2)

for (arg in args) {
  kv <- strsplit(arg, "=")[[1]]
  key <- kv[1]
  val <- as.numeric(kv[2])
  if (key %in% names(params)) {
    params[[key]] <- val
  }
}

# 检查种子
if (is.na(params$seed)) stop("You must provide seed=xxx")

# ---------- Step 2: 生成数据 ----------
set.seed(params$seed)

sim_data <- prepare_simulation_data(
  K = params$K,
  n_ave = params$n_ave,
  m = 50,
  k = params$k,
  l = params$l,
  scale_sigma = 5,
  seed_sigma = 311,
  seed_mu = 311,
  seed_X = 69
)

# ---------- Step 3: 估计 μ 和 σ² ----------
mu_res <- estimate_mu_from_data(sim_data)
sigma_res <- inference_sigma_from_data(sim_data)

# ---------- Step 4: 误差统计 ----------
rmse_mu      <- sqrt(mean((mu_res$mu_hat - sim_data$mu_true)^2))
rmse_m_mu    <- sqrt(mean((mu_res$m_mu_hat - rowMeans(sim_data$mu_true))^2))
rmse_G_mu    <- sqrt(mean((mu_res$G_mu_hat - cov(t(sim_data$mu_true)))^2))
rmse_sigma2  <- sqrt(mean((sigma_res$sigma2_hat - sim_data$sigma2_true)^2))

# ---------- Step 5: 输出结果 ----------
res <- data.frame(
  seed = params$seed,
  K = params$K,
  n_ave = params$n_ave,Rscript summarize_results.R
  k = params$k,
  l = params$l,
  rmse_mu = rmse_mu,
  rmse_m_mu = rmse_m_mu,
  rmse_G_mu = rmse_G_mu,
  rmse_sigma2 = rmse_sigma2
)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
output_file <- sprintf("output/result_seed%03d_K%d_n%d_k%d_l%d.csv",
                       params$seed, params$K, params$n_ave, params$k, params$l)

write.csv(res, file = output_file, row.names = FALSE)
