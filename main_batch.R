# =========================================
# main_batch.R  ——  单次仿真+误差输出
# =========================================
# 用法示例：
#   Rscript main_batch.R seed=1 K=50 n_ave=50 k=1 l=1 \
#                      output_dir=output/K50_n50_k1_l1

# ---------- 载入代码模块 ----------
source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/mu_sigma_est.R")

# ---------- 解析命令行参数 ----------
args   <- commandArgs(trailingOnly = TRUE)
params <- list(seed = NA, K = 1000, n_ave = 1, k = 2, l = 2)
output_dir <- "output"   # 缺省目录

for (arg in args) {
  kv <- strsplit(arg, "=")[[1]]
  if (length(kv) != 2) next
  key <- kv[1]; val <- kv[2]
  if (key %in% names(params)) {
    params[[key]] <- as.numeric(val)
  } else if (key == "output_dir") {
    output_dir <- val
  }
}

if (is.na(params$seed)) stop("Please provide seed=xxx")

# ---------- 数据生成 ----------
set.seed(params$seed)
sim_data <- prepare_simulation_data(
  K           = params$K,
  n_ave       = params$n_ave,
  m           = 50,
  k           = params$k,
  l           = params$l,
  scale_sigma = 5,
  seed_sigma  = 311,
  seed_mu     = 311,
  seed_X      = 69
)

# ---------- 估计 ----------
mu_res    <- estimate_mu_from_data(sim_data)
sigma_res <- inference_sigma_from_data(sim_data)

# ---------- 误差评估 ----------
rmse_mu      <- sqrt(mean((mu_res$mu_hat - sim_data$mu_true)^2))
rmse_m_mu    <- sqrt(mean((mu_res$m_mu_hat - rowMeans(sim_data$mu_true))^2))
rmse_G_mu    <- sqrt(mean((mu_res$G_mu_hat - cov(t(sim_data$mu_true)))^2))
rmse_sigma2  <- sqrt(mean((sigma_res$sigma2_hat - sim_data$sigma2_true)^2))

# ---------- 保存结果 ----------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

res <- data.frame(
  seed         = params$seed,
  K            = params$K,
  n_ave        = params$n_ave,
  k            = params$k,
  l            = params$l,
  rmse_mu      = rmse_mu,
  rmse_m_mu    = rmse_m_mu,
  rmse_G_mu    = rmse_G_mu,
  rmse_sigma2  = rmse_sigma2
)

outfile <- sprintf("%s/result_seed%03d_K%d_n%d_k%d_l%d.csv",
                   output_dir, params$seed,
                   params$K, params$n_ave, params$k, params$l)
write.csv(res, outfile, row.names = FALSE)
