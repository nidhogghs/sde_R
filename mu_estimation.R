# ========= 模拟主程序入口 =========

source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/one_run.R")

library("np")
library("plotly")
dir.create("output", showWarnings = FALSE, recursive = TRUE)

run_main <- function(K = NULL, n_ave = NULL, Nsim_local = NULL, m = NULL) {
  if (is.null(K)) K <- K_grid[1]
  if (is.null(n_ave)) n_ave <- n_ave_grid[1]
  if (is.null(Nsim_local)) Nsim_local <- Nsim
  if (is.null(m)) m <- m_default

  res <- list()
  cat("Running K =", K, "n_ave =", n_ave, "\n")
  res_N <- sapply(1:Nsim_local, function(sd) {
    cat("  Simulation:", sd, "\n")
    one_run(K, n_ave, m, sd)
  })
  res_mean <- rowMeans(res_N)
  res[[paste0("K", K, "_n", n_ave)]] <- res_mean

  res_df <- do.call(cbind, res)
  rownames(res_df) <- c("err_m", "err_G", "err_lams", "err_PCs", "err_mu")
  print("=== 平均误差汇总 ===")
  print(round(res_df, 4))
}


# ========== 自动执行逻辑 ==========

# 从命令行提取参数（Rscript 方式）
args <- commandArgs(trailingOnly = TRUE)
args_list <- list()
if (length(args) > 0) {
  for (arg in args) {
    kv <- strsplit(arg, "=")[[1]]
    key <- kv[1]
    value <- as.numeric(kv[2])
    args_list[[key]] <- value
  }
}


# 判断当前运行环境
if (sys.nframe() == 0) {
  if (length(args_list) > 0) {
    # 有命令行参数时，只运行一个组合
    do.call(run_main, args_list)
  } else {
    # 无参数时，枚举所有参数组合
    for (K_val in K_grid) {
      for (n_val in n_ave_grid) {
        run_main(K = K_val, n_ave = n_val)
      }
    }
  }
}

