# main.R
# ========= 模拟主程序入口 =========

source("config.R")
source("R/simulation.R")
source("R/estimation.R")
source("R/recovery.R")
source("R/utils.R")
source("R/one_run.R")

library("np")
library("plotly")

# 单次运行示例
res_example <- one_run(K_default, n_ave_default, m_default, sd = 1)

# 多次仿真，统计误差均值
K_s <- K_grid
n_aves <- n_ave_grid

res <- list()
for (K in K_s) {
  for (n_ave in n_aves) {
    cat("Running K =", K, "n_ave =", n_ave, "\n")
    res_N <- sapply(1:Nsim, function(sd) {
      cat("  Simulation:", sd, "\n")
      one_run(K, n_ave, m_default, sd)
    })
    res_mean <- rowMeans(res_N)
    res[[paste0("K", K, "_n", n_ave)]] <- res_mean
  }
}

res_df <- do.call(cbind, res)
print("=== 平均误差汇总 ===")
print(round(res_df, 4))
