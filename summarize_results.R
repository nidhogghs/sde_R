# =========================================
# summarize_results.R —— 汇总单组合 500 次结果
# =========================================
# 用法示例：
#   Rscript summarize_results.R output/K50_n50_k1_l1

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "output"

files <- list.files(output_dir, pattern = "^result_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop(sprintf("没有找到任何 result_*.csv 文件，请检查目录: %s", output_dir))
}

# 合并
results <- do.call(rbind, lapply(files, read.csv))
write.csv(results, file.path(output_dir, "all_results.csv"), row.names = FALSE)

# 取参数
params <- results[1, c("K", "n_ave", "k", "l")]

summary_avg <- data.frame(
  K               = params$K,
  n_ave           = params$n_ave,
  k               = params$k,
  l               = params$l,
  avg_rmse_mu     = mean(results$rmse_mu),
  avg_rmse_m_mu   = mean(results$rmse_m_mu),
  avg_rmse_G_mu   = mean(results$rmse_G_mu),
  avg_rmse_sigma2 = mean(results$rmse_sigma2)
)

summary_file <- sprintf("output/summary_K%d_n%d_k%d_l%d.csv",
                        params$K, params$n_ave, params$k, params$l)
write.csv(summary_avg, summary_file, row.names = FALSE)
cat("汇总完成 =>", summary_file, "\n")
