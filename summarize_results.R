# ========== summarize_results.R ==========
# 用法: Rscript summarize_results.R output/K500_n300_k1_l2

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "output"

files <- list.files(output_dir, pattern = "^result_.*\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop(sprintf("没有找到任何 result_*.csv 文件，请检查目录: %s", output_dir))
}

# 合并结果
results <- do.call(rbind, lapply(files, read.csv))

# 保存完整记录
write.csv(results, file.path(output_dir, "all_results.csv"), row.names = FALSE)

# 提取参数
params <- results[1, c("K", "n_ave", "k", "l")]

# 计算平均误差
summary_avg <- data.frame(
  K = params$K,
  n_ave = params$n_ave,
  k = params$k,
  l = params$l,
  avg_rmse_mu     = mean(results$rmse_mu),
  avg_rmse_m_mu   = mean(results$rmse_m_mu),
  avg_rmse_G_mu   = mean(results$rmse_G_mu),
  avg_rmse_sigma2 = mean(results$rmse_sigma2)
)

# 写入汇总文件到统一目录（output/summary_K*_n*_k*_l*.csv）
summary_file <- sprintf("output/summary_K%d_n%d_k%d_l%d.csv",
                        summary_avg$K, summary_avg$n_ave,
                        summary_avg$k, summary_avg$l)

write.csv(summary_avg, summary_file, row.names = FALSE)

cat("汇总完成 =>", summary_file, "\n")
