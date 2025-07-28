# ========== summarize_results.R ==========

files <- list.files("output", pattern = "^result_.*\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop("没有找到任何 result_*.csv 文件，请检查是否运行成功！")
}

# 读取所有结果
results <- do.call(rbind, lapply(files, read.csv))

# 保存完整结果表
write.csv(results, "output/all_results_summary.csv", row.names = FALSE)

# 取首行参数（假设所有实验使用相同参数）
param_values <- results[1, c("K", "n_ave", "k", "l")]

# 计算平均值
summary_avg <- data.frame(
  K        = param_values$K,
  n_ave    = param_values$n_ave,
  k        = param_values$k,
  l        = param_values$l,
  avg_rmse_mu     = mean(results$rmse_mu),
  avg_rmse_m_mu   = mean(results$rmse_m_mu),
  avg_rmse_G_mu   = mean(results$rmse_G_mu),
  avg_rmse_sigma2 = mean(results$rmse_sigma2)
)

# 根据参数生成输出文件名，防止覆盖
summary_file <- sprintf("output/summary_K%d_n%d_k%d_l%d.csv",
                        summary_avg$K, summary_avg$n_ave,
                        summary_avg$k, summary_avg$l)

write.csv(summary_avg, summary_file, row.names = FALSE)

cat("实验汇总完成，写入:", summary_file, "\n")
print(summary_avg)
