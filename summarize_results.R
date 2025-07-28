# summarize_results.R
# 汇总所有 result_*.csv 文件，并计算平均误差

files <- list.files("output", pattern = "^result_.*\\.csv$", full.names = TRUE)

# 合并所有结果
results <- do.call(rbind, lapply(files, read.csv))

# 保存完整汇总表
write.csv(results, "output/all_results_summary.csv", row.names = FALSE)

# 计算平均值
summary_avg <- data.frame(
  avg_rmse_mu      = mean(results$rmse_mu),
  avg_rmse_m_mu    = mean(results$rmse_m_mu),
  avg_rmse_G_mu    = mean(results$rmse_G_mu),
  avg_rmse_sigma2  = mean(results$rmse_sigma2)
)

write.csv(summary_avg, "output/summary_average.csv", row.names = FALSE)
cat("汇总完成，共", nrow(results), "次仿真。\n")
print(summary_avg)
