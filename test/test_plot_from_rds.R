# test_plot_from_rds.R
# ========= 从 RDS 文件加载数据并绘图，支持参数化输出格式 ==========

library(ggplot2)

# === 读取命令行参数 ===
args <- commandArgs(trailingOnly = TRUE)
args_list <- list()
for (arg in args) {
  keyval <- strsplit(arg, "=")[[1]]
  if (length(keyval) == 2) {
    key <- keyval[1]
    val <- keyval[2]
    args_list[[key]] <- val
  }
}

# 默认输入文件路径
file_path <- if (!is.null(args_list$file)) args_list$file else "output/result_latest.rds"
# 默认输出格式为 png
fmt <- if (!is.null(args_list$fmt)) args_list$fmt else "png"

# === 加载数据 ===
result <- readRDS(file_path)

mu <- result$mu
mu_hat <- result$mu_hat
m_mu_true <- result$m_mu_true
m_mu_hat <- result$m_mu_hat
PCs_true <- result$PCs_true
PCs_hat <- result$PCs_hat
lams_true <- result$lams_true
lams_hat <- result$lams_hat
ts <- result$ts
K <- if (!is.null(result$K)) result$K else NA
n_ave <- if (!is.null(result$n_ave)) result$n_ave else NA
L <- if (!is.null(result$L)) result$L else NA

prefix <- sprintf("K%d_n%d_L%d", K, n_ave, L)
dir.create("figures", showWarnings = FALSE)

# 构造文件名工具
outname <- function(name) {
  sprintf("figures/plot_%s_%s.%s", name, prefix, fmt)
}

# === 图像保存 ===
dpi_val <- if (fmt %in% c("png", "jpeg", "jpg", "tiff")) 300 else NA

p1 <- ggplot(data.frame(t = ts, True = m_mu_true, Estimated = m_mu_hat), aes(x = t)) +
  geom_line(aes(y = True, color = "True")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "Mean Function m_mu(t)", y = "Value") +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  theme_minimal()
ggsave(outname("m_mu"), p1, width = 6, height = 4, dpi = dpi_val)

p2 <- ggplot(data.frame(t = ts,
                        True = PCs_true[,1] * sign(PCs_true[1,1]),
                        Estimated = PCs_hat[,1] * sign(PCs_hat[1,1])),
             aes(x = t)) +
  geom_line(aes(y = True, color = "True")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "First Principal Component phi_1(t)", y = "Value") +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  theme_minimal()
ggsave(outname("pc1"), p2, width = 6, height = 4, dpi = dpi_val)

p3 <- ggplot(data.frame(t = ts, True = mu[-1,1], Estimated = mu_hat[,1]), aes(x = t)) +
  geom_line(aes(y = True, color = "True")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "First Sample mu_1(t)", y = "Value") +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  theme_minimal()
ggsave(outname("mu1"), p3, width = 6, height = 4, dpi = dpi_val)

p4 <- ggplot(data.frame(Index = 1:length(lams_true), True = lams_true, Estimated = lams_hat), aes(x = Index)) +
  geom_point(aes(y = True, color = "True")) +
  geom_point(aes(y = Estimated, color = "Estimated")) +
  labs(title = "Eigenvalue Spectrum lambda_j", y = "Value") +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  theme_minimal()
ggsave(outname("eigen"), p4, width = 6, height = 4, dpi = dpi_val)

# # 保存为 PNG
# Rscript test_plot_from_rds.R file=output/result_K100_n100_L2.rds fmt=png

# # 保存为 JPEG
# Rscript test_plot_from_rds.R file=output/result_K100_n100_L2.rds fmt=jpeg

# # 保存为 SVG（可用于网页或矢量编辑）
# Rscript test_plot_from_rds.R file=output/result_K100_n100_L2.rds fmt=svg
