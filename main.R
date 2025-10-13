# ======= 主程序入口 main.R =======

# 载入配置与函数模块
source("config.R")               # 模型参数配置
source("R/simulation.R")         # 数据生成函数
source("R/estimation.R")         # μ 和 G_μ 的估计
source("R/recovery.R")           # 主成分重建函数
source("R/utils.R")              # 工具函数：如 cluster_mean, compute_qk_mc 等
source("R/mu_sigma_est.R")       # 封装的 estimate_mu_from_data / inference_sigma_from_data
source("R/traditional.R")
# 加载绘图和核回归相关库
library(gridExtra)               # 用于组合多个 ggplot 图
library(np)                      # 非参数核平滑包
library(ggplot2)                 # 基础 ggplot 绘图
library(reshape2)                # 数据长宽转换

# ======= Step 1: 模拟数据 =======
sim_data <- prepare_simulation_data(
  K = 500,           # 过程个数
  n_ave = 30,       # 每个过程的子样本数
  m = 50,            # 时间网格数
  k = 4,
  l = 4,              # μ 的主成分个数（真实结构）
  scale_sigma = 5
)

# ======= Step 2: 估计 μ_k(t) =======
mu_res <- estimate_mu_from_data(sim_data)
# ======= Step 2.5: 输出 μ 的 RMSE =======
rmse_mu <- sqrt(mean((mu_res$mu_hat - sim_data$mu_true)^2))
cat(sprintf("RMSE(mu): %.6f\n", rmse_mu))


# ======= Step 3: 推断 σ²_k(t) =======
sigma_res <- inference_sigma_from_data(sim_data)

# ======= Step 4: 可视化第 k 条曲线的估计效果 =======

# ==============================
# 可运行版：单独展示并保存所有图
# ==============================

# 创建输出目录（若不存在则创建）
if (!dir.exists("output")) dir.create("output", recursive = TRUE)

# --- 绘制 μ_k(t) 与 σ²_k(t) ---
k_show <- 1
p_mu <- plot_compare_single(sim_data$ts, sim_data$mu_true, mu_res$mu_hat,
                            k = k_show, label = "mu_i(t)")
p_sigma <- plot_compare_single(sim_data$ts, sim_data$sigma2_true, sigma_res$sigma2_hat,
                               k = k_show, label = "sigma²_i(t)")

# 展示并保存
print(p_mu)
ggsave(sprintf("output/mu_k%d.png", k_show), p_mu, width = 6, height = 4, dpi = 300)

print(p_sigma)
ggsave(sprintf("output/sigma2_k%d.png", k_show), p_sigma, width = 6, height = 4, dpi = 300)

# --- 计算置信区间 ---
ci_mu <- compute_mu_ci(
  mu_hat = mu_res$mu_hat,
  PCs = mu_res$PCs_hat,
  sigma2_hat = sigma_res$sigma2_hat,
  ts = sim_data$ts,
  n_vec = sim_data$n_vec,
  alpha = 0.05,
  L = mu_res$L
)

p_ci <- plot_mu_with_ci(
  ts = sim_data$ts,
  mu_true = sim_data$mu_true,
  mu_hat = mu_res$mu_hat,
  ci_lower = ci_mu$lower,
  ci_upper = ci_mu$upper,
  k = k_show
)
print(p_ci)
ggsave(sprintf("output/mu_ci_k%d.png", k_show), p_ci, width = 6, height = 4, dpi = 300)

# --- 传统估计法 ---
mu_hat_traditional <- estimate_mu_traditional_from_data(sim_data)
res_traditional <- evaluate_traditional_mu_estimation(sim_data, mu_hat_traditional, show_plot_k = k_show)
# --- 三方法对比（黑白线型区分） ---
#（保持 df_mu_methods_long 不变）
# 为了图例顺序稳定：True -> FDA -> Traditional
df_mu_methods <- data.frame(
  t = sim_data$ts,
  True = sim_data$mu_true[, k_show],
  FDA = mu_res$mu_hat[, k_show],
  Traditional = mu_hat_traditional[, k_show]
)
df_mu_methods_long <- melt(df_mu_methods, id.vars = "t",
                           variable.name = "Method", value.name = "Value")
df_mu_methods_long$Method <- factor(df_mu_methods_long$Method,
                                    levels = c("True", "FDA", "Traditional"))

p_mu_methods <- ggplot(df_mu_methods_long, aes(x = t, y = Value, linetype = Method)) +
  geom_line(linewidth = 1, color = "black") +
  scale_linetype_manual(
    name = "Legend",
    values = c("True" = "solid", "FDA" = "longdash", "Traditional" = "dotdash")
  ) +
  labs(x = "t", y = expression(mu[i](t))) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()  # 去掉标题
  )

print(p_mu_methods)
ggsave(sprintf("output/mu_methods_i%d.png", k_show), p_mu_methods, width = 6, height = 4, dpi = 300)

# p_mu_methods <- ggplot(df_mu_methods_long, aes(x = t, y = Value, color = Method)) +
#   geom_line(linewidth = 1) +
#   labs(title = sprintf("Comparison of mu_k(t) Estimates (k = %d)", k_show),
#        x = "t", y = "mu_k(t)") +
#   theme_minimal()

# print(p_mu_methods)
# ggsave(sprintf("output/mu_methods_k%d.png", k_show), p_mu_methods, width = 6, height = 4, dpi = 300)

# --- Population mean m_mu(t) ---
m_mu_true <- rowMeans(sim_data$mu_true)
m_mu_hat  <- mu_res$m_mu_hat

df_m_mu <- data.frame(
  t = sim_data$ts,
  True = m_mu_true,
  Estimated = m_mu_hat
)
df_m_mu_long <- melt(df_m_mu, id.vars = "t",
                     variable.name = "Type", value.name = "Value")

# p_m_mu <- ggplot(df_m_mu_long, aes(x = t, y = Value, color = Type)) +
#   geom_line(linewidth = 1) +
#   labs(title = "Comparison of m_mu(t): True vs Estimated",
#        x = "t", y = "m_mu(t)") +
#   theme_minimal()

# print(p_m_mu)
# ggsave("output/m_mu_comparison.png", p_m_mu, width = 6, height = 4, dpi = 300)
# --- Population mean m_mu(t)（黑白线型区分） ---
#（保持 df_m_mu_long 不变）
# 统一顺序：True -> Estimated
df_m_mu_long$Type <- factor(df_m_mu_long$Type, levels = c("True", "Estimated"))

p_m_mu <- ggplot(df_m_mu_long, aes(x = t, y = Value, linetype = Type)) +
  geom_line(linewidth = 1, color = "black") +
  scale_linetype_manual(
    name = "Legend",
    values = c("True" = "solid", "Estimated" = "longdash")
  ) +
  labs(x = "t", y = expression(m[mu](t))) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()  # 去掉标题
  )

print(p_m_mu)
ggsave("output/m_mu_comparison.png", p_m_mu, width = 6, height = 4, dpi = 300)

# --- 同时生成一个多页 PDF ---
pdf("output/all_plots.pdf", width = 7, height = 5)
print(p_mu)
print(p_sigma)
print(p_ci)
print(p_mu_methods)
print(p_m_mu)
dev.off()

cat("[Saved] 所有图像已输出至 output/ 目录下\n")
