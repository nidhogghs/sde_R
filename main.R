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
  n_ave = 300,       # 每个过程的子样本数
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

# 绘制 μ_k(t) 与 σ²_k(t) 的真实值与估计值对比图
#展示第k_show个过程
k_show <- 1
p_mu <- plot_compare_single(sim_data$ts, sim_data$mu_true, mu_res$mu_hat, k = k_show, label = "mu_k(t)")
p_sigma <- plot_compare_single(sim_data$ts, sim_data$sigma2_true, sigma_res$sigma2_hat, k = k_show, label = "sigma²_k(t)")

# 并排展示两个图
# grid.arrange(p_mu, p_sigma, ncol = 2)

# 计算置信区间
ci_mu <- compute_mu_ci(
  mu_hat = mu_res$mu_hat,
  PCs = mu_res$PCs_hat,
  sigma2_hat = sigma_res$sigma2_hat,
  ts = sim_data$ts,
  n_vec = sim_data$n_vec,
  alpha = 0.05,
  L = mu_res$L
)
# 绘图
p_ci <- plot_mu_with_ci(
  ts = sim_data$ts,
  mu_true = sim_data$mu_true,
  mu_hat = mu_res$mu_hat,
  ci_lower = ci_mu$lower,
  ci_upper = ci_mu$upper,
  k = k_show
)
grid.arrange(p_mu, p_sigma, p_ci, ncol = 3)



mu_hat_traditional <- estimate_mu_traditional_from_data(sim_data)
res_traditional <- evaluate_traditional_mu_estimation(sim_data, mu_hat_traditional, show_plot_k = k_show)

# 组合 True / FDA / Traditional 三种估计进行对比
df_mu_methods <- data.frame(
  t = sim_data$ts,
  True = sim_data$mu_true[, k_show],
  FDA = mu_res$mu_hat[, k_show],
  Traditional = mu_hat_traditional[, k_show]
)
df_mu_methods_long <- reshape2::melt(df_mu_methods, id.vars = "t",
                                     variable.name = "Method", value.name = "Value")
p_mu_methods <- ggplot(df_mu_methods_long, aes(x = t, y = Value, color = Method)) +
  geom_line(linewidth = 1) +
  labs(title = sprintf("Comparison of mu_k(t) Estimates (k = %d)", k_show),
       x = "t", y = "mu_k(t)") +
  theme_minimal()

# # ======= End of main.R =======
# ===== Visualize and compare m_mu(t) =====
m_mu_true <- rowMeans(sim_data$mu_true)       # true mean function
m_mu_hat  <- mu_res$m_mu_hat                  # estimated mean

df_m_mu <- data.frame(
  t = sim_data$ts,
  True = m_mu_true,
  Estimated = m_mu_hat
)
df_m_mu_long <- reshape2::melt(df_m_mu, id.vars = "t",
                               variable.name = "Type", value.name = "Value")

p_m_mu <- ggplot(df_m_mu_long, aes(x = t, y = Value, color = Type)) +
  geom_line(linewidth = 1) +
  labs(title = "Comparison of m_mu(t): True vs Estimated",
       x = "t", y = "m_mu(t)") +
  theme_minimal()

grid.arrange(  p_ci, p_m_mu, ncol = 2)
