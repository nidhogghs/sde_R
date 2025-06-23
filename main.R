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
  K = 300,           # 过程个数
  n_ave = 500,       # 每个过程的子样本数
  m = 50,            # 时间网格数
  l = 4,              # μ 的主成分个数（真实结构）
  scale_sigma = 10
)

# ======= Step 2: 估计 μ_k(t) =======
mu_res <- estimate_mu_from_data(sim_data)

# ======= Step 3: 推断 σ²_k(t) =======
sigma_res <- inference_sigma_from_data(sim_data, L = 2)

# ======= Step 4: 可视化第 k 条曲线的估计效果 =======

# 绘制 μ_k(t) 与 σ²_k(t) 的真实值与估计值对比图
#展示第k_show个过程
k_show <- 3
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
# # === Step 5: 传统方法估计 μ_k(t) 并可视化 ===
# mu_hat_traditional <- estimate_mu_traditional_from_data(sim_data)
# res_traditional <- evaluate_traditional_mu_estimation(sim_data, mu_hat_traditional, show_plot_k = c(1, 2, 3))
# # ======= End of main.R =======
