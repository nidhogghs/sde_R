# config.R
# ========= 参数配置文件 ==========

# 时间范围
T <- 1            # 观察区间长度 [0, T]
X0 <- 1           # 初始值 log(X0) 起始点

# sigma(t) 模拟参数
a <- 1
alpha <- 1
k <- 2

# mu(t) 模拟参数
b <- 1.5
beta <- 1
l <- 2

# Default simulation grids
K_grid <- c(100, 300, 500)        # 样本规模列表
n_ave_grid <- c(100, 300)         # 每条曲线的观测数列表
m_default <- 50                   # 默认时间网格数
Nsim <- 10                        # 仿真重复次数
h_min <- 0.021                    # 核带宽下限

# 最小核带宽（用于二维协方差估计）
h_min <- 0.021
