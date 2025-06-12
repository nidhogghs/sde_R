# config.R
# ========= 参数配置文件 ==========

# 时间范围
T <- 1            # 观察区间长度 [0, T]
X0 <- 1           # 初始值 log(X0) 起始点

# sigma(t) 模拟参数
a <- 1
alpha <- 1
k <- 4

# mu(t) 模拟参数
b <- 1.5
beta <- 1
l <- 4

# 默认仿真设置
K_default <- 500
n_ave_default <- 300
m_default <- 50
Nsim <- 10

# 最小核带宽（用于二维协方差估计）
h_min <- 0.021
