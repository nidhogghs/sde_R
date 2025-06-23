# Functional Estimation for Drift and Diffusion in SDE

本项目基于函数型数据分析（FDA）框架，实现对随机微分方程中漂移项 $\mu_k(t)$ 与扩散项 $\sigma_k^2(t)$ 的非参数估计。项目支持完整的模拟、估计、置信区间推断与可视化流程，并提供传统估计方法对比。

---

## 使用方式

### 1. 运行主流程

在 R 环境中执行：

```r
source("main.R")
```

这将依次完成以下任务：

- 生成模拟数据；
- 估计漂移项 μ 与扩散项 σ²；
- 构造 μ 的置信区间；
- 可视化 μ_k(t)、σ²_k(t) 估计与 μ 的置信区间。

默认参数：K = 300，n_ave = 500，m = 50，L = 2。

---

## 模块说明

### μ 与 σ² 的估计函数入口

两者均封装于 `R/mu_sigma_est.R`，使用如下函数：

```r
mu_res <- estimate_mu_from_data(sim_data, L = 2)
sigma_res <- inference_sigma_from_data(sim_data, L = 2)
```

其中：

- `mu_res$mu_hat`：估计的 μ_k(t)
- `sigma_res$sigma2_hat`：估计的 σ²_k(t)

---

## μ 的置信区间构造

根据函数型主成分重构后的置信区间理论，估计值具有如下性质：

$$
\sqrt{n_k}(\hat{\mu}_k(t) - \mu_k(t)) \overset{d}{\to} \mathcal{N}(0, V_\mu(t))
$$

基于此，置信区间表达为：

$$
\hat{\mu}_k(t) \pm z_{1 - \alpha/2} \sqrt{V_\mu(t) / n_k}
$$

调用接口如下：

```r
ci_mu <- compute_mu_ci(
  mu_hat = mu_res$mu_hat,
  PCs = mu_res$PCs_hat,
  sigma2_hat = sigma_res$sigma2_hat,
  ts = sim_data$ts,
  n_vec = sim_data$n_vec,
  alpha = 0.05,
  L = mu_res$L
)
```

---

## 可视化模块

图像均为 ggplot 生成，可使用如下方式绘制并排列展示：

```r
p_mu <- plot_compare_single(...)
p_sigma <- plot_compare_single(...)
p_ci <- plot_mu_with_ci(...)
grid.arrange(p_mu, p_sigma, p_ci, ncol = 3)
```

---

## 传统方法对比

对每个过程 $k$，取平均轨迹后核平滑估计：

```r
res_traditional <- evaluate_traditional_mu_estimation(
  sim_data, 
  mu_hat = estimate_mu_traditional_from_data(sim_data), 
  show_plot_k = c(1, 2, 3)
)
```

输出包括每个过程 RMSE 与平均误差。

---

## 误差精度说明

- 若先生成 $\log \sigma^2$ 再估计 μ，可得误差 $O(10^{-3})$；
- 若直接生成 σ 估计 μ，误差升至 $O(10^{-1}) \sim O(1)$；
- 理论上更推荐建模 $\log \sigma^2$ 后取指数；
- FDA 估计融合整体协方差结构，对单个过程能显著提升鲁棒性。

---

## 文件结构总览

| 文件                | 功能描述                                       |
|---------------------|-----------------------------------------------|
| `main.R`            | 主控流程入口：仿真、估计、CI、绘图            |
| `config.R`          | 参数配置，包括常数、主成分个数、核带宽等       |
| `simulation.R`      | 生成 μ 和 log σ² 的轨迹及观测数据              |
| `estimation.R`      | FDA 方法估计 μ 与协方差结构                    |
| `recovery.R`        | 主成分重建个体 μ_k(t) 或 V_k(t)                |
| `mu_sigma_est.R`    | μ 与 σ² 的估计与置信区间封装模块              |
| `traditional.R`     | 传统方法估计 μ 与误差分析                      |
| `utils.R`           | 工具函数：如 cluster_mean, compute_qk_mc 等   |
| `testonerun.R`      | 各测试函数封装入口                             |

