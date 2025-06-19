使用说明
对于μ的估计：
在R环境中使用source("main.R")加载μ的估计部分主程序，使用run_main()指令运行(此处可指定K=xxx,n_ave=xxx,目前模型仅支持所有的k轨迹数n_ave相等，等合并跑通了再改成可不同)。
对于σ的估计：
在R环境中使用source("inference.R")加载并运行sigma的估计部分。
模块化与合并：
目前正在对μ和σ部分的代码解耦并合并，但是目前遇到的问题是。如果先生成指数，则μ的估计会低超过两个数量级(e-3到5e-1)，几乎失效。直接生成sigma,对sigma的估计很差。
模块化测试需加载脚本 source("testonerun")  随后运行
运行test_mu_estimation() 先生成指数 对μ的估计
运行test_mu_estimation2() 直接生成sigma,对μ的估计
运行test_mu_estimation2() 先生成指数 对sigma的估计

