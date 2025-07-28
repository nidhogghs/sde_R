#!/bin/bash

# 参数（可调整）
K=500
n_ave=500
k=4
l=4
CORES=32

# 创建输出目录
mkdir -p output

# 并行仿真
echo ">>> Running simulations..."
parallel -j $CORES "Rscript main_batch.R seed={} K=$K n_ave=$n_ave k=$k l=$l" ::: {1..500}

# 汇总结果
echo ">>> Summarizing results..."
Rscript summarize_results.R

echo ">>> Done. Output in output/ directory."
