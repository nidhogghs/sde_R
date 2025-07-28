#!/bin/bash

# 参数网格
K_list=(50 300 500)
n_ave_list=(50 300 500)
k_list=(1 2 4)
l_list=(1 2 4)

# 资源限制
CORES=32   # 可调并发数

# 主循环
for K in "${K_list[@]}"; do
  for n_ave in "${n_ave_list[@]}"; do
    for k in "${k_list[@]}"; do
      for l in "${l_list[@]}"; do

        echo "=== Running K=$K n_ave=$n_ave k=$k l=$l ==="

        OUTDIR="output/K${K}_n${n_ave}_k${k}_l${l}"
        mkdir -p "$OUTDIR"

        # 并行执行 500 次仿真
        parallel -j $CORES "Rscript main_batch.R seed={} K=$K n_ave=$n_ave k=$k l=$l > $OUTDIR/log_seed_{}.txt 2>&1" ::: {1..500}

        # 汇总当前组结果
        Rscript summarize_results.R "$OUTDIR"

        echo "=== Done K=$K n_ave=$n_ave k=$k l=$l ==="
      done
    done
  done
done
