#!/bin/bash
# =========================================
# grid_run_all.sh —— 遍历 3×3×3×3 参数组合，每组 500 次仿真
# =========================================

# 参数网格
K_list=(50 300 500)
n_ave_list=(50 300 500)
k_list=(1 2 4)
l_list=(1 2 4)

CORES=32   # 并发核数，请按机器能力调整

# 关闭 citation 提示（已运行过可注释）
parallel --citation >/dev/null 2>&1

for K in "${K_list[@]}"; do
  for n_ave in "${n_ave_list[@]}"; do
    for k in "${k_list[@]}"; do
      for l in "${l_list[@]}"; do
        echo "=== Running K=$K n_ave=$n_ave k=$k l=$l ==="

        OUTDIR="output/K${K}_n${n_ave}_k${k}_l${l}"
        mkdir -p "$OUTDIR"

        # 并行 500 次
        parallel -j $CORES \
          "Rscript main_batch.R seed={} K=$K n_ave=$n_ave k=$k l=$l output_dir=$OUTDIR" \
          ::: {1..500}

        # 汇总当前组合
        Rscript summarize_results.R "$OUTDIR"

        echo "=== Done  K=$K n_ave=$n_ave k=$k l=$l ==="
      done
    done
  done
done

# 可选：最终整合全部 summary_*.csv
echo ">>> Creating master grid summary..."
Rscript - <<'EOF'
files <- list.files("output", pattern="^summary_K\\d+_n\\d+_k\\d+_l\\d+\\.csv$", full.names=TRUE)
if (length(files)>0) {
  all <- do.call(rbind, lapply(files, read.csv))
  write.csv(all, "output/grid_summary_all.csv", row.names=FALSE)
  cat("Master summary saved to output/grid_summary_all.csv\n")
}
EOF
