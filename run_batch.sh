#!/bin/bash

mkdir -p results

K=300
n_ave=500
m=50
l=4
k=4
Nsim=1  # 每个子任务仅运行1次（总共运行500次）

for seed_id in $(seq 1 500); do
  echo "Launching task for seed $seed_id"
  Rscript run_single_experiment.R $K $n_ave $m $seed_id $Nsim $l $k &
done

wait
echo "All 500 jobs completed."
