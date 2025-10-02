#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_sector_means.py — visualize averaged ΔlogP / sqrt(Δ) trajectories by sector.

Usage examples:
  python plot_sector_means.py data_dir=merged out_dir=output/sector_means
  python plot_sector_means.py data_dir=merged out_dir=output/test_sector_means \
      row_start=10 row_end=50 col_start=1 col_end=100 per_sector=1

Notes:
- row_start / row_end: 基于“原始价格表”的 1-based 行号；映射到 ΔlogP 索引为 ts_idx = raw - 1。
- col_start / col_end: 
    * 默认 (per_sector=0)：全局拼接矩阵上的列号。
    * per_sector=1：每个行业内部的局部列号。
- 输出包括：CSV（sector_means_rows_*_cols_*.csv）和 PNG 图。
"""

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------- CLI parsing ------------------------------------
def parse_kv_args(argv):
    par = dict(
        data_dir="merged",
        out_dir="output/sector_means",
        row_start=None,
        row_end=None,
        col_start=None,
        col_end=None,
        T=1.0,
        per_sector=0,   # 0=全局拼接列号，1=每个行业局部列号
    )
    for a in argv:
        if "=" not in a:
            continue
        key, val = a.split("=", 1)
        if key not in par:
            continue
        if key in ("row_start", "row_end", "col_start", "col_end"):
            par[key] = int(val) if val.strip() != "" else None
        elif key in ("T",):
            par[key] = float(val)
        elif key == "per_sector":
            par[key] = int(val)
        else:
            par[key] = val
    return par

# ----------------------------- I/O helpers ------------------------------------
def read_sector_csv(fp):
    df = pd.read_csv(fp)
    date_col = None
    for c in df.columns:
        if str(c).lower() == "date":
            date_col = c
            break
    if date_col is None:
        raise ValueError(f"No 'Date' column in file: {fp}")
    df[date_col] = pd.to_datetime(df[date_col])
    df = df.sort_values(by=date_col).set_index(date_col)
    price_df = df.select_dtypes(include=[np.number]).copy()
    price_df = price_df.dropna(axis=1, how="all")
    if price_df.shape[1] == 0:
        raise ValueError(f"No numeric price columns in: {fp}")
    return price_df.index, price_df

def to_log_diff(price_df):
    price = price_df.copy()
    price[price <= 0] = np.nan
    logp = np.log(price)
    dlogp = logp.diff().dropna(axis=0, how="all")
    return dlogp

# ----------------------------- Core pipeline ----------------------------------
def load_folder_as_blocks(data_dir):
    sector_names, n_vec, date_list, dlogp_list = [], [], [], []
    for fname in sorted(os.listdir(data_dir)):
        if not fname.lower().endswith(".csv"):
            continue
        sector = os.path.splitext(fname)[0]
        fp = os.path.join(data_dir, fname)
        dates, price_df = read_sector_csv(fp)
        dlogp = to_log_diff(price_df)
        sector_names.append(sector)
        n_vec.append(dlogp.shape[1])
        date_list.append(dlogp.index)
        dlogp_list.append(dlogp)
    if not sector_names:
        raise RuntimeError(f"No CSV found in data_dir: {data_dir}")
    return sector_names, n_vec, date_list, dlogp_list

def intersect_dates(date_list):
    common = None
    for idx in date_list:
        common = idx if common is None else common.intersection(idx)
    if common is None or len(common) == 0:
        raise RuntimeError("No common dates across sectors.")
    return common.sort_values()

def align_blocks_to_common_dates(dlogp_list, common_dates):
    return [df.reindex(index=common_dates).sort_index() for df in dlogp_list]

def build_concatenated_matrix(aligned_list):
    return pd.concat(aligned_list, axis=1)

def restrict_rows_delta_index(dlogp_concat, row_start_raw, row_end_raw):
    m_total = dlogp_concat.shape[0]
    used = [1, m_total]
    def convert(idx_raw, default_val, ub):
        if idx_raw is None: return default_val
        ts_idx = idx_raw - 1
        ts_idx = max(1, min(int(ts_idx), ub))
        return ts_idx
    used[0] = convert(row_start_raw, used[0], m_total)
    used[1] = convert(row_end_raw,   used[1], m_total)
    if used[0] > used[1]:
        raise ValueError(f"Invalid row restriction after conversion: start({used[0]}) > end({used[1]}).")
    if used != [1, m_total]:
        dlogp_concat = dlogp_concat.iloc[used[0]-1:used[1]]
        print(f"[Info] Using ΔlogP rows {used[0]}-{used[1]} (ts indices) out of {m_total}.")
    else:
        print(f"[Info] Using all ΔlogP rows (1-{m_total}).")
    return dlogp_concat, (used[0], used[1])

def restrict_cols_global_blocks(dlogp_concat, n_vec, sector_names, col_start, col_end):
    total_cols = int(sum(n_vec))
    used = [1, total_cols]
    def clamp(idx, default_val, lo, hi):
        if idx is None: return default_val
        return max(lo, min(int(idx), hi))
    used[0] = clamp(col_start, used[0], 1, total_cols)
    used[1] = clamp(col_end,   used[1], 1, total_cols)
    if used[0] > used[1]:
        raise ValueError(f"Invalid column restriction after clamping: start({used[0]}) > end({used[1]}).")

    if used == [1, total_cols]:
        print(f"[Info] Using all stock columns (1-{total_cols}) across {len(n_vec)} sectors.")
        return dlogp_concat, n_vec, sector_names, (used[0], used[1])

    new_blocks, new_n_vec, new_names, kept_global = [], [], [], []
    cursor = 0
    for nk, name in zip(n_vec, sector_names):
        block_range = range(cursor + 1, cursor + nk + 1)  # 1-based
        cursor += nk
        keep_global = [g for g in block_range if used[0] <= g <= used[1]]
        if not keep_global: continue
        keep_0 = [g - 1 for g in keep_global]
        new_blocks.append(dlogp_concat.iloc[:, keep_0])
        new_n_vec.append(len(keep_0))
        new_names.append(name)
        kept_global.extend(keep_global)

    if not new_blocks:
        raise ValueError(f"Column restriction [{used[0]}, {used[1]}] removed all stock series.")

    dlogp_new = pd.concat(new_blocks, axis=1)
    real_span = (min(kept_global), max(kept_global))
    print(f"[Info] Using stock columns {real_span[0]}-{real_span[1]} (global indices). "
          f"Retained {len(new_n_vec)} sectors / {sum(new_n_vec)} stocks.")
    return dlogp_new, new_n_vec, new_names, real_span

def restrict_cols_per_sector(dlogp_concat, n_vec, sector_names, per_start, per_end):
    if per_start is None and per_end is None:
        print("[Info] per-sector col slice not set; keeping all columns.")
        return dlogp_concat, n_vec, sector_names, (None, None)

    new_blocks, new_n_vec, new_names = [], [], []
    kept_total = 0
    cursor = 0
    for nk, name in zip(n_vec, sector_names):
        if nk <= 0:
            cursor += nk
            continue
        s = 1 if per_start is None else max(1, min(per_start, nk))
        e = nk if per_end is None else max(s, min(per_end, nk))
        if e < s:
            cursor += nk
            continue
        g0 = cursor + (s - 1)
        g1 = cursor + (e - 1)
        keep_0 = list(range(g0, g1 + 1))
        if keep_0:
            new_blocks.append(dlogp_concat.iloc[:, keep_0])
            new_n_vec.append(len(keep_0))
            new_names.append(name)
            kept_total += len(keep_0)
        cursor += nk

    if not new_blocks:
        raise ValueError("Per-sector slice removed all sectors/columns.")
    dlogp_new = pd.concat(new_blocks, axis=1)
    print(f"[Info] Per-sector col slice [{per_start},{per_end}] -> kept {len(new_n_vec)} sectors / {kept_total} stocks.")
    return dlogp_new, new_n_vec, new_names, (per_start, per_end)

def sector_means_from_blocks(dlogp_concat, n_vec, sector_names):
    means = {}
    cursor = 0
    for name, nk in zip(sector_names, n_vec):
        cols = dlogp_concat.columns[cursor:cursor + nk]
        cursor += nk
        means[name] = dlogp_concat[cols].mean(axis=1, skipna=True)
    return pd.DataFrame(means, index=dlogp_concat.index)

# ----------------------------- Main -------------------------------------------
def main(argv):
    par = parse_kv_args(argv)
    os.makedirs(par["out_dir"], exist_ok=True)

    # 1) Load & build blocks
    sector_names, n_vec, date_list, dlogp_list = load_folder_as_blocks(par["data_dir"])
    common_dates = intersect_dates(date_list)
    aligned_list = align_blocks_to_common_dates(dlogp_list, common_dates)
    dlogp_concat = build_concatenated_matrix(aligned_list)

    # 2) Normalize by sqrt(Δ)
    m_all = dlogp_concat.shape[0]
    Delta = par["T"] / float(m_all)
    sqrt_Delta = math.sqrt(Delta)
    dlogp_concat = dlogp_concat / sqrt_Delta
    print(f"[Info] Normalized by sqrt(Δ): T={par['T']}, m_all={m_all}, Δ={Delta:.6g}")

    # 3) Row restriction
    dlogp_concat, row_used = restrict_rows_delta_index(dlogp_concat, par["row_start"], par["row_end"])

    # 4) Column restriction
    if par["per_sector"] == 1:
        dlogp_concat, n_vec_new, names_new, col_used = restrict_cols_per_sector(
            dlogp_concat, n_vec, sector_names, par["col_start"], par["col_end"])
        col_tag = f"each_{col_used[0]}-{col_used[1]}"
    else:
        dlogp_concat, n_vec_new, names_new, col_used = restrict_cols_global_blocks(
            dlogp_concat, n_vec, sector_names, par["col_start"], par["col_end"])
        col_tag = f"{col_used[0]}-{col_used[1]}"

    # 5) Sector means
    Z = sector_means_from_blocks(dlogp_concat, n_vec_new, names_new)
    m_used = Z.shape[0]
    t = pd.Series(np.arange(1, m_used + 1), index=Z.index, name="t")
    out_df = pd.concat([t, Z], axis=1)
    out_df.insert(1, "date", out_df.index)

    # Save CSV
    row_tag = f"{row_used[0]}-{row_used[1]}"
    csv_path = os.path.join(par["out_dir"], f"sector_means_rows_{row_tag}_cols_{col_tag}.csv")
    out_df.reset_index(drop=True).to_csv(csv_path, index=False)
    print(f"[Write] Sector means saved: {csv_path}")

    # Plot
    plt.figure(figsize=(10, 6))
    for name in Z.columns:
        plt.plot(out_df["t"], Z[name].values, label=name, linewidth=1.0, alpha=0.9)
    plt.axhline(0.0, linestyle="--", linewidth=0.8)
    plt.title(f"Sector mean ΔlogP / sqrt(Δ) (rows {row_tag}, cols {col_tag})")
    plt.xlabel("t (ΔlogP index)")
    plt.ylabel("mean ΔlogP / sqrt(Δ)")
    plt.legend(ncol=2, fontsize=9)
    plt.tight_layout()
    png_path = os.path.join(par["out_dir"], f"sector_means_rows_{row_tag}_cols_{col_tag}.png")
    plt.savefig(png_path, dpi=150)
    plt.close()
    print(f"[Plot] Figure saved: {png_path}")
    print("[Done] Sector mean computation complete.")

if __name__ == "__main__":
    main(sys.argv[1:])
