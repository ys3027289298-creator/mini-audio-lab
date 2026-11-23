#!/usr/bin/env python3
"""全网格扫描 + 临界窗长/移动率求解。

输出（results/ 下）：
  grid_fp.csv / grid_fn.csv      4 窗 x 9 窗长 x 37 间距 的误检/漏检计数
  feasible.csv                   每个单元的 clean 标志与计数（热力图源）
  lcrit.csv                      分辨临界/严格临界窗长、理论值、移动率
  gate_sensitivity.csv           单音旁瓣假峰数随门限变化
  summary.json                   汇总指标与校验
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (FS, WIN_NAMES, cell_totals, gate_sensitivity,
                    measured_lcrit)

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")

L_GRID = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
DF_GRID = list(range(40, 401, 10))                       # Hz
# 理论临界间距（bin，无限过采样拐点法，theory_study.py 结果）
CRIT_BIN_THEORY = {"rect": 1.324, "hann": 1.885,
                   "hamming": 1.670, "blackman": 2.068}


def main():
    os.makedirs(RES, exist_ok=True)
    rows = []
    for win in WIN_NAMES:
        for L in L_GRID:
            for df in DF_GRID:
                t = cell_totals(win, L, df)
                rows.append({"window": win, "L": L, "df": df, **t})
                print(f"{win:8s} L={L:5d} df={df:4d} "
                      f"FP={t['fp_sum']:2d} FN={t['fn_sum']:2d}", flush=True)

    with open(os.path.join(RES, "grid_fp.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # 可行性谓词单调性校验：
    #  - clean（零误检且零漏检）随 L 必须全局 False*->True*；
    #  - zero_fn 只在 d_bin>=1 的物理可分上分支单调；d_bin<1 时宽主瓣窗
    #    可能在亚 bin 间距仍撑出两个局部极大（见 README），不构成回退判据。
    import math
    violations = []
    for win in WIN_NAMES:
        for df in DF_GRID:
            sub = sorted((r for r in rows if r["window"] == win and r["df"] == df),
                         key=lambda r: r["L"])
            cflags = [r["clean"] for r in sub]
            if any(cflags[i] and not cflags[i + 1]
                   for i in range(len(cflags) - 1)):
                violations.append((win, df, "clean"))
            up = [(r["L"], r["fn_sum"] == 0) for r in sub
                  if df * r["L"] / FS >= 1.0]
            if any(up[i][1] and not up[i + 1][1]
                   for i in range(len(up) - 1)):
                violations.append((win, df, "zero_fn_upper_branch"))
    assert not violations, f"可行性谓词非单调: {violations[:5]}"

    # 临界窗长：网格夹住后在相邻 2 的幂间二分
    lcrit_rows = []
    for win in WIN_NAMES:
        c_th = CRIT_BIN_THEORY[win]
        for df in DF_GRID:
            sub = sorted((r for r in rows if r["window"] == win and r["df"] == df),
                         key=lambda r: r["L"])
            res_grid = [r["L"] for r in sub if r["fn_sum"] == 0]
            full_grid = [r["L"] for r in sub if r["clean"]]

            def refine(mode, grid_list):
                if not grid_list:
                    return None
                Lg = min(grid_list)
                if mode == "res":
                    # 限制在 d_bin>=1 上分支，避开亚 bin 双极大区
                    lo = max(L_GRID[max(0, L_GRID.index(Lg) - 1)],
                             math.ceil(FS / df))
                else:
                    lo = L_GRID[max(0, L_GRID.index(Lg) - 1)]
                return measured_lcrit(win, df, mode=mode,
                                      L_lo=min(lo, Lg), L_hi=Lg)

            L_res = refine("res", res_grid)
            L_full = refine("full", full_grid)
            lcrit_rows.append({
                "window": win, "df": df,
                "lcrit_res": L_res if L_res is not None else "",
                "lcrit_full": L_full if L_full is not None else "",
                "lcrit_theory": c_th * FS / df,
                "crit_bin_theory": c_th,
            })
            print(f"lcrit {win:8s} df={df:4d} res={L_res} full={L_full} "
                  f"theory={c_th * FS / df:8.1f}", flush=True)

    # 移动率：Lc 对 fs/df 过原点拟合 c_hat（窗型常数，单位 bin），
    # dLc/ddf = -c_hat*fs/df^2；另给整数窗长原始差分对照。
    for mode_col in ("lcrit_res", "lcrit_full"):
        for win in WIN_NAMES:
            sub = [r for r in lcrit_rows
                   if r["window"] == win and r[mode_col] != ""]
            if not sub:
                continue
            x = np.array([FS / r["df"] for r in sub])
            y = np.array([r[mode_col] for r in sub], dtype=float)
            c_hat = float(np.dot(x, y) / np.dot(x, x))
            d = np.array([r["df"] for r in sub], dtype=float)
            order = np.argsort(d)
            raw = np.gradient(y[order], d[order])
            tag = mode_col.split("_")[1]
            for r, rr in zip((sub[i] for i in order), raw):
                r[f"rate_fit_{tag}"] = -c_hat * FS / r["df"] ** 2
                r[f"rate_raw_{tag}"] = rr
                r[f"c_hat_{tag}"] = c_hat
    for r in lcrit_rows:
        for k in ("rate_fit_res", "rate_raw_res", "c_hat_res",
                  "rate_fit_full", "rate_raw_full", "c_hat_full"):
            r.setdefault(k, "")

    with open(os.path.join(RES, "lcrit.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lcrit_rows[0].keys()))
        w.writeheader()
        w.writerows(lcrit_rows)

    with open(os.path.join(RES, "feasible.csv"), "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["window", "L", "df", "clean", "fp_sum", "fn_sum"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in
                        ("window", "L", "df", "clean", "fp_sum", "fn_sum")})

    # 门限敏感度
    gs = gate_sensitivity()
    with open(os.path.join(RES, "gate_sensitivity.csv"), "w", newline="") as f:
        gates = sorted(next(iter(gs.values())).keys())
        w = csv.writer(f)
        w.writerow(["window"] + [f"fp_gate_{abs(g)}dB" for g in gates])
        for win in WIN_NAMES:
            w.writerow([win] + [gs[win][g] for g in gates])

    summary = {
        "fs": FS, "gate_db": 40.0, "zp": 4,
        "n_trials_per_cell": 16,
        "L_grid": L_GRID,
        "df_grid_hz": [DF_GRID[0], DF_GRID[-1], DF_GRID[1] - DF_GRID[0]],
        "crit_bin_theory": CRIT_BIN_THEORY,
        "total_cells": len(rows),
        "cells_with_any_error": sum(1 for r in rows if not r["clean"]),
        "total_fp": sum(r["fp_sum"] for r in rows),
        "total_fn": sum(r["fn_sum"] for r in rows),
        "feasibility_predicate_monotone": True,
        "gate_sensitivity_fp": gs,
        "c_hat_res_bin": {
            w: next((r["c_hat_res"] for r in lcrit_rows
                     if r["window"] == w and r["c_hat_res"] != ""), None)
            for w in WIN_NAMES},
        "c_hat_full_bin": {
            w: next((r["c_hat_full"] for r in lcrit_rows
                     if r["window"] == w and r["c_hat_full"] != ""), None)
            for w in WIN_NAMES},
        "full_feasible_windows": sorted({
            r["window"] for r in lcrit_rows if r["lcrit_full"] != ""}),
    }
    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
