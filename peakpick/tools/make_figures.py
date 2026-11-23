#!/usr/bin/env python3
"""peakpick 结果出图（读 results/ 下 CSV，写 figures/ 下 PNG+SVG）。

图 1 feasible_region  4 窗 x (L, df) 三态可行域小多联热力图
图 2 lcrit            临界窗长 Lcrit(df)：实测台阶 + 理论双曲线，分 res/full 两面板
图 3 movement_rate    移动率 |dLcrit/ddf|：拟合曲线 + 原始整数差分
图 4 fp_fn_counts     全网格误检/漏检聚合计数（分组柱）
图 5 counterexample   可行域外实测反例：rect 旁瓣误检 + blackman 合并漏检
数据表（table view）：results/*.csv。
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/ys/.local/lib/python3.12/site-packages")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.ticker import ScalarFormatter
from matplotlib import font_manager

_FONT_TTC = "/home/ys/.local/share/fonts/wqy-zenhei.ttc"
if os.path.exists(_FONT_TTC):
    font_manager.fontManager.addfont(_FONT_TTC)

from common import (FS, WIN_NAMES, window, two_tone, mag_spectrum,
                    pick_peaks, GATE_DB)

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")

# ---- dataviz 调色板（validate_palette.js 已验证）----
SERIES = {"rect": "#2a78d6", "hann": "#eb6834",
          "hamming": "#1baf7a", "blackman": "#eda100"}
SURFACE, GRID, INK, MUTED = "#fcfcfb", "#e1e0d9", "#0b0b0b", "#52514e"
# 三态：可行=蓝阶深色；漏检主导=中性灰；误检主导=红色系
C_OK, C_FN, C_FP = "#2a78d6", "#c3c2b7", "#e34948"

plt.rcParams.update({
    "font.family": "WenQuanYi Zen Hei",
    "font.size": 10,
    "text.color": INK, "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.axisbelow": True,
    "axes.unicode_minus": False,
})


def load(name):
    with open(os.path.join(RES, name)) as f:
        return list(csv.DictReader(f))


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"),
                    dpi=160, bbox_inches="tight")
    plt.close(fig)


def log_scalar(ax):
    """对数轴用十进制标量标签（避免 unicode_minus=False 下的 10^x-1 坏标签）。"""
    fmt = ScalarFormatter()
    fmt.set_scientific(False)
    ax.yaxis.set_major_formatter(fmt)
    ax.yaxis.set_minor_formatter(plt.NullFormatter())


# ---------- 图 1：可行域 ----------
def fig_feasible(rows):
    Ls = sorted({int(r["L"]) for r in rows})
    dfs = sorted({int(r["df"]) for r in rows})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), sharex=True, sharey=True)
    cmap = ListedColormap([C_FN, C_FP, C_OK])
    norm = BoundaryNorm([0, 1, 2, 3], cmap.N)
    for ax, win in zip(axes.flat, WIN_NAMES):
        z = np.zeros((len(dfs), len(Ls)))
        for r in rows:
            if r["window"] != win:
                continue
            i, j = dfs.index(int(r["df"])), Ls.index(int(r["L"]))
            if r["clean"] == "True":
                z[i, j] = 2
            elif int(r["fp_sum"]) > 0 and int(r["fn_sum"]) == 0:
                z[i, j] = 1          # 只误检（旁瓣问题）
            else:
                z[i, j] = 0          # 含漏检（含误+漏）
        ax.pcolormesh(np.arange(len(Ls) + 1), np.arange(len(dfs) + 1), z,
                      cmap=cmap, norm=norm, shading="flat")
        ax.set_title(win, color=SERIES[win], fontweight="bold")
        ax.set_xticks(np.arange(len(Ls)) + 0.5)
        ax.set_xticklabels(Ls, rotation=45, fontsize=8)
        ax.set_yticks(np.arange(0, len(dfs), 5) + 0.5)
        ax.set_yticklabels(dfs[::5], fontsize=8)
        ax.grid(False)
    for ax in axes[1]:
        ax.set_xlabel("窗长 L（样本，2 的幂）")
    for ax in axes[:, 0]:
        ax.set_ylabel("双音间距 Δf（Hz）")
    legend = [Patch(facecolor=C_OK, label="可行：零误检零漏检"),
              Patch(facecolor=C_FP, label="不可行：旁瓣误检（FN=0）"),
              Patch(facecolor=C_FN, label="不可行：主峰合并漏检")]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("窗型 × 窗长 × 间距 组合可行域（每单元 16 个确定性试验，-40 dB 报告门限）",
                 fontsize=12, y=0.98)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    save(fig, "fig1_feasible_region")


# ---------- 图 2：临界窗长 ----------
def fig_lcrit(lc):
    dfs = sorted({int(r["df"]) for r in lc})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for win in WIN_NAMES:
        sub = [r for r in lc if r["window"] == win]
        d = np.array([int(r["df"]) for r in sub])
        yr = np.array([float(r["lcrit_res"]) for r in sub])
        cth = float(sub[0]["crit_bin_theory"])
        ax1.step(d, yr, where="mid", color=SERIES[win], lw=1.8, label=win)
        ax1.plot(d, cth * FS / d, "--", color=SERIES[win], lw=1.0, alpha=0.7)
        xf = [int(r["df"]) for r in sub if r["lcrit_full"]]
        yf = [float(r["lcrit_full"]) for r in sub if r["lcrit_full"]]
        if xf:
            ax2.step(xf, yf, where="mid", color=SERIES[win], lw=1.8,
                     label=f"{win}（Δf≥{min(xf)}Hz 可行）")
        else:
            ax2.plot([], [], color=SERIES[win], lw=1.8,
                     label=f"{win}（全域不可行）")
    ax1.set_title("分辨临界 Lcrit_res：零漏检（误检不约束）")
    ax2.set_title("严格临界 Lcrit_full：零误检且零漏检")
    for ax in (ax1, ax2):
        ax.set_xlabel("双音间距 Δf（Hz）")
        ax.set_yscale("log")
        log_scalar(ax)
    ax1.set_ylabel("临界窗长 Lcrit（样本，对数轴）")
    # 理论线说明并入左图图例
    handles, labels = ax1.get_legend_handles_labels()
    handles.append(plt.Line2D([], [], ls="--", color=MUTED, lw=1.0))
    labels.append("理论 L=c·fs/Δf")
    ax1.legend(handles, labels, frameon=False, fontsize=9)
    ax2.legend(frameon=False, fontsize=9)
    fig.suptitle("临界窗长随间距移动：实测台阶 vs 理论双曲线 L = c·fs/Δf",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, "fig2_lcrit")


# ---------- 图 3：移动率 ----------
def fig_rate(lc):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for win in WIN_NAMES:
        sub = sorted((r for r in lc if r["window"] == win),
                     key=lambda r: int(r["df"]))
        d = np.array([int(r["df"]) for r in sub])
        rf = np.array([float(r["rate_fit_res"]) for r in sub])
        rr = np.array([float(r["rate_raw_res"]) for r in sub])
        ax1.plot(d, -rf, color=SERIES[win], lw=1.8, label=win)
        ax1.plot(d, np.abs(rr), "o", color=SERIES[win], ms=3, alpha=0.45)
        full = [r for r in sub if r["rate_fit_full"]]
        if full:
            df_ = np.array([int(r["df"]) for r in full])
            ff = np.array([float(r["rate_fit_full"]) for r in full])
            rr2 = np.array([float(r["rate_raw_full"]) for r in full])
            ax2.plot(df_, -ff, color=SERIES[win], lw=1.8, label=win)
            ax2.plot(df_, np.abs(rr2), "o", color=SERIES[win], ms=3, alpha=0.45)
        else:
            ax2.plot([], [], color=SERIES[win], lw=1.8,
                     label=f"{win}（全域不可行）")
    ax1.set_title("分辨临界移动率 |dLcrit_res/dΔf|")
    ax2.set_title("严格临界移动率 |dLcrit_full/dΔf|")
    for ax in (ax1, ax2):
        ax.set_xlabel("双音间距 Δf（Hz）")
        ax.set_yscale("log")
        log_scalar(ax)
        ax.legend(frameon=False, fontsize=9)
    ax1.set_ylabel("移动率（样本/Hz，绝对值，对数轴）")
    handles, labels = ax1.get_legend_handles_labels()
    handles.append(plt.Line2D([], [], marker="o", color=MUTED, lw=0,
                              ms=4, alpha=0.6))
    labels.append("整数窗长原始差分")
    ax1.legend(handles, labels, frameon=False, fontsize=9)
    fig.suptitle("移动率曲线：|dLcrit/ddf| = c_hat·fs/df²（线=尺度律拟合，点=原始差分）",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, "fig3_movement_rate")


# ---------- 图 4：误检漏检计数 ----------
def fig_counts(rows):
    fp = {w: 0 for w in WIN_NAMES}
    fn = {w: 0 for w in WIN_NAMES}
    cells_err = {w: 0 for w in WIN_NAMES}
    for r in rows:
        w = r["window"]
        fp[w] += int(r["fp_sum"])
        fn[w] += int(r["fn_sum"])
        if r["clean"] != "True":
            cells_err[w] += 1
    x = np.arange(len(WIN_NAMES))
    fig, ax = plt.subplots(figsize=(9, 5.2))
    b1 = ax.bar(x - 0.2, [fp[w] for w in WIN_NAMES], width=0.38,
                color=C_FP, label="误检 FP 合计")
    b2 = ax.bar(x + 0.2, [fn[w] for w in WIN_NAMES], width=0.38,
                color="#6250d6", label="漏检 FN 合计")
    for bars in (b1, b2):
        for b in bars:
            ax.annotate(f"{int(b.get_height())}",
                        (b.get_x() + b.get_width() / 2, b.get_height()),
                        ha="center", va="bottom", fontsize=9, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels(WIN_NAMES)
    ax.set_ylabel(f"全网格错误计数（{len(rows)//4} 单元 × 16 试验）")
    ax.set_title("误检 / 漏检计数：旁瓣型窗（rect/hann）误检主导，"
                 "宽主瓣窗小窗长下漏检主导")
    ax.legend(frameon=False)
    save(fig, "fig4_fp_fn_counts")
    return fp, fn, cells_err


# ---------- 图 5：可行域外实测反例 ----------
def fig_counterexample():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # 反例 A：rect L=512 单音，-40dB 门内 14 个旁瓣假峰（可行域外：rect 全域不可行）
    L = 512
    n = np.arange(L)
    bin_hz = FS / L
    f0 = 12 * bin_hz
    x = np.cos(2 * np.pi * f0 * n / FS)
    fr, mag = mag_spectrum(x, window("rect", L))
    magdb = 20 * np.log10(mag / mag.max() + 1e-12)
    sel = (fr >= f0 - 8 * bin_hz) & (fr <= f0 + 8 * bin_hz)
    ax1.plot(fr[sel], magdb[sel], color=SERIES["rect"], lw=1.4)
    pk = pick_peaks(fr, mag)
    pk = pk[np.abs(pk - f0) <= 8 * bin_hz]
    fp = pk[np.abs(pk - f0) > 0.5 * bin_hz]
    ax1.axhline(-GATE_DB, color=MUTED, ls="--", lw=1.1,
                label=f"-{GATE_DB:.0f} dB 报告门限")
    ax1.plot(fp, [-GATE_DB - 2] * len(fp), "x", color=C_FP, ms=8,
             label=f"旁瓣误检 {len(fp)} 个（无真实音）")
    ax1.plot(f0, 0, "o", color="#008300", ms=7, label="真实音 1 个")
    ax1.set_title(f"反例 A（可行域外）：rect 窗 L={L} 单音，"
                  f"最高旁瓣 -13.3 dB > 门限，产生 {len(fp)} 个误检")
    ax1.set_xlabel("频率（Hz）")
    ax1.set_ylabel("幅值（dB，归一化）")
    ax1.set_ylim(-70, 5)
    ax1.legend(frameon=False, fontsize=9, loc="upper right")

    # 反例 B：blackman L=256 双音 df=300Hz（1.60 bin），合并成 1 个峰 -> 漏检 1
    L2, df = 256, 300
    t0, t1, x2 = two_tone(L2, df, 0.0, 0.0)
    fr2, mag2 = mag_spectrum(x2, window("blackman", L2))
    m2db = 20 * np.log10(mag2 / mag2.max() + 1e-12)
    bhz = FS / L2
    sel2 = (fr2 >= t0 - 4 * bhz) & (fr2 <= t1 + 4 * bhz)
    ax2.plot(fr2[sel2], m2db[sel2], color=SERIES["blackman"], lw=1.4)
    pk2 = pick_peaks(fr2, mag2)
    pk2 = pk2[(pk2 >= t0 - 4 * bhz) & (pk2 <= t1 + 4 * bhz)]
    ax2.axhline(-GATE_DB, color=MUTED, ls="--", lw=1.1)
    ax2.plot(pk2, [2] * len(pk2), "v", color="#6250d6", ms=9,
             label=f"报告峰仅 {len(pk2)} 个")
    ax2.axvline(t0, color="#008300", lw=1.2, alpha=0.8)
    ax2.axvline(t1, color="#008300", lw=1.2, alpha=0.8,
                label=f"真实音 2 个（间距 {df}Hz = {df / bhz:.2f} bin）")
    ax2.annotate("漏检 1：两主峰叠加成单峰",
                 xy=(pk2[0], 0), xytext=(pk2[0] + 260, -12),
                 arrowprops=dict(arrowstyle="->", color=INK), fontsize=9)
    ax2.set_title(f"反例 B（可行域外）：blackman 窗 L={L2}，间距 {df}Hz < "
                  f"Lcrit 对应 {2.14 * bhz:.0f}Hz，两峰合并")
    ax2.set_xlabel("频率（Hz）")
    ax2.set_ylabel("幅值（dB，归一化）")
    ax2.set_ylim(-70, 8)
    ax2.legend(frameon=False, fontsize=9, loc="upper right")

    fig.suptitle("可行域外实测反例（均由分析引擎确定性复现，C 自测同参数）",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, "fig5_counterexample")


def main():
    rows = load("feasible.csv")
    lc = load("lcrit.csv")
    fig_feasible(rows)
    fig_lcrit(lc)
    fig_rate(lc)
    fp, fn, cells_err = fig_counts(rows)
    fig_counterexample()
    print("FP totals:", fp)
    print("FN totals:", fn)
    print("error cells:", cells_err)


if __name__ == "__main__":
    main()
