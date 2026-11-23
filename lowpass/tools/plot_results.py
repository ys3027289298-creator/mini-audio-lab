#!/usr/bin/env python3
"""Plot the Pareto front and the resampler-cascade results.

Reads sweep.csv / cascade.csv produced by the lowpass binary and writes
pareto_front.png and cascade.png into this directory.
Run from lowpass/: python3 tools/plot_results.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")

# dataviz reference palette (validated separately)
BLUE = "#2a78d6"      # non-dominated / Pareto front
RED = "#e34948"       # the three selected non-dominated points
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"
GRAY = "#898781"      # dominated candidates (ink, not a data hue)

for path in ("/home/ys/.local/share/fonts/wqy-zenhei.ttc",):
    if os.path.exists(path):
        font_manager.fontManager.addfont(path)
plt.rcParams.update({
    "font.family": "WenQuanYi Zen Hei",
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.unicode_minus": False,
    "font.size": 11,
})


def load_sweep():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "sweep.csv"))))
    fc = [float(r["fc_frac"]) for r in rows]
    a = [float(r["alias_residual_db"]) for r in rows]
    p = [float(r["phase_distortion_rad"]) for r in rows]
    front = [int(r["pareto"]) == 1 for r in rows]
    return fc, a, p, front


def plot_pareto():
    fc, a, p, front = load_sweep()
    fa = [v for v, ok in zip(a, front) if ok]
    fp = [v for v, ok in zip(p, front) if ok]
    order = sorted(range(len(fa)), key=lambda i: fa[i])
    fa = [fa[i] for i in order]
    fp = [fp[i] for i in order]

    fig, ax = plt.subplots(figsize=(9, 6.2), dpi=150)

    # front connecting line first (under the markers), 2px
    ax.plot(fa, fp, color=BLUE, linewidth=2, zorder=2,
            solid_capstyle="round", label="帕累托前沿")

    # dominated candidates: muted gray, smaller
    da = [v for v, ok in zip(a, front) if not ok]
    dp = [v for v, ok in zip(p, front) if not ok]
    ax.scatter(da, dp, s=55, color=GRAY, zorder=3, marker="o",
               edgecolor=SURFACE, linewidth=1.5, label="被支配候选")

    # non-dominated front points: blue, >=8px, surface ring
    ax.scatter(fa, fp, s=95, color=BLUE, zorder=4, marker="o",
               edgecolor=SURFACE, linewidth=1.5, label="非支配候选")

    # three selected points: red, larger, on top
    picks = [order[0], order[len(order) // 2], order[-1]]
    pa = [fa[i] for i in picks]
    pp = [fp[i] for i in picks]
    names = ["P1 混叠最小", "P2 折中", "P3 相位最小"]
    ax.scatter(pa, pp, s=170, color=RED, zorder=5, marker="o",
               edgecolor=SURFACE, linewidth=1.5, label="选取的三个非支配点")

    # selective direct labels for the three picks only, with offsets to
    # avoid covering marks
    offsets = [(14, 12), (-58, 22), (-8, 16)]
    for x, y, name, (dx, dy) in zip(pa, pp, names, offsets):
        ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                    color=INK, fontsize=10.5,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))

    ax.set_xlabel("目标 A：混叠残留 (dB，越小越好)")
    ax.set_ylabel("目标 P：相位失真 (rad，越小越好)")
    ax.set_title("低通截止频率扫描的帕累托前沿\n(4 阶 Butterworth，fs=1000 Hz，15 个候选点)",
                 fontsize=13, pad=12)
    ax.legend(frameon=False, loc="upper right", fontsize=10)
    ax.margins(0.10)
    fig.tight_layout()
    out = os.path.join(ROOT, "pareto_front.png")
    fig.savefig(out)
    plt.close(fig)
    return out, list(zip(pa, pp))


def plot_cascade():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "cascade.csv"))))
    labels = [r["label"] for r in rows]
    fold = [float(r["cascade_fold_vs_pass_db"]) for r in rows]
    gain = [float(r["passband_gain_db"]) for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.2), dpi=150)
    colors = [GRAY, RED, RED, RED]

    bars1 = ax1.bar(labels, fold, color=colors, width=0.62, zorder=3,
                    edgecolor=SURFACE, linewidth=1.5)
    ax1.axhline(0, color=AXIS, linewidth=1.0, zorder=2)
    ax1.set_ylabel("折叠能量 / 通带能量 (dB)")
    ax1.set_title("串联 2:1 线性重采样后的混叠折叠", fontsize=12)
    ax1.tick_params(axis="x", labelrotation=20)
    lo1 = min(fold)
    ax1.set_ylim(lo1 - 6, max(fold) + 5)
    for b, v in zip(bars1, fold):
        ax1.annotate(f"{v:.1f}", (b.get_x() + b.get_width() / 2, v),
                     xytext=(0, 5 if v >= 0 else -13),
                     textcoords="offset points", ha="center",
                     color=INK2, fontsize=9.5)

    bars2 = ax2.bar(labels, gain, color=colors, width=0.62, zorder=3,
                    edgecolor=SURFACE, linewidth=1.5)
    ax2.axhline(0, color=AXIS, linewidth=1.0, zorder=2)
    ax2.set_ylabel("通带平均增益 (dB)")
    ax2.set_title("通带保真度（越接近 0 越好）", fontsize=12)
    ax2.tick_params(axis="x", labelrotation=20)
    lo2 = min(gain)
    ax2.set_ylim(lo2 - 5, 3)
    for b, v in zip(bars2, gain):
        ax2.annotate(f"{v:.2f}", (b.get_x() + b.get_width() / 2, v),
                     xytext=(0, 5 if v >= 0 else -13),
                     textcoords="offset points", ha="center",
                     color=INK2, fontsize=9.5)

    legend_items = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY,
               markersize=9, markeredgecolor=SURFACE, label="无滤波基线"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=RED,
               markersize=9, markeredgecolor=SURFACE, label="前沿选取点 P1/P2/P3"),
    ]
    fig.legend(handles=legend_items, frameon=False, loc="lower center",
               ncol=2, fontsize=10, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("低通 + 线性重采样串联实测", fontsize=13, y=0.99)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    out = os.path.join(ROOT, "cascade.png")
    fig.savefig(out)
    plt.close(fig)
    return out


if __name__ == "__main__":
    p1, picks = plot_pareto()
    p2 = plot_cascade()
    print("wrote", p1)
    print("wrote", p2)
    print("selected points (alias_dB, phase_rad):")
    for pt in picks:
        print("  %.4f, %.6f" % pt)
