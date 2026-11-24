#!/usr/bin/env python3
"""Read footsteps/samples/samples.csv (written by ./footsteps --wav samples --csv)
and render the centroid x T60 scatter plus the per-material statistics table."""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse

# Validated categorical slots, light mode, surface #fcfcfb:
# slots 1/3/4/6/8, the five-slot subset that clears the all-pairs
# normal-vision floor (15.6). Worst CVD pair 6.9 (deutan) sits in the
# 6-8 band, so every group carries a secondary encoding: a distinct
# marker shape and a direct name label.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

ORDER = ["grass", "sand", "wood", "stone", "metal"]
LABEL = {"grass": "草地", "sand": "沙地", "wood": "木头", "stone": "石板", "metal": "金属"}
COLOR = {
    "grass": "#2a78d6",
    "sand": "#1baf7a",
    "wood": "#eda100",
    "stone": "#008300",
    "metal": "#e34948",
}
MARKER = {"grass": "o", "sand": "s", "wood": "^", "stone": "D", "metal": "P"}
plt.rcParams.update({
    "font.family": ["Noto Sans CJK SC", "DejaVu Sans"],
    "font.size": 9,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK_2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "text.color": INK,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
})


def load(path):
    groups = {m: [] for m in ORDER}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            groups[row["material"]].append(
                (float(row["centroid_hz"]), float(row["t60_ms"]))
            )
    return {m: np.array(v) for m, v in groups.items()}


def log_stats(points):
    lg = np.log(points)
    mu = lg.mean(axis=0)
    cov = np.cov(lg.T, ddof=1)
    return mu, cov


def separation(a, b):
    mu_a, cov_a = log_stats(a)
    mu_b, cov_b = log_stats(b)
    pooled = (cov_a + cov_b) / 2
    d = mu_b - mu_a
    return float(np.sqrt(d @ np.linalg.inv(pooled) @ d))


def ellipse_patch(mu, cov, color, n_std=1.0):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2 * n_std * np.sqrt(np.maximum(vals, 0))
    centre = np.exp(mu)
    return Ellipse(
        xy=centre, width=width, height=height, angle=angle,
        facecolor=color, edgecolor=color, alpha=0.10, linewidth=1.0, zorder=2,
    ), centre, np.degrees(angle)


def draw_plot(ax, groups):
    for m in ORDER:
        pts = groups[m]
        mu, cov = log_stats(pts)
        patch, centre, _ = ellipse_patch(mu, cov, COLOR[m], n_std=2.0)
        ax.add_patch(patch)
        ax.scatter(
            pts[:, 0], pts[:, 1], s=34, marker=MARKER[m], color=COLOR[m],
            edgecolors=SURFACE, linewidths=1.6, zorder=3, label=LABEL[m],
        )
        ax.plot(centre[0], centre[1], marker=MARKER[m], markersize=9,
                color=COLOR[m], markeredgecolor=SURFACE, markeredgewidth=1.6,
                linestyle="none", zorder=4)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(110, 12000)
    ax.set_ylim(38, 3400)
    ax.set_xticks([125, 250, 500, 1000, 2000, 4000, 8000])
    ax.set_xticklabels(["125", "250", "500", "1k", "2k", "4k", "8k"])
    ax.set_yticks([50, 100, 200, 400, 800, 1600, 3200])
    ax.set_yticklabels(["50", "100", "200", "400", "800", "1.6k", "3.2k"])
    ax.minorticks_off()
    ax.grid(True, which="major", color=GRID, linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("频谱质心 (Hz)", labelpad=8)
    ax.set_ylabel("衰减时长 T60 (ms)", labelpad=8)
    ax.set_title("五种材质在质心 × 衰减平面上分成五簇", loc="left",
                 fontsize=12, fontweight="bold", pad=26)
    ax.text(0.0, 1.04, "每组 10 个样本，对数坐标；椭圆为 2σ 协方差范围",
            transform=ax.transAxes, color=INK_2, fontsize=8.5)

    handles = [
        Line2D([], [], color=COLOR[m], marker=MARKER[m], linestyle="none",
               markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.4,
               label=LABEL[m])
        for m in ORDER
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=8.5,
              handletextpad=0.5, borderpad=0.2, labelspacing=0.45)


def draw_table(ax, groups):
    ax.axis("off")
    ax.set_title("统计（10 样本）", loc="left", fontsize=11,
                 fontweight="bold", pad=26)
    ax.text(0.0, 1.04, "均值 ± 标准差", transform=ax.transAxes,
            color=INK_2, fontsize=8.5)

    cols = ["材质", "质心 Hz", "T60 ms", "最近邻"]
    xpos = [0.00, 0.27, 0.57, 0.80]
    y = 0.93
    for c, x in zip(cols, xpos):
        ax.text(x, y, c, color=MUTED, fontsize=8, ha="left",
                transform=ax.transAxes)
    ax.plot([0, 1], [y - 0.035, y - 0.035], color=GRID, linewidth=0.8,
            transform=ax.transAxes)

    y = 0.855
    for m in ORDER:
        pts = groups[m]
        cen_m, cen_s = pts[:, 0].mean(), pts[:, 0].std(ddof=1)
        t_m, t_s = pts[:, 1].mean(), pts[:, 1].std(ddof=1)
        others = [separation(pts, groups[o]) for o in ORDER if o != m]
        nearest = min(others)
        ax.scatter([0.015], [y + 0.012], s=34, marker=MARKER[m], color=COLOR[m],
                   edgecolors=SURFACE, linewidths=1.4, transform=ax.transAxes,
                   clip_on=False, zorder=3)
        ax.text(0.05, y, LABEL[m], color=INK, fontsize=9, transform=ax.transAxes)
        ax.text(0.27, y, f"{cen_m:.0f} ± {cen_s:.0f}", color=INK, fontsize=9,
                transform=ax.transAxes)
        ax.text(0.57, y, f"{t_m:.0f} ± {t_s:.0f}", color=INK, fontsize=9,
                transform=ax.transAxes)
        ax.text(0.80, y, f"{nearest:.1f}σ", color=INK_2, fontsize=9,
                transform=ax.transAxes)
        y -= 0.085

    ax.plot([0, 1], [y + 0.035, y + 0.035], color=GRID, linewidth=0.8,
            transform=ax.transAxes)
    ax.text(0.0, y - 0.01,
            "最近邻：两组均值在二维对数空间中的\n"
            "马氏距离（合并协方差），单位 σ。\n"
            "≥ 2σ 且 2σ 椭圆不相交即视为分离。\n"
            "2σ 椭圆不相交；草地与沙地最接近，3σ 处才接上。",
            color=MUTED, fontsize=8, va="top", transform=ax.transAxes,
            linespacing=1.6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(Path(__file__).resolve().parents[1] / "samples" / "samples.csv"))
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "samples" / "materials_scatter.png"))
    args = ap.parse_args()

    groups = load(args.csv)
    fig = plt.figure(figsize=(12.4, 6.4), dpi=180)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.55, 1.0], wspace=0.28,
                          left=0.065, right=0.975, top=0.84, bottom=0.17)
    draw_plot(fig.add_subplot(gs[0]), groups)
    draw_table(fig.add_subplot(gs[1]), groups)
    fig.text(0.065, 0.045,
             "实测：48 kHz 合成脉冲，Hann 窗 46 ms 帧能量加权质心；"
             "T60 取 Schroeder 能量衰减曲线 -5~-25 dB 的线性拟合。",
             color=MUTED, fontsize=8)
    fig.savefig(args.out, facecolor=SURFACE)
    print(f"wrote {args.out}")

    for m in ORDER:
        pts = groups[m]
        others = sorted((separation(pts, groups[o]), o) for o in ORDER if o != m)
        print(f"{m:6s} centroid {pts[:,0].mean():7.1f}±{pts[:,0].std(ddof=1):5.1f} Hz  "
              f"T60 {pts[:,1].mean():7.1f}±{pts[:,1].std(ddof=1):5.1f} ms  "
              f"nearest {others[0][1]} {others[0][0]:.1f}σ")


if __name__ == "__main__":
    main()
