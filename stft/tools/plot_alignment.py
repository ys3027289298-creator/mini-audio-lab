#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STFT framing alignment — analytic region / relaxation / phase / budget figures.

Inputs  : artifacts/{summary,phase_sweep,timestamps_pass,timestamps_fail}.csv
Outputs : artifacts/fig_{region,relax,phase,budget}.png
          artifacts/analytic_region.csv  (rho-solved N/H limits + chosen point)

Palette validated with the dataviz validator (light mode):
  blue #2a78d6 (series 1), orange #eb6834 (series 2), aqua #1baf7a (series 3),
  critical red #d03b3b (out-of-budget), surface #fcfcfb.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FS = 48000.0
DELTA_MS = 5.0
DELTA_S = DELTA_MS / 1000.0
A = FS * DELTA_S                      # 240 samples: (H+N)/2 <= A
RHO_MAX = 8.0

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RED, SURFACE, INK, MUTED, GRID = "#d03b3b", "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "font.size": 10, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "legend.frameon": False,
})


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def load_summary(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            for k in ("N", "H"):
                r[k] = int(r[k])
            for k in ("rho", "grid_ms", "smear_ms", "budget_ms",
                      "meas_nearest_ms", "meas_onset_ms", "delta_max_ms"):
                r[k] = float(r[k])
            rows.append(r)
    return rows


def fig_region(rows, out_csv, out_png):
    # analytic rho sweep -> N_max, H_max in samples and ms
    rho = np.linspace(1.0, RHO_MAX, 281)
    N_max = 2.0 * A * rho / (rho + 1.0)                 # samples
    H_max = 2.0 * A / (rho + 1.0)
    N_ms, H_ms = N_max / FS * 1000.0, H_max / FS * 1000.0

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rho", "N_max_samples", "H_max_samples",
                    "N_max_ms", "H_max_ms", "budget_ms"])
        for r, n, h, nm, hm in zip(rho, N_max, H_max, N_ms, H_ms):
            w.writerow([f"{r:.4f}", f"{n:.3f}", f"{h:.3f}",
                        f"{nm:.4f}", f"{hm:.4f}", f"{(n+h)/(2*FS)*1000:.4f}"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # left: feasible region in the (H, N) plane
    Hv = np.linspace(0.0, 2 * A, 300)
    Nline = 2 * A - Hv                                   # H + N <= 2A
    ax1.fill_between(Hv, 0, Nline, color=BLUE, alpha=0.14, label="feasible: H+N ≤ 480")
    ax1.plot(Hv, Nline, color=BLUE, lw=2, label="boundary H+N = 480 (δ* = 5 ms)")
    ax1.plot(Hv, Hv, color=MUTED, lw=1.4, ls="--", label="N = H (ρ = 1, no overlap)")
    # per-point label offsets tuned so the six annotations do not collide
    label_dx = {"chosen": (12, 12), "dense": (10, 34), "small": (12, -18),
                "critical": (12, -20), "wide": (-138, 36), "critbig": (-150, -18)}
    for r in rows:
        c = AQUA if r["verdict"] == "PASS" else RED
        ax1.scatter(r["H"], r["N"], s=64, color=c, zorder=5,
                    edgecolor=SURFACE, linewidth=1.5)
        dx, dy = label_dx[r["tag"].split()[0]]
        ax1.annotate(f'N={r["N"]}, H={r["H"]}, ρ={r["rho"]:g}',
                     (r["H"], r["N"]), textcoords="offset points",
                     xytext=(dx, dy), fontsize=8, color=MUTED)
    ax1.set_xlim(-40, 580); ax1.set_ylim(0, 600)
    ax1.set_xlabel("hop H (samples @48 kHz)")
    ax1.set_ylabel("window N (samples @48 kHz)")
    ax1.set_title("Feasible (H, N) region under δ* = (H+N)/(2·fs) ≤ 5 ms")
    ax1.legend(loc="lower right", fontsize=8)
    style(ax1)

    # right: rho-solved limits
    ax2.plot(rho, N_ms, color=ORANGE, lw=2, label="N_max(ρ) = 2·fs·Δ·ρ/(ρ+1)")
    ax2.plot(rho, H_ms, color=BLUE, lw=2, label="H_max(ρ) = 2·fs·Δ/(ρ+1)")
    ax2.axhline(DELTA_MS * 2, color=MUTED, lw=1, ls=":")
    ax2.annotate("N_max → 2·fs·Δ = 10 ms", (4.6, 9.6), fontsize=8, color=MUTED)
    # chosen operating point rho=2
    i2 = np.argmin(abs(rho - 2.0))
    ax2.scatter([2, 2], [N_ms[i2], H_ms[i2]], s=64,
                color=[ORANGE, BLUE], zorder=5, edgecolor=SURFACE, linewidth=1.5)
    ax2.annotate(f"chosen ρ=2: N≤{N_ms[i2]:.2f} ms, H≤{H_ms[i2]:.2f} ms",
                 (2.05, (N_ms[i2] + H_ms[i2]) / 2), fontsize=8, color=MUTED)
    ax2.set_xlabel("oversampling factor ρ = N/H")
    ax2.set_ylabel("limit (ms)")
    ax2.set_title("Parameter limits solved per oversampling factor")
    ax2.legend(loc="center right", fontsize=8)
    style(ax2)

    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def fig_relax(out_png):
    rho = np.linspace(1.0, RHO_MAX, 400)
    R = DELTA_MS * (rho - 1.0) / (rho + 1.0)             # slack in ms
    N_budget = DELTA_MS * 2.0 * rho / (rho + 1.0)       # window half-width budget (ms)

    fig, ax1 = plt.subplots(figsize=(8.6, 4.6))
    ax1.plot(rho, R, color=AQUA, lw=2, label="relaxation R(ρ) = Δ·(ρ−1)/(ρ+1)")
    ax1.plot(rho, N_budget, color=ORANGE, lw=2,
             label="window budget N/(2fs) = 2Δ·ρ/(ρ+1)")
    ax1.axhline(DELTA_MS, color=BLUE, lw=1.6, ls="--", label="total budget Δ = 5 ms")
    ax1.axhline(2 * DELTA_MS, color=MUTED, lw=1, ls=":")
    ax1.annotate("window budget → 2Δ = 10 ms", (1.15, 9.55), fontsize=8, color=MUTED)
    for rv, txt in ((2.0, "ρ=2: R=1.667 ms"), (4.0, "ρ=4: R=3.000 ms")):
        y = DELTA_MS * (rv - 1) / (rv + 1)
        ax1.scatter([rv], [y], s=56, color=AQUA, zorder=5,
                    edgecolor=SURFACE, linewidth=1.5)
        ax1.annotate(txt, (rv + 0.08, y - 0.45), fontsize=8, color=MUTED)
    ax1.set_xlabel("oversampling factor ρ = N/H")
    ax1.set_ylabel("ms")
    ax1.set_title("Oversampling relaxation of the window-length constraint (Δ = 5 ms)")
    ax1.set_ylim(0, 10.6)
    ax1.legend(loc="lower right", fontsize=8)
    style(ax1)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def fig_phase(phase_csv, out_png):
    ph, near, onset = [], [], []
    with open(phase_csv, newline="") as f:
        for r in csv.DictReader(f):
            ph.append(float(r["phase_frac"]))
            near.append(float(r["nearest_ms"]))
            onset.append(float(r["onset_ms"]))
    ph, near, onset = map(np.array, (ph, near, onset))

    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    ax.plot(ph, near, color=BLUE, lw=2, label="A nearest-center deviation")
    ax.plot(ph, onset, color=ORANGE, lw=2, label="B onset-frame center deviation")
    ax.axhline(0, color="#c3c2b7", lw=1)
    ax.axhline(128 / (2 * FS) * 1000, color=BLUE, lw=1, ls="--",
               label="grid bound H/(2fs) = 1.333 ms")
    ax.axhline(256 / (2 * FS) * 1000, color=ORANGE, lw=1, ls="--",
               label="smear bound N/(2fs) = 2.667 ms")
    ax.set_xlabel("event phase within one hop cell (fraction of H)")
    ax.set_ylabel("stamp − true event (ms)")
    ax.set_title("Phase sweep, chosen config N=256 H=128 (ρ=2): signed deviation")
    ax.legend(loc="lower center", fontsize=8, ncol=2)
    style(ax)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def fig_budget(rows, out_png):
    tags = [r["tag"].split()[0] for r in rows]
    grid = np.array([r["grid_ms"] for r in rows])
    smear = np.array([r["smear_ms"] for r in rows])
    meas = np.maximum(np.array([r["meas_nearest_ms"] for r in rows]),
                      np.array([r["meas_onset_ms"] for r in rows]))
    x = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(9.4, 4.8))
    b1 = ax.bar(x - 0.2, grid, width=0.38, color=BLUE, label="grid term H/(2fs)")
    b2 = ax.bar(x - 0.2, smear, width=0.38, bottom=grid, color=ORANGE,
                label="smear term N/(2fs)")
    ax.bar(x + 0.2, meas, width=0.38, color=AQUA,
           label="measured worst deviation")
    ax.axhline(DELTA_MS, color=RED, lw=2, label="budget Δ = 5 ms")
    for xi, r in zip(x, rows):
        c = AQUA if r["verdict"] == "PASS" else RED
        ax.text(xi - 0.2, r["budget_ms"] + 0.18, f'{r["budget_ms"]:.2f}',
                ha="center", fontsize=8, color=MUTED)
        ax.text(xi + 0.2, meas[xi] + 0.18, f'{meas[xi]:.2f}',
                ha="center", fontsize=8, color=MUTED)
        ax.text(xi, -0.9, r["verdict"], ha="center", fontsize=9,
                color=c, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([f'{t}\nN={r["N"]},H={r["H"]},ρ={r["rho"]:g}'
                        for t, r in zip(tags, rows)], fontsize=8)
    ax.set_ylabel("ms")
    ax.set_ylim(-1.4, 11.4)
    ax.set_title("Analytic budget (stacked) vs measured worst deviation, 6 configs")
    ax.legend(loc="upper left", fontsize=8)
    style(ax)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main(art):
    rows = load_summary(os.path.join(art, "summary.csv"))
    fig_region(rows, os.path.join(art, "analytic_region.csv"),
               os.path.join(art, "fig_region.png"))
    fig_relax(os.path.join(art, "fig_relax.png"))
    fig_phase(os.path.join(art, "phase_sweep.csv"),
              os.path.join(art, "fig_phase.png"))
    fig_budget(rows, os.path.join(art, "fig_budget.png"))
    print("figures written to", art)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "../artifacts")
