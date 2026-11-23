#!/usr/bin/env python3
"""Render the doppler figures from the CSVs under artifacts/doppler/.

Reads sweep.csv and alias.csv, writes trajectory.svg and alias.svg next to
them. Colors are the validated light-mode categorical slots 1 (blue, slow)
and 2 (orange, fast) from the dataviz reference palette.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = {"slow": "#2a78d6", "fast": "#eb6834"}


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelsize=8, length=3, width=0.8)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def trajectory(rows, out):
    fig, ax = plt.subplots(figsize=(8.0, 3.7), dpi=170)
    fig.patch.set_facecolor(SURFACE)
    style(ax)

    for prof, label in (("slow", "slow  v = 30 m/s, a = 0"),
                        ("fast", "fast  v = 120 m/s, a = 20 m/s²")):
        r = [x for x in rows if x["profile"] == prof]
        t = [float(x["t"]) for x in r]
        f = [float(x["f_analytic"]) for x in r]
        ax.plot(t, f, color=SERIES[prof], linewidth=2.0, label=label,
                solid_capstyle="round")
        # direct label at the end of each trace
        ax.annotate(label.split()[0], xy=(t[-1], f[-1]), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=8,
                    color=INK2)

    ax.axhline(440.0, color=MUTED, linewidth=1.0, linestyle=(0, (4, 3)))
    ax.annotate("rest pitch f₀ = 440 Hz", xy=(0.01, 440.0), xycoords=("axes fraction", "data"),
                xytext=(0, 4), textcoords="offset points", fontsize=8, color=MUTED)
    ax.axvline(0.0, color=AXIS, linewidth=0.8, linestyle=(0, (2, 3)))

    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(300, 690)
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.set_xlabel("listener time t (s)", fontsize=9, color=INK2)
    ax.set_ylabel("received frequency (Hz)", fontsize=9, color=INK2)
    ax.set_title("Instantaneous received frequency, analytic f = f₀c/(c + vᵣ)",
                 fontsize=10, color=INK, loc="left", pad=10)
    leg = ax.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK2,
                    ncol=2, handlelength=1.6, columnspacing=1.2)
    for line in leg.get_lines():
        line.set_linewidth(2.0)
    fig.tight_layout(pad=0.6)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)


def alias(rows, out):
    fig, ax = plt.subplots(figsize=(8.0, 3.3), dpi=170)
    fig.patch.set_facecolor(SURFACE)
    style(ax)

    fs = 1000
    spans = []
    for prof, label in (("slow", "slow, fs = 1 kHz"), ("fast", "fast, fs = 1 kHz")):
        r = [x for x in rows if x["profile"] == prof and int(x["fs"]) == fs]
        t = [float(x["t"]) for x in r]
        s = [float(x["phase_step_pi"]) for x in r]
        ax.plot(t, s, color=SERIES[prof], linewidth=2.0, label=label,
                solid_capstyle="round")
        # contiguous runs where the true step passes pi
        run = None
        for ti, row in zip(t, r):
            if row["aliased"] == "1":
                if run is None:
                    run = [ti, ti]
                else:
                    run[1] = ti
            elif run is not None:
                spans.append((prof, run[0], run[1]))
                run = None
        if run is not None:
            spans.append((prof, run[0], run[1]))

    for prof, a, b in spans:
        ax.axvspan(a, b, color=SERIES[prof], alpha=0.10, linewidth=0)

    ax.axhline(1.0, color=INK2, linewidth=1.0)
    ax.annotate("Nyquist edge: one-sample step = π  (f = fs/2)",
                xy=(0.01, 1.0), xycoords=("axes fraction", "data"),
                xytext=(0, 4), textcoords="offset points", fontsize=8, color=INK2)
    ax.annotate("shaded: aliased (step ≥ π)", xy=(0.02, 0.92), xycoords="axes fraction",
                fontsize=8, color=INK2)

    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0.0, 1.75)
    ax.set_xlabel("listener time t (s)", fontsize=9, color=INK2)
    ax.set_ylabel("phase step per sample / π", fontsize=9, color=INK2)
    ax.set_title("Aliasing boundary at fs = 1 kHz (fast exceeds it, slow stays below)",
                 fontsize=10, color=INK, loc="left", pad=10)
    leg = ax.legend(loc="upper right", frameon=False, fontsize=8, labelcolor=INK2,
                    handlelength=1.6)
    for line in leg.get_lines():
        line.set_linewidth(2.0)
    fig.tight_layout(pad=0.6)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "artifacts/doppler"
    sweep = read_rows(os.path.join(d, "sweep.csv"))
    al = read_rows(os.path.join(d, "alias.csv"))
    trajectory(sweep, os.path.join(d, "trajectory.svg"))
    alias(al, os.path.join(d, "alias.svg"))
    print("wrote", os.path.join(d, "trajectory.svg"))
    print("wrote", os.path.join(d, "alias.svg"))


if __name__ == "__main__":
    main()
