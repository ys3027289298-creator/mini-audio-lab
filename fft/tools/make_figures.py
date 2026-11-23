#!/usr/bin/env python3
"""make_figures.py — render the spectral-averaging evidence plots and the
deviation/SNR tables from the Monte-Carlo CSVs.

Inputs : --in  directory written by spectral_mc.py
Outputs: PNG figures and tables.md into --out.

Run: python3 -I make_figures.py --in <csv dir> --out <fig dir>
"""
import argparse
import csv
import os
import sys

# Allow `python3 -I tools/make_figures.py`: -I drops user site-packages.
try:
    import numpy  # noqa: F401
    import matplotlib  # noqa: F401
except ImportError:
    import site
    sys.path.insert(0, site.getusersitepackages())

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# dataviz validated categorical palette (light surface). Aqua/yellow sit below
# the 3:1 contrast line, so every series is direct-labelled and tables.md
# carries the same numbers.
BLUE, ORANGE, AQUA, YELLOW, RED = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e34948"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "font.size": 11, "axes.titlesize": 12,
    "axes.titleweight": "bold", "legend.frameon": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "axes.spines.top": False, "axes.spines.right": False,
})


def load(d, name):
    with open(os.path.join(d, name)) as f:
        return list(csv.DictReader(f))


def fnum(rows, key):
    return np.array([float(r[key]) for r in rows])


def style(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.tick_params(length=3)


def save(fig, out, name):
    fig.tight_layout()
    fig.savefig(os.path.join(out, name), dpi=110)
    plt.close(fig)


def figure_law(indir, out):
    rows = load(indir, "expA_law.csv")
    ks = sorted({int(r["k"]) for r in rows})
    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    snr_curves = [-10.0, 0.0, 10.0]
    colors = {-10.0: BLUE, 0.0: ORANGE, 10.0: AQUA}
    for snr in snr_curves:
        sub = [r for r in rows if float(r["snr_db"]) == snr]
        sub.sort(key=lambda r: int(r["k"]))
        k = fnum(sub, "k")
        ax.plot(k, fnum(sub, "cv_pred"), "-", color=colors[snr], lw=2,
                label=f"theory {snr:g} dB")
        ax.plot(k, fnum(sub, "cv_meas"), "o", color=colors[snr], ms=6,
                mfc=SURFACE, mew=1.8, label=f"measured {snr:g} dB")
    # reference 1/sqrt(K) guide through the 0 dB first point
    sub0 = [r for r in rows if float(r["snr_db"]) == 0.0]
    sub0.sort(key=lambda r: int(r["k"]))
    k0 = fnum(sub0, "k"); c0 = float(sub0[0]["cv_pred"])
    ax.plot(k0, c0 * np.sqrt(k0[0] / k0), "--", color=MUTED, lw=1.4,
            label="1/√K guide")
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks(ks); ax.set_xticklabels([str(k) for k in ks])
    style(ax, "frames K (log2)", "jitter  CV = σ/μ (log)",
          "Fig.1  Power-average jitter vs frame count — theory C(ρ)/K vs measured")
    ax.legend(loc="upper right", ncol=2, fontsize=9)
    save(fig, out, "fig1_law.png")


def figure_refined(indir, out):
    rows = load(indir, "expA2_refined_law.csv")
    rows.sort(key=lambda r: float(r["snr_db"]))
    snr = fnum(rows, "snr_db")
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    ax.plot(snr, fnum(rows, "cv_meas"), "o", color=INK, ms=7,
            label="measured CV (K=32, 50% overlap)")
    ax.plot(snr, fnum(rows, "j_const15"), "s", color=ORANGE, ms=7,
            label="textbook law: constant IF = 1.5")
    ax.plot(snr, fnum(rows, "j_refined"), "^", color=BLUE, ms=8,
            label="refined law: bin-independent IF(ρ), g=1/6")
    style(ax, "frame SNR (dB)", "jitter CV",
          "Fig.2  Constant-IF law over-predicts; the ρ-dependent law tracks measurement")
    ax.legend(loc="upper right", fontsize=9)
    save(fig, out, "fig2_refined_law.png")


def figure_overlap(indir, out):
    g_rows = load(indir, "expB2_perbin_g.csv")
    b = fnum(g_rows, "bin").astype(int)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 5.2))
    ax1.plot(b, fnum(g_rows, "gamma_pred"), "-", color=BLUE, lw=2,
             label="theory: flat g² = 1/36 (interior bins)")
    edge = np.array([r["edge"] == "True" for r in g_rows])
    ax1.plot(b[~edge], fnum(g_rows, "gamma_meas")[~edge], "o", color=ORANGE,
             ms=5, mfc=SURFACE, mew=1.6, label="measured interior bins")
    ax1.plot(b[edge], fnum(g_rows, "gamma_meas")[edge], "D", color=YELLOW,
             ms=7, mew=1.4, label="DC/Nyquist (single-DOF, separate)")
    style(ax1, "FFT bin", "lag-1 frame correlation γ(b)",
          "Fig.3a  Interior-bin overlap correlation is frequency-flat")
    ax1.legend(fontsize=9)
    # inflation factor vs rho: analytic curve and measured points (bin-independent)
    rows = load(indir, "expA2_refined_law.csv")
    rho = fnum(rows, "rho")
    k = 32
    g = float(g_rows[0]["g"])
    grid = np.geomspace(rho.min() / 2, rho.max() * 2, 200)
    infl = 1.0 + 2.0 * (1 - 1.0 / k) * (g ** 2 + 2.0 * grid * g) / (1 + 2.0 * grid)
    ax2.plot(grid, infl, "-", color=BLUE, lw=2,
             label="refined law IF(ρ), g=1/6")
    ax2.plot(rho, fnum(rows, "if_refined"), "o", color=ORANGE, ms=7,
             mfc=SURFACE, mew=1.7, label="IF implied by measured CV")
    ax2.axhline(1.5, color=MUTED, ls="--", lw=1.4, label="textbook bound IF=1.5")
    ax2.axhline(4.0 / 3.0, color=AQUA, ls=":", lw=1.6,
                label="high-ρ asymptote IF→4/3")
    ax2.set_xscale("log")
    ax2.set_ylim(1.0, 1.6)
    style(ax2, "local bin SNR ρ (log)", "inflation factor IF(ρ)",
          "Fig.3b  Inflation rises from 1+g²=1.028 to the 4/3 ceiling")
    ax2.legend(fontsize=8.5)
    save(fig, out, "fig3_overlap.png")


def figure_leakage(indir, out):
    rows = load(indir, "expC1_leakage.csv")
    b = fnum(rows, "bin").astype(int)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 5.2))
    ax1.semilogy(b, np.maximum(fnum(rows, "mean"), 1e-4), "o-", color=BLUE,
                 ms=4, label="measured mean power")
    ax1.semilogy(b, np.maximum(1.0 + fnum(rows, "rho_leak"), 1e-4), "^",
                 color=ORANGE, ms=6, label="analytic 1 + ρ_leak(b)")
    ax1.axhline(1.0, color=MUTED, ls="--", lw=1.3, label="noise floor = 1")
    style(ax1, "FFT bin (tone at 8.5, half-bin off-grid, 0 dB)",
          "mean power (log)", "Fig.4a  Leakage raises the floor across sidelobe bands")
    ax1.legend(fontsize=9)
    # bias vs SNR in four regions
    lrows = load(indir, "expC1b_leakage_snr.csv")
    snr = sorted({float(r["snr_db"]) for r in lrows})
    cols = {"main-lobe": ORANGE, "main-lobe-edge": YELLOW,
            "near-sidelobe": AQUA, "mid-sidelobe": BLUE}
    for region, col in cols.items():
        sub = sorted((r for r in lrows if r["region"] == region),
                     key=lambda r: float(r["snr_db"]))
        ax2.plot([float(r["snr_db"]) for r in sub],
                 [float(r["bias_meas"]) for r in sub], "o-", color=col,
                 ms=5, lw=1.8, label=region)
    style(ax2, "tone SNR (dB)", "leakage bias  (mean−floor)/mean",
          "Fig.4b  Leakage bias vs SNR by band")
    ax2.text(0.03, 0.55, "bias grows ~10^(SNR/10);\nmain lobe saturates first",
             transform=ax2.transAxes, fontsize=9.5, color=MUTED)
    ax2.legend(fontsize=9, loc="upper right")
    save(fig, out, "fig4_leakage.png")


def figure_floors_design(indir, out):
    c3 = load(indir, "expC3_nonstat.csv")
    d = load(indir, "expD_design.csv")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 5.2))
    for depth, col, name in ((0.0, MUTED, "stationary"),
                             (0.25, AQUA, "modulation depth 25%"),
                             (0.5, RED, "modulation depth 50%")):
        sub = sorted((r for r in c3 if float(r["depth"]) == depth),
                     key=lambda r: int(r["k"]))
        k = fnum(sub, "k")
        ax1.plot(k, fnum(sub, "cv_meas"), "o", color=col, ms=5, label=name)
        ax1.plot(k, fnum(sub, "cv_pred_with_floor"), "-", color=col, lw=1.6,
                 alpha=0.85, label=name + " model" if depth else "stationary law")
    ax1.set_xscale("log", base=2); ax1.set_yscale("log")
    style(ax1, "frames K (log2)", "jitter CV (log)",
          "Fig.5a  Non-stationarity floor")
    ax1.text(0.03, 0.92, "jitter stops following 1/√K at the modulation floor",
             transform=ax1.transAxes, fontsize=9.5, color=MUTED)
    ax1.legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    for jmax, col in ((0.20, AQUA), (0.10, BLUE), (0.05, ORANGE)):
        sub = sorted((r for r in d if float(r["j_max"]) == jmax),
                     key=lambda r: float(r["snr_db"]))
        ax2.plot([float(r["snr_db"]) for r in sub],
                 [int(r["k_design"]) for r in sub], "o-", color=col,
                 lw=1.8, ms=6, label=f"budget J_max={jmax:g}")
    ax2.set_yscale("log")
    style(ax2, "frame SNR (dB)", "designed frames K (log)",
          "Fig.5b  Frame count from the design formula")
    ax2.text(0.03, 0.05, r"K = ⌈1.5·C(ρ)/J_max²⌉, all 15 designs within budget",
             transform=ax2.transAxes, fontsize=9.5, color=MUTED)
    ax2.legend(fontsize=9)
    fig.subplots_adjust(right=0.80)
    save(fig, out, "fig5_floors_design.png")


def write_tables(indir, out):
    lines = ["# Spectral averaging — deviation tables (Monte-Carlo evidence)\n"]

    a2 = load(indir, "expA2_refined_law.csv")
    lines.append("## Table 1 — Law departure vs SNR (K=32, 50% overlap, tone bin)\n")
    lines.append("| SNR dB | ρ (local) | measured CV | textbook IF=1.5 | q (vs 1.5) | refined IF(ρ) | q (vs refined) |")
    lines.append("|--:|--:|--:|--:|--:|--:|--:|")
    for r in sorted(a2, key=lambda r: float(r["snr_db"])):
        lines.append("| {snr:g} | {rho:.2f} | {m:.4f} | {c:.4f} | {qc:.3f} | {ir:.3f} | {qr:.3f} |".format(
            snr=float(r["snr_db"]), rho=float(r["rho"]), m=float(r["cv_meas"]),
            c=float(r["j_const15"]), qc=float(r["q_vs_const"]),
            ir=float(r["if_refined"]), qr=float(r["q_vs_refined"])))

    c1b = load(indir, "expC1b_leakage_snr.csv")
    lines.append("\n## Table 2 — Leakage bias vs SNR (half-bin off-grid tone, K=128)\n")
    lines.append("| SNR dB | region (bin) | measured bias | analytic bias | ρ_leak |")
    lines.append("|--:|:--|--:|--:|--:|")
    for r in c1b:
        rl = float(r["rho_leak"])
        rl_s = f"{rl:.4f}" if rl >= 1e-3 else f"{rl:.1e}"
        lines.append("| {snr:g} | {reg} (b{b}) | {m:.4f} | {p:.4f} | {rl} |".format(
            snr=float(r["snr_db"]), reg=r["region"], b=r["bin"],
            m=float(r["bias_meas"]), p=float(r["bias_pred"]), rl=rl_s))

    c2b = load(indir, "expC2b_colored_perbin.csv")
    lines.append("\n## Table 3 — Coloured AR(1) noise per-bin inflation q (K=64)\n")
    lines.append("| φ | low bins 1-4 (mean q) | mid bins 12-18 | high bins 28-32 | max-q bin |")
    lines.append("|--:|--:|--:|--:|:--|")
    for phi in sorted({float(r["phi"]) for r in c2b}):
        sub = [r for r in c2b if float(r["phi"]) == phi]
        def qmean(lo, hi):
            v = [float(r["q"]) for r in sub if lo <= int(r["bin"]) <= hi]
            return sum(v) / len(v)
        worst = max(sub, key=lambda r: float(r["q"]))
        lines.append("| {phi:g} | {lo:.3f} | {mid:.3f} | {hi:.3f} | b{b} ({q:.3f}) |".format(
            phi=phi, lo=qmean(1, 4), mid=qmean(12, 18), hi=qmean(28, 32),
            b=worst["bin"], q=float(worst["q"])))

    c3 = load(indir, "expC3_nonstat.csv")
    lines.append("\n## Table 4 — Non-stationarity floor (0 dB, depth 50%)\n")
    lines.append("| K | measured CV | 1/√K law (IF=1.5) | analytic floor | q vs law | q vs law+floor |")
    lines.append("|--:|--:|--:|--:|--:|--:|")
    for r in [r for r in c3 if float(r["depth"]) == 0.5]:
        lines.append("| {k} | {m:.4f} | {p:.4f} | {fl:.4f} | {q:.3f} | {qf:.3f} |".format(
            k=r["k"], m=float(r["cv_meas"]), p=float(r["cv_pred"]),
            fl=float(r["floor_pred"]), q=float(r["q"]), qf=float(r["q_floor"])))

    c4 = load(indir, "expC4_magnitude.csv")
    lines.append("\n## Table 5 — Magnitude (Rice) vs power averaging (selected SNR)\n")
    lines.append("| SNR dB | K | power CV | magnitude CV | mag/power | magnitude bias |")
    lines.append("|--:|--:|--:|--:|--:|--:|")
    for r in c4:
        if float(r["snr_db"]) in (-15.0, -6.0, 0.0, 10.0) and r["k"] == "256":
            lines.append("| {snr:g} | {k} | {p:.4f} | {mg:.4f} | {q:.3f} | {bias:+.4f} |".format(
                snr=float(r["snr_db"]), k=r["k"], p=float(r["cv_power"]),
                mg=float(r["cv_mag"]), q=float(r["q_mag"]),
                bias=float(r["mag_bias"])))

    d = load(indir, "expD_design.csv")
    lines.append("\n## Table 6 — Frame-count design K=⌈1.5·C(ρ)/J_max²⌉ and measured margin\n")
    lines.append("| J_max | SNR dB | ρ | K design | predicted jitter | measured jitter | margin | verdict |")
    lines.append("|--:|--:|--:|--:|--:|--:|--:|:--|")
    for r in d:
        verdict = "PASS" if float(r["margin"]) > 0 else "FAIL"
        lines.append("| {j:g} | {snr:g} | {rho:.2f} | {k} | {jp:.4f} | {jm:.4f} | {mg:+.4f} | {v} |".format(
            j=float(r["j_max"]), snr=float(r["snr_db"]), rho=float(r["rho"]),
            k=r["k_design"], jp=float(r["j_pred"]), jm=float(r["j_meas"]),
            mg=float(r["margin"]), v=verdict))

    with open(os.path.join(out, "tables.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    figure_law(args.indir, args.out)
    figure_refined(args.indir, args.out)
    figure_overlap(args.indir, args.out)
    figure_leakage(args.indir, args.out)
    figure_floors_design(args.indir, args.out)
    write_tables(args.indir, args.out)
    print("figures and tables ->", args.out)


if __name__ == "__main__":
    main()
