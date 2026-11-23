#!/usr/bin/env python3
"""Render the delay-estimation figures from ../artifacts/delay/measurements.txt.

Run `make csv` in this directory first.  Writes one light + one dark SVG (and a
PNG twin for review) per figure into ../artifacts/delay/.

Palette: the validated reference instance, categorical slots 1-3.  Blue alone
carries the estimator, orange carries the integer/coarse grid, gray carries a
reference floor.  Direct labels accompany every series (secondary encoding), as
the light-mode aqua slot needs.
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "artifacts", "delay"))
SRC = os.path.join(OUT, "measurements.txt")

# --- palette (reference instance; validated with validate_palette.js) ------
LIGHT = dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
             grid="#e1e0d9", axis="#c3c2b7",
             s1="#2a78d6", s2="#eb6834", s3="#1baf7a")
DARK = dict(surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", muted="#898781",
            grid="#2c2c2a", axis="#383835",
            s1="#3987e5", s2="#d95926", s3="#199e70")

FONT = next((c for c in ("DejaVu Sans", "Liberation Sans", "Arial")
             if any(c in f.name for f in font_manager.fontManager.ttflist)),
            "sans-serif")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [FONT],
    "font.size": 10,
    "axes.linewidth": 1.0,
    "svg.fonttype": "none",
})


def read_data(path):
    d = dict(xcorr=[], coarse=[], sweep=[], estimate=None)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split(",")
            if p[0] == "xcorr":
                d["xcorr"].append((float(p[1]), float(p[2])))
            elif p[0] == "coarse":
                d["coarse"].append((int(p[1]), float(p[2])))
            elif p[0] == "estimate":
                d["estimate"] = dict(k0=int(p[1]), tau_est=float(p[2]),
                                     tau_true=float(p[3]))
            elif p[0] == "sweep":
                d["sweep"].append(dict(snr=float(p[1]), crlb=float(p[2]),
                                       rmse=float(p[3]), eff=float(p[4]),
                                       frac=float(p[5]), gross=float(p[6])))
    return d


def style(ax, P, xgrid=True, ygrid=True):
    ax.set_facecolor(P["surface"])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(P["axis"])
        ax.spines[s].set_linewidth(1.0)
    ax.tick_params(colors=P["muted"], labelsize=9, length=0)
    ax.grid(xgrid, axis="x", color=P["grid"], linewidth=1.0, linestyle="-")
    ax.grid(ygrid, axis="y", color=P["grid"], linewidth=1.0, linestyle="-")
    ax.xaxis.set_visible(True)
    ax.set_axisbelow(True)


def title(ax, P, main, sub):
    ax.set_title(main, color=P["ink"], fontsize=12.5, fontweight="bold",
                 loc="left", pad=18)
    ax.text(0.0, 1.02, sub, transform=ax.transAxes, color=P["ink2"],
            fontsize=9.5, ha="left", va="bottom")


def save(fig, stem, mode, P):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(OUT, f"{stem}-{mode}.{ext}"),
                    facecolor=P["surface"], bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)


def fig_xcorr(d, P, mode):
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=120)
    fig.patch.set_facecolor(P["surface"])
    xs = [p[0] for p in d["xcorr"]]
    ys = [p[1] for p in d["xcorr"]]
    tau_true = d["estimate"]["tau_true"]
    ax.plot(xs, ys, color=P["s1"], linewidth=2.0, solid_capstyle="round",
            label="interpolated r(τ)")
    cx = [p[0] for p in d["coarse"]]
    cy = [p[1] for p in d["coarse"]]
    ax.scatter(cx, cy, s=48, color=P["s2"], edgecolor=P["surface"],
               linewidth=2.0, zorder=5, label="integer lags (coarse argmax)")
    ax.axvline(tau_true, color=P["muted"], linewidth=1.0, zorder=1)
    ax.set_xlabel("lag τ (samples)", color=P["ink2"])
    ax.set_ylabel("normalised cross-correlation", color=P["ink2"])
    ax.set_xlim(min(xs), max(xs))
    ax.set_ylim(-0.6, 1.12)
    style(ax, P)
    # direct labels, placed clear of each other
    ax.text(tau_true + 0.12, -0.42, f"true τ = {tau_true:.2f}", color=P["ink2"],
            fontsize=9, ha="left", va="center")
    ax.annotate("peak at 7.37", xy=(tau_true, max(ys)),
                xytext=(tau_true + 0.45, 1.03), color=P["ink2"], fontsize=9,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=P["muted"], linewidth=1.0))
    ax.text(7.08, cy[cx.index(7)] - 0.16, "integer samples", color=P["ink2"],
            fontsize=9, ha="left", va="top")
    leg = ax.legend(loc="upper left", frameon=False, fontsize=9,
                    bbox_to_anchor=(0.0, 1.0))
    for t in leg.get_texts():
        t.set_color(P["ink2"])
    title(ax, P, "Cross-correlation peak locates the delay",
          "the interpolated curve resolves a fraction of a sample; the integer grid cannot")
    save(fig, "xcorr", mode, P)


def fig_fractional(d, P, mode):
    fig, ax = plt.subplots(figsize=(6.6, 4.2), dpi=120)
    fig.patch.set_facecolor(P["surface"])
    est = d["estimate"]
    tau_true = est["tau_true"]
    tau_est = est["tau_est"]
    int_est = float(est["k0"])
    rows = [("true delay", tau_true, P["muted"]),
            ("integer argmax", int_est, P["s2"]),
            ("cross-correlation", tau_est, P["s1"])]
    for y, (name, val, col) in zip([2, 1, 0], rows):
        ax.scatter([val], [y], s=120, color=col, edgecolor=P["surface"],
                   linewidth=2.0, zorder=4)
        ax.text(val, y + 0.24, f"{val:.3f}", color=P["ink"], fontsize=10,
                ha="center", va="bottom")
    ax.errorbar([int_est], [1], xerr=[[0.5], [0.5]], fmt="none",
                ecolor=P["s2"], elinewidth=1.5, capsize=5, alpha=0.6, zorder=2)
    ax.text(int_est, 1 - 0.34, "±0.5 sample rounding floor", color=P["ink2"],
            fontsize=8.5, ha="center", va="top")
    ax.text(6.52, 2, "true delay", color=P["ink2"], fontsize=9.5, ha="right", va="center")
    ax.text(6.52, 1, "integer argmax", color=P["ink2"], fontsize=9.5, ha="right", va="center")
    ax.text(6.52, 0, "cross-correlation", color=P["ink2"], fontsize=9.5, ha="right", va="center")
    ax.set_yticks([])
    ax.set_ylim(-0.7, 2.7)
    ax.set_xlim(6.9, 7.9)
    ax.set_xlabel("estimated delay τ (samples)", color=P["ink2"])
    style(ax, P, ygrid=False)
    err_int = abs(int_est - tau_true)
    err_frac = abs(tau_est - tau_true)
    title(ax, P, "Fractional estimate beats integer rounding",
          f"error {err_int:.3f} → {err_frac:.3f} samples: "
          "no longer quantised to the sample grid")
    save(fig, "fractional", mode, P)


def fig_efficiency(d, P, mode):
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=120)
    fig.patch.set_facecolor(P["surface"])
    sw = [s for s in d["sweep"] if s["frac"] >= 0.8]      # bound is meaningful here
    snr = [s["snr"] for s in sw]
    eff = [s["eff"] for s in sw]
    ax.plot(snr, eff, color=P["s1"], linewidth=2.0, solid_capstyle="round", zorder=3)
    ax.scatter(snr, eff, s=48, color=P["s1"], edgecolor=P["surface"],
               linewidth=2.0, zorder=4)
    ax.axhline(1.0, color=P["muted"], linewidth=1.0, zorder=1)
    ax.set_ylim(0.6, 1.25)
    ax.set_xlabel("SNR (dB)", color=P["ink2"])
    ax.set_ylabel("efficiency  CRLB / measured variance", color=P["ink2"])
    style(ax, P)
    ax.text(max(snr), 1.02, "efficiency = 1 (CRLB achieved)", color=P["ink2"],
            fontsize=9, ha="right", va="bottom")
    ax.text(min(snr) + 0.2, min(eff) - 0.05, f"{min(eff):.2f}", color=P["ink"],
            fontsize=9.5, ha="left", va="top")
    title(ax, P, "Estimator reaches the CRLB across SNR",
          "measured variance within ~10% of the bound; plotted only where ≥80% of trials detect")
    save(fig, "efficiency", mode, P)


def fig_detection(d, P, mode):
    """Error vs SNR (log), plus the fraction of trials that are measurable."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.6, 4.0), dpi=120,
                                  gridspec_kw=dict(width_ratios=[1.25, 1.0], wspace=0.42))
    fig.patch.set_facecolor(P["surface"])

    ok = [s for s in d["sweep"] if not math.isnan(s["rmse"]) and s["frac"] >= 0.8]
    snr = [s["snr"] for s in ok]
    ax.plot(snr, [s["crlb"] for s in ok], color=P["muted"], linewidth=2.0,
            zorder=2, label="CRLB σ (lower bound)")
    ax.plot(snr, [s["rmse"] for s in ok], color=P["s1"], linewidth=2.0,
            zorder=3, label="measured RMSE")
    ax.scatter(snr, [s["rmse"] for s in ok], s=46, color=P["s1"],
               edgecolor=P["surface"], linewidth=2.0, zorder=4)
    ax.set_yscale("log")
    ax.set_xlabel("SNR (dB)", color=P["ink2"])
    ax.set_ylabel("delay error (samples, log)", color=P["ink2"])
    style(ax, P)
    leg = ax.legend(loc="upper right", frameon=False, fontsize=8.5)
    for t in leg.get_texts():
        t.set_color(P["ink2"])

    # detection fraction: the "unmeasurable" boundary made explicit
    sw = d["sweep"]
    ax2.plot([s["snr"] for s in sw], [100 * s["frac"] for s in sw],
             color=P["s1"], linewidth=2.0, solid_capstyle="round", zorder=3)
    ax2.scatter([s["snr"] for s in sw], [100 * s["frac"] for s in sw],
                s=42, color=P["s1"], edgecolor=P["surface"], linewidth=2.0, zorder=4)
    ax2.set_ylim(-4, 104)
    ax2.set_xlabel("SNR (dB)", color=P["ink2"])
    ax2.set_ylabel("trials measured (%)", color=P["ink2"])
    style(ax2, P)
    # below the crossing the estimate is refused, not merely noisy
    threshold = next((s["snr"] for s in sw if s["frac"] < 0.5), None)
    if threshold is not None:
        ax2.axvspan(min(s["snr"] for s in sw) - 1, threshold,
                    color=P["grid"], alpha=0.5, zorder=0, linewidth=0)
        ax2.axvline(threshold, color=P["s2"], linewidth=1.5, zorder=2)
        ax2.text(threshold - 0.6, 50, "half detectable", color=P["ink2"],
                 fontsize=8.5, ha="right", va="center", rotation=90)
        ax2.text(min(s["snr"] for s in sw) + 0.4, 12, "unmeasurable",
                 color=P["ink2"], fontsize=8.5, ha="left", va="bottom")
    title(ax, P, "Measured error sits on the CRLB",
          "both fall as 1/√SNR; shown where ≥80% of trials detect")
    ax2.set_title("Detectability gate", color=P["ink"], fontsize=10.5,
                  fontweight="bold", loc="left", pad=8)
    save(fig, "detection", mode, P)


def main():
    if not os.path.exists(SRC):
        sys.exit(f"no data in {SRC}; run `make csv` here first")
    d = read_data(SRC)
    for mode, P in (("light", LIGHT), ("dark", DARK)):
        fig_xcorr(d, P, mode)
        fig_fractional(d, P, mode)
        fig_efficiency(d, P, mode)
        fig_detection(d, P, mode)
    sw = [s for s in d["sweep"] if not math.isnan(s["eff"]) and s["frac"] >= 0.8]
    print(f"wrote figures to {OUT}")
    print("efficiency: " + ", ".join(f"{s['snr']:.0f}dB={s['eff']:.3f}" for s in sw))
    thr = next((s["snr"] for s in d["sweep"] if s["frac"] < 0.5), None)
    print(f"half-detectable SNR ≈ {thr} dB")


if __name__ == "__main__":
    main()
