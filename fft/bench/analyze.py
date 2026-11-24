#!/usr/bin/env python3
"""Analyze the block-length vs. error/latency trade-off.

Reads bench/results_blocking1.csv, fits the two mutually constraining scaling
laws (block error decreasing with block length, latency increasing), builds the
trade-off curve, evaluates the three candidate block lengths and emits:
  - bench/scaling_fit.txt   fitted laws and R^2
  - bench/tradeoff.png      four-panel static figure
  - bench/candidates.csv    three-point evaluation table
  - bench/timing.csv        two sample-rate timing summary
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# dataviz reference palette (light surface)
BLUE, ORANGE, AQUA, VIOLET, RED = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e34948"
SURFACE, GRID, INK, MUTED = "#fcfcfb", "#e1e0d9", "#0b0b0b", "#898781"
CANDIDATES = [128, 512, 1024]

rows = list(csv.DictReader(open("bench/results_blocking1.csv")))


def load(rate):
    B, err, lat, peak, resolved, fft_us, total = [], [], [], [], [], [], []
    for r in rows:
        if int(r["rate_hz"]) != rate or r["mode"] != "block":
            continue
        B.append(int(r["block"])); err.append(float(r["band_err_db"]))
        lat.append(float(r["latency_ms"])); peak.append(float(r["peak_err_hz"]))
        resolved.append(int(r["resolved"])); fft_us.append(float(r["fft_us"]))
        total.append(float(r["total_ms"]))
    return map(np.array, (B, np.array(err), lat, peak, resolved, fft_us, total))


def fit_error_dB(B, err_db):
    """Fit the dB-domain law E_dB = c + s*log10(B). This is the natural form
    for an error made of leakage sidelobes (decaying ~1/B in linear power),
    and corresponds to the power law E_lin = 10^(c/10) * B^(s/10).
    Returns (s_dB_per_decade, c, r2_in_dB, predicted_dB)."""
    lb = np.log10(B)
    s, c = np.polyfit(lb, err_db, 1)
    yhat = s * lb + c
    ss_res = np.sum((err_db - yhat) ** 2); ss_tot = np.sum((err_db - np.mean(err_db)) ** 2)
    return s, c, 1 - ss_res / ss_tot, yhat


def fit_power_law(x, y):
    """y = a * x^p via log-log least squares; returns (a, p, r2, yhat)."""
    lx, ly = np.log(x), np.log(y)
    p, lna = np.polyfit(lx, ly, 1)
    yhat = np.exp(lna) * x ** p
    ss_res = np.sum((y - yhat) ** 2); ss_tot = np.sum((y - np.mean(y)) ** 2)
    return np.exp(lna), p, 1 - ss_res / ss_tot, yhat


report, tables, timing = [], [], []
fig, axes = plt.subplots(2, 3, figsize=(18, 9.5), facecolor=SURFACE)
for ax in axes.flat:
    ax.set_facecolor(SURFACE); ax.grid(True, color=GRID, lw=0.8)
    for s in ax.spines.values(): s.set_color("#c3c2b7")
    ax.tick_params(colors=MUTED); ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)
    ax.title.set_color(INK)

for idx, sr in enumerate((48000, 96000)):
    B, err, lat, peak, resolved, fft_us, total = load(sr)

    # --- scaling law 1: error decreases with block length (fit in dB) ---
    s_e, c_e, r2_e, err_hat = fit_error_dB(B, err)
    p_lin = s_e / 10.0                                   # implied E_lin = a*B^p exponent
    # --- scaling law 2: latency increases with block length ---
    a_l, p_l, r2_l, lat_hat = fit_power_law(B, lat)

    report.append(f"=== {sr} Hz ===")
    report.append(f"error   law: E_dB(B) = {c_e:.3f} + ({s_e:.3f})*log10(B)   R2(dB)={r2_e:.4f}")
    report.append(f"          => linear-power form E = {10**(c_e/10):.3e} * B^({p_lin:.3f})")
    report.append(f"latency law: L(B) = {a_l:.3e} * B^({p_l:.3f})   R2={r2_l:.4f}  "
                  f"(theory L = B/{sr} = {1000.0/sr:.4f} ms * B, slope exponent 1.0)")
    # direct linear latency fit for the exact per-sample coefficient
    k_lat = np.polyfit(B, lat, 1)
    report.append(f"latency linear: L = {k_lat[0]:.5f}*B + {k_lat[1]:.4f} ms "
                  f"(expected {1000.0/sr:.5f} ms/sample)")

    # panel A: error decrease law (log x, dB y)
    ax = axes[idx, 0]
    ax.semilogx(B, err, "o-", color=ORANGE, lw=2, ms=8, label="block error E (dB)")
    ax.semilogx(B, err_hat, "--", color=ORANGE, lw=1.5, alpha=0.7,
                label=f"fit: {s_e:.1f} dB/decade (R²={r2_e:.3f})")
    ax.set_xlabel("block length B (samples)"); ax.set_ylabel("error (dB RMS band-power)")
    ax.set_title(f"{sr} Hz — error decrease law"); ax.legend(fontsize=9, framealpha=0.9)

    # panel B: latency increase law (log-log, single axis)
    ax = axes[idx, 1]
    ax.loglog(B, lat, "s-", color=BLUE, lw=2, ms=8, label="latency L (ms)")
    ax.loglog(B, lat_hat, "--", color=BLUE, lw=1.5, alpha=0.7,
              label=f"fit L∝B$^{{{p_l:.2f}}}$ (R²={r2_l:.3f})")
    ax.set_xlabel("block length B (samples)"); ax.set_ylabel("latency (ms)")
    ax.set_title(f"{sr} Hz — latency increase law"); ax.legend(fontsize=9, framealpha=0.9)

    # panel C: trade-off frontier error vs latency (the actual design curve)
    ax = axes[idx, 2]
    order = np.argsort(lat)
    ax.plot(lat[order], err[order], "o-", color=VIOLET, lw=2, ms=8)
    n_pts = len(B)
    for j, (b, x_, y_) in enumerate(zip(B, lat, err)):
        # last (rightmost) label goes left to avoid clipping; others alternate above/below
        if j == n_pts - 1:
            off = (-46, 4)
        else:
            off = (8, 6 if j % 2 == 0 else -11)
        ax.annotate(f"B={b}", (x_, y_), textcoords="offset points", xytext=off,
                    fontsize=8, color=MUTED)
    for cb, col in zip(CANDIDATES, (AQUA, RED, ORANGE)):
        m = B == cb
        ax.scatter(lat[m], err[m], s=220, facecolors="none", edgecolors=col, lw=2.5, zorder=5)
    ax.set_xscale("log"); ax.set_xlabel("latency (ms)"); ax.set_ylabel("error (dB)")
    ax.set_title(f"{sr} Hz — trade-off frontier (rings = candidates)")
    from matplotlib.ticker import FixedLocator, FixedFormatter
    xt = np.unique(np.round(lat, 1))
    ax.xaxis.set_major_locator(FixedLocator(xt))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{v:g}" for v in xt]))
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.tick_params(axis="x", labelsize=8)

    # candidate table + timing rows
    for cb in CANDIDATES:
        m = np.nonzero(B == cb)[0][0]
        tables.append([sr, cb, round(float(lat[m]), 3), round(float(err[m]), 3),
                       round(float(peak[m]), 3), "yes" if resolved[m] else "no",
                       round(float(fft_us[m]), 3), round(float(total[m]), 3)])
    for b, fu in zip(B, fft_us):
        # 50% overlap => one transform every hop=B/2 samples, i.e. fs/hop per second
        us_per_audio_sec = fu * (sr / (b / 2.0))
        timing.append([sr, b, round(1000.0 * b / sr, 3), round(float(fu), 3),
                       round(float(us_per_audio_sec), 1),
                       round(float(us_per_audio_sec) / 1e6, 6)])

# write scaling report
open("bench/scaling_fit.txt", "w").write("\n".join(report) + "\n")

# candidate evaluation CSV
with open("bench/candidates.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rate_hz", "block", "latency_ms", "band_err_db", "peak_err_hz",
                "pair_resolved", "fft_us", "total_ms"])
    w.writerows(tables)

# timing CSV: fft_us = per-transform cost; cpu_per_s = transforms/sec * us = us of FFT per audio second
with open("bench/timing.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rate_hz", "block", "latency_ms", "fft_us_per_transform",
                "fft_us_per_audio_second", "cpu_fraction"])
    w.writerows(timing)

fig.suptitle("Blocked FFT spectrum: error-latency trade-off "
             "(Hann, 50%-overlap Welch vs full-segment FFT reference)",
             fontsize=14, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("bench/tradeoff.png", dpi=140, facecolor=SURFACE)
print("\n".join(report))
print("\nwrote bench/scaling_fit.txt, candidates.csv, timing.csv, tradeoff.png")
