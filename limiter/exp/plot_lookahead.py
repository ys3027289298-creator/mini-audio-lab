#!/usr/bin/env python3
"""Plot lookahead study results: retention curve + latency, 3 waveforms."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = sys.argv[1] if len(sys.argv) > 1 else "."

sw = np.genfromtxt(f"{d}/sweep.csv", delimiter=",", names=True)

# fit retention(L) = r_inf + (r0 - r_inf) * exp(-L/tau)
from scipy.optimize import curve_fit
def model(L, r0, r_inf, tau):
    return r_inf + (r0 - r_inf) * np.exp(-L / tau)
try:
    p, _ = curve_fit(model, sw["lookahead_ms"], sw["retention"],
                     p0=[1.0, 0.8, 3.0], maxfev=10000)
    fit_txt = f"r(L) = {p[1]:.3f} + {p[0]-p[1]:.3f}·e^(-L/{p[2]:.2f}ms)"
    print("fit:", fit_txt)
except Exception:
    p, fit_txt = None, "fit failed"

fig, ax1 = plt.subplots(figsize=(8, 5))
ax1.plot(sw["lookahead_ms"], sw["retention"], "o-", color="tab:blue", label="retention (measured)")
if p is not None:
    Lf = np.linspace(0, sw["lookahead_ms"].max(), 200)
    ax1.plot(Lf, model(Lf, *p), "--", color="tab:blue", alpha=0.5, label=fit_txt)
ax1.axhline(0.60, color="tab:red", ls=":", label="acceptance limit r >= 0.60")
ax1.set_xlabel("lookahead L (ms)")
ax1.set_ylabel("transient peak retention", color="tab:blue")
ax1.tick_params(axis="y", labelcolor="tab:blue")
ax1.set_ylim(0.5, 1.0)
ax2 = ax1.twinx()
ax2.plot(sw["lookahead_ms"], sw["latency_ms"], "s--", color="tab:orange", label="latency")
ax2.set_ylabel("added latency (ms)", color="tab:orange")
ax2.tick_params(axis="y", labelcolor="tab:orange")
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, loc="center right", fontsize=9)
ax1.set_title("Transient peak retention & latency vs lookahead length (48 kHz, attack 3 ms)")
fig.tight_layout()
fig.savefig(f"{d}/retention_curve.png", dpi=140)
plt.close(fig)

fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
for ax, L in zip(axes, ["0.5", "3.0", "12.0"]):
    w = np.genfromtxt(f"{d}/wave_{L}ms.csv", delimiter=",", names=True)
    ax.plot(w["t_ms"], w["input"], color="0.6", lw=0.8, label="input")
    ax.plot(w["t_ms"], w["output"], color="tab:blue", lw=1.0, label="output")
    ax.axhline(0.8, color="tab:red", ls=":", lw=0.8)
    ax.axhline(-0.8, color="tab:red", ls=":", lw=0.8)
    ax.set_ylabel("amp")
    ax.set_title(f"lookahead = {L} ms", fontsize=10)
    ax.legend(fontsize=8, loc="upper right")
axes[-1].set_xlabel("time (ms)")
fig.suptitle("Drum transient through lookahead limiter (ceiling 0.8)")
fig.tight_layout()
fig.savefig(f"{d}/waves_3lengths.png", dpi=140)
plt.close(fig)
print("plots written")
