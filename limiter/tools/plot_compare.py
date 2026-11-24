#!/usr/bin/env python3
"""Plot old vs new limiter gain curves and print the segment table.

Reads out/material_gain.csv and out/impulse_gain.csv (written by ./compare),
writes out/gain_compare.png and out/segments.md.
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

mat = np.genfromtxt("out/material_gain.csv", delimiter=",", names=True)
imp = np.genfromtxt("out/impulse_gain.csv", delimiter=",", names=True)

fig, (axm, axi) = plt.subplots(2, 1, figsize=(11, 7.5), height_ratios=[2, 1])

axm.plot(mat["t"], mat["gain_old"], color="#d9544d", lw=1.0,
         label="old: single release 300 ms")
axm.plot(mat["t"], mat["gain_new"], color="#1f77b4", lw=1.0,
         label="new: fast 60 ms / slow 400 ms, split 0.25")
axenv = axm.twinx()
axenv.plot(mat["t"], np.abs(mat["in"]), color="gray", lw=0.4, alpha=0.35)
axenv.set_ylabel("|input|", color="gray")
axenv.set_ylim(0, 1.3)
axenv.tick_params(axis="y", labelcolor="gray")
axm.axvline(0.4, color="k", ls=":", lw=0.8)
axm.text(0.405, 0.05, "loud -> quiet", fontsize=8, color="k")
axm.set_ylabel("gain envelope")
axm.set_ylim(0, 1.05)
axm.set_title("limiter gain: old vs envelope-split release (same material)")
axm.legend(loc="lower right")

axi.plot(imp["t"], imp["gain_old"], color="#d9544d", lw=1.0, label="old")
axi.plot(imp["t"], imp["gain_new"], color="#1f77b4", lw=1.0, label="new")
axi.axhline(1.0, color="k", ls="--", lw=0.8)
axi.set_xlabel("time (s)")
axi.set_ylabel("gain envelope")
axi.set_title("pure impulse @ 10 ms (zoomed: overshoot would show above 1.0)")
axi.set_ylim(0.999, 1.0005)
axi.legend(loc="lower right")

fig.tight_layout()
fig.savefig("out/gain_compare.png", dpi=110)

# segment table over the material
seg = 0.2
rows = []
t = mat["t"]
t0 = 0.0
while t0 < t[-1] - 1e-9:
    m = (t >= t0) & (t < t0 + seg)
    rows.append((t0, min(t0 + seg, t[-1]),
                 mat["gain_old"][m].mean(), mat["gain_new"][m].mean(),
                 mat["gain_old"][m].min(), mat["gain_new"][m].min()))
    t0 += seg

lines = [
    "| 时段 (s) | 旧均值 | 新均值 | 旧最小 | 新最小 |",
    "|---|---|---|---|---|",
]
for r in rows:
    lines.append("| %.1f–%.1f | %.3f | %.3f | %.3f | %.3f |" % r)
table = "\n".join(lines)
with open("out/segments.md", "w") as fh:
    fh.write(table + "\n")
print(table)

print("impulse max gain: old=%.6f new=%.6f (overshoot if > 1.0)"
      % (imp["gain_old"].max(), imp["gain_new"].max()))
print("wrote out/gain_compare.png, out/segments.md")
