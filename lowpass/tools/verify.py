#!/usr/bin/env python3
"""Independent cross-check of the C lowpass measurements.

Recomputes the 4th-order Butterworth design with scipy (bilinear transform,
not sharing any code with src/lowpass.c), re-derives both objectives for the
same cutoff sweep, and compares against sweep.csv written by the C binary.
Run: python3 -I tools/verify.py
"""
import csv
import os
import sys

import numpy as np
from scipy import signal

FS = 1000.0
NFFT = 1024
NBIN = NFFT // 2
F_SIG = 125.0
F_ALIAS = 250.0
HERE = os.path.dirname(os.path.abspath(__file__))


def objectives(fc_frac):
    # scipy's own Butterworth design -> digital biquads via bilinear transform.
    # scipy normalizes Wn to the Nyquist rate, while fc_frac here is a
    # fraction of fs, hence the factor of 2.
    b, a = signal.butter(4, fc_frac * 2.0, btype="low", analog=False)
    h = np.zeros(NFFT)
    h[0] = 1.0
    h = signal.lfilter(b, a, h)
    H = np.fft.fft(h, NFFT)
    mag = np.abs(H)
    phase = np.angle(H)

    b_alias_lo = int(F_ALIAS / FS * NFFT) + 1
    e_alias = np.sum(mag[b_alias_lo:NBIN] ** 2)
    alias_db = 10.0 * np.log10(e_alias / mag[0] ** 2 + 1e-20)

    b_hi = int(F_SIG / FS * NFFT)
    unw = np.unwrap(phase[: b_hi + 1])
    w = mag[: b_hi + 1] ** 2
    bins = np.arange(b_hi + 1, dtype=float)
    # weighted least squares with w entering linearly (same normal equations
    # as the C code); np.polyfit would square the weights instead
    sw = w.sum()
    sf = np.sum(w * bins)
    sp = np.sum(w * unw)
    sff = np.sum(w * bins * bins)
    sfp = np.sum(w * bins * unw)
    slope = (sw * sfp - sf * sp) / (sw * sff - sf * sf)
    p0 = (sp - slope * sf) / sw
    resid = unw - (slope * bins + p0)
    phase_rad = np.sqrt(np.sum(w * resid ** 2) / sw)
    return alias_db, phase_rad


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, "..", "sweep.csv"))))
    max_a, max_p, n = 0.0, 0.0, 0
    for r in rows:
        fc = float(r["fc_frac"])
        a_ref, p_ref = objectives(fc)
        a_c = float(r["alias_residual_db"])
        p_c = float(r["phase_distortion_rad"])
        max_a = max(max_a, abs(a_ref - a_c))
        max_p = max(max_p, abs(p_ref - p_c))
        n += 1

    # Pareto re-ranking from the CSV alone
    pts = [(float(r["alias_residual_db"]),
            float(r["phase_distortion_rad"]),
            int(r["pareto"])) for r in rows]
    mism = 0
    for i, (ai, pi, flag) in enumerate(pts):
        dom = any(j != i and aj <= ai and pj <= pi and (aj < ai or pj < pi)
                  for j, (aj, pj, _) in enumerate(pts))
        if (not dom) != bool(flag):
            mism += 1

    print(f"rows={n} max_alias_diff_db={max_a:.4f} "
          f"max_phase_diff_rad={max_p:.2e} pareto_mismatches={mism}")
    ok = max_a < 1e-3 and max_p < 1e-6 and mism == 0
    print("verify_ok=%d" % ok)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
