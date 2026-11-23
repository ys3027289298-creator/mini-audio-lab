#!/usr/bin/env python3
"""spectral_mc.py — Welch power-spectrum averaging: Monte-Carlo verification of the
variance-vs-frames law and of its low-SNR departures.

Model
-----
Signal x[n] = tone + Gaussian noise, framed with a Hann window w, periodogram of
each frame P_k(f) = |FFT(w * frame)|^2 / S2, S2 = sum w^2 (unbiased for white
noise). Averaged estimator Pbar = (1/K) sum_k P_k.

Reproduces (see fft/SPECTRAL_AVERAGING.md for the derivations):
  A. law      CV[Pbar] = C(rho)/K   (independent frames, stationary Gaussian)
     C(rho)  = (1 + 2 rho)/(1 + rho)^2, rho = local bin SNR (signal power / noise power)
  B. overlap inflation factor IF for Hann 50%/75%, analytic = 1.5 / 2.015...
  C. low-SNR departures:
       C1 off-bin spectral leakage floor (Hann sidelobes), analytic envelope
       C2 frame correlation from coloured AR(1) noise
       C3 non-stationarity floor (noise level modulates during the block)
       C4 magnitude (Rice) averaging instead of power averaging
  D. frame-count design K_design = ceil(IF * C(rho) / J_max^2) and its measured
     residual jitter vs the budget.

Outputs (CSV + markdown tables + PNG) land in the directory given by --out.

Run:  python3 -I spectral_mc.py --out <dir>
"""
import argparse
import csv
import math
import os
import sys

# Allow `python3 -I tools/spectral_mc.py`: -I drops user site-packages where
# numpy may live, so add it back explicitly (stdlib `site` stays importable).
try:
    import numpy as np  # noqa: F401
except ImportError:
    import site
    sys.path.insert(0, site.getusersitepackages())

import numpy as np

# ----------------------------------------------------------------------------
# Shared simulation constants — kept in lock-step with fft/src/spectral.c
# ----------------------------------------------------------------------------
NFFT = 64
HOP_RATIO_DEFAULT = 0.5          # 50 % overlap
KS = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256], dtype=int)
SNR_DB = np.array([-15.0, -10.0, -6.0, -3.0, 0.0, 3.0, 6.0, 10.0, 20.0], dtype=float)
TRIALS_LAW = 4000                # main-law trials -> ~2.2 % 1-sigma error on CV
TRIALS_FLOOR = 2000              # heavy K=256 sweeps
CHUNK = 250                      # trials processed per chunk to bound memory
RNG_SEED = 20261124

# dataviz validated categorical palette (light surface), direct labels + tables
# provide the relief required for the two sub-3:1 contrast slots.
COL_BLUE = "#2a78d6"
COL_ORANGE = "#eb6834"
COL_AQUA = "#1baf7a"
COL_YELLOW = "#eda100"
COL_RED = "#e34948"
COL_INK = "#0b0b0b"
COL_MUTED = "#898781"
COL_GRID = "#e1e0d9"
COL_SURFACE = "#fcfcfb"


def hann(n):
    """Periodic Hann window (matches the C side: 0.5*(1-cos(2*pi*n/N)))."""
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(n) / n))


W = hann(NFFT)
S2 = float(np.sum(W ** 2))


# ----------------------------------------------------------------------------
# Signal generation
# ----------------------------------------------------------------------------
def gen_frames(rng, n_real, n_frames, nfft, hop, kind, snr_db, tone_bin, ar_phi=0.0,
               mod_rate=0.0, mod_depth=0.0):
    """Generate n_real independent realizations, each n_frames overlapping frames.

    kind: "white" | "ar1" | "nonstat" | "leak" (leak == white + off-grid tone)
    Returns array (n_real, n_frames, nfft).
    """
    total = (n_frames - 1) * hop + nfft
    t = np.arange(total)
    a = math.sqrt(2.0 * 10.0 ** (snr_db / 10.0))   # unit-variance noise reference
    if kind == "ar1":
        eps = rng.standard_normal((n_real, total))
        noise = np.empty((n_real, total))
        noise[:, 0] = eps[:, 0]
        scale = math.sqrt(1.0 - ar_phi ** 2)       # renormalise to unit variance
        for i in range(1, total):
            noise[:, i] = ar_phi * noise[:, i - 1] + eps[:, i]
        noise *= scale
    elif kind == "nonstat":
        env = 1.0 + mod_depth * np.sin(2.0 * np.pi * mod_rate * np.arange(total))
        noise = rng.standard_normal((n_real, total)) * env[None, :]
    else:
        noise = rng.standard_normal((n_real, total))
    tone = a * np.cos(2.0 * np.pi * tone_bin * t / nfft)[None, :]
    x = tone + noise
    idx = np.arange(n_frames)[None, :, None] * hop + np.arange(nfft)[None, None, :]
    return x[:, idx[0]]


def periodograms(frames, window=W, s2=S2):
    """Periodogram bins 0..N/2, scaled by S2. Accepts (frames, nfft) or
    (realizations, frames, nfft)."""
    win = window.reshape((1,) * (frames.ndim - 1) + (-1,))
    spec = np.fft.rfft(frames * win, axis=-1)
    return (spec.real ** 2 + spec.imag ** 2) / s2


def avg_estimate(frames, k, mode="power"):
    """Average of the first k frame periodograms (power or magnitude)."""
    p = periodograms(frames[..., :k, :])
    axis = -2 if p.ndim == 3 else 0
    if mode == "power":
        return p.mean(axis=axis)
    return np.sqrt(p).mean(axis=axis)


# ----------------------------------------------------------------------------
# Theory
# ----------------------------------------------------------------------------
def c_of_rho(rho):
    """Squared relative variance of one-frame periodogram at local SNR rho."""
    return (1.0 + 2.0 * rho) / (1.0 + rho) ** 2


def cv_law(k, rho, inflate=1.0):
    """Predicted CV of the K-frame power average."""
    return np.sqrt(inflate * c_of_rho(np.asarray(rho, dtype=float)) / k)


def design_frames(j_max, rho, inflate):
    """Smallest K (integer) meeting the jitter budget J_max."""
    return int(math.ceil(inflate * c_of_rho(rho) / j_max ** 2))


def hann_dtft_sq(omega, n=NFFT):
    """|DTFT of the periodic Hann window|^2, dense-grid interpolation."""
    w = 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(n) / n))
    v = np.fft.fftshift(np.fft.fft(w, 4096))
    freqs = np.fft.fftshift(2.0 * np.pi * np.fft.fftfreq(4096, d=1.0))
    return np.interp(omega, freqs, np.abs(v) ** 2, period=2.0 * np.pi)


def leakage_envelope(rho_tone, delta, bin_idx):
    """Analytic per-bin leakage fraction L_b / noise power.

    Tone at non-integer bin (bin_idx_of_tone = nearest integer + delta).
    Returns rho_b^leak = tone power delivered into each rFFT bin (one-sided
    scaling) divided by the noise floor of that bin. Noise floor per rFFT bin
    with periodogram scaling is 1.0 (unit-variance input noise).
    """
    n = NFFT
    tone_pos = bin_idx + delta
    omegas = 2.0 * np.pi * (np.arange(n // 2 + 1) - tone_pos) / n
    # one-sided tone coefficient A -> complex amplitude A/2; power A^2/4.
    a_sq_over_4 = 0.5 * 10.0 ** (rho_tone / 10.0)
    w_sq = np.array([hann_dtft_sq(o, n) for o in omegas])
    return a_sq_over_4 * w_sq / S2


# ----------------------------------------------------------------------------
# Monte-Carlo driver
# ----------------------------------------------------------------------------
def measure_cv(rng, kind, k, snr_db, tone_bin, trials, mode="power", **kw):
    """CV (std/|mean|) of the averaged estimate at the tone bin."""
    hop = kw.pop("hop", NFFT // 2)
    vals = np.empty(trials)
    bin_sel = int(round(tone_bin))
    for c0 in range(0, trials, CHUNK):
        m = min(CHUNK, trials - c0)
        frames = gen_frames(rng, m, k, NFFT, hop, kind, snr_db, tone_bin, **kw)
        vals[c0:c0 + m] = avg_estimate(frames, k, mode)[:, bin_sel]
    mean = vals.mean()
    return vals.std(ddof=1) / abs(mean), mean, vals.var(ddof=1)


def measure_local_rho(rng, snr_db, tone_bin, trials=1500):
    """Measured mean periodogram at the tone bin minus the noise floor,
    over noise-only floor, giving the empirical local SNR rho."""
    bin_sel = int(round(tone_bin))
    frames_s = gen_frames(rng, trials, 1, NFFT, NFFT // 2, "white", snr_db, tone_bin)
    p_s = periodograms(frames_s)[:, 0, bin_sel].mean()
    frames_n = gen_frames(rng, trials, 1, NFFT, NFFT // 2, "white", -200.0, tone_bin)
    p_n = periodograms(frames_n)[:, 0, bin_sel].mean()
    return (p_s - p_n) / p_n, p_n


# ----------------------------------------------------------------------------
# Experiment A: the 1/K law across SNR
# ----------------------------------------------------------------------------
def experiment_a(rng):
    rows = []
    for snr in SNR_DB:
        rho, floor = measure_local_rho(rng, snr, tone_bin=8.0)
        for k in KS:
            tr = TRIALS_LAW if k <= 64 else TRIALS_FLOOR
            cv, mean, var = measure_cv(rng, "white", k, snr, 8.0, tr)
            pred = cv_law(k, rho, inflate=1.0)
            rows.append(dict(snr_db=snr, k=k, rho=rho, cv_meas=cv,
                             cv_pred=pred, q=cv / pred, mean=mean))
    return rows


# ----------------------------------------------------------------------------
# Experiment B: Hann overlap inflation factor (noise-only, tone off)
# ----------------------------------------------------------------------------
def experiment_b(rng):
    rows = []
    for hop_ratio in (1.0, 0.5, 0.25):
        hop = max(1, int(round(NFFT * hop_ratio)))
        for k in KS:
            tr = TRIALS_LAW if k <= 64 else TRIALS_FLOOR
            vals = np.empty(tr)
            for c0 in range(0, tr, CHUNK):
                m = min(CHUNK, tr - c0)
                frames = gen_frames(rng, m, k, NFFT, hop, "white", -200.0, 8.0)
                vals[c0:c0 + m] = avg_estimate(frames, k)[:, 16]
            cv = vals.std(ddof=1) / vals.mean()
            rows.append(dict(hop_ratio=hop_ratio, k=k, cv_meas=cv,
                             cv2=cv ** 2))
    # fit IF from CV^2 * K (white noise, C(0)=1)
    summary = {}
    for hop_ratio in (1.0, 0.5, 0.25):
        sub = [r for r in rows if r["hop_ratio"] == hop_ratio and r["k"] >= 4]
        ift = np.mean([r["cv2"] * r["k"] for r in sub])
        summary[hop_ratio] = ift
    return rows, summary


# ----------------------------------------------------------------------------
# Experiment C1: off-grid tone leakage floor, per band
# ----------------------------------------------------------------------------
def experiment_c1(rng):
    snr = 0.0
    delta = 0.5                      # worst-case half-bin offset
    tone_bin = 8.0 + delta
    k = 256
    trials = TRIALS_FLOOR
    # measured per-bin mean and CV of the K-frame average
    acc_p = np.zeros(NFFT // 2 + 1)
    acc_p2 = np.zeros(NFFT // 2 + 1)
    for c0 in range(0, trials, CHUNK):
        m = min(CHUNK, trials - c0)
        frames = gen_frames(rng, m, k, NFFT, NFFT // 2, "leak", snr, tone_bin)
        est = avg_estimate(frames, k)
        acc_p += est.sum(axis=0)
        acc_p2 += (est ** 2).sum(axis=0)
    mean = acc_p / trials
    var = acc_p2 / trials - mean ** 2
    cv = np.sqrt(np.clip(var, 0, None)) / mean
    # analytic leakage and analytic stochastic-law CV at each bin's local rho
    rho_leak = leakage_envelope(snr, delta, 8)
    rho_local = np.maximum(mean - 1.0, 0.0) / 1.0      # floor ~ 1 by scaling
    cv_pred = cv_law(k, rho_local, inflate=1.5)
    bands = []
    for b in range(NFFT // 2 + 1):
        d = abs(b - 8.5)
        if d <= 1.5:
            band = "main-lobe (|d|<=1.5 bin)"
        elif d <= 5.5:
            band = "near sidelobe (2..5 bins)"
        elif d <= 15.5:
            band = "mid sidelobe (6..15 bins)"
        else:
            band = "far sidelobe (>15 bins)"
        bias_frac = max(mean[b] - 1.0, 0.0) / mean[b]
        bands.append(dict(bin=b, band=band, mean=mean[b], cv_meas=cv[b],
                          cv_pred=cv_pred[b], rho_leak=rho_leak[b],
                          rho_local=rho_local[b], bias_frac=bias_frac,
                          excess=cv[b] / max(cv_pred[b], 1e-12)))
    return bands, dict(snr_db=snr, delta=delta, k=k)


# ----------------------------------------------------------------------------
# Experiment C2: coloured AR(1) noise, law departure vs phi and SNR
# ----------------------------------------------------------------------------
def experiment_c2(rng):
    rows = []
    for phi in (0.0, 0.5, 0.9):
        for snr in (-10.0, 0.0, 10.0):
            rho, _ = measure_local_rho(rng, snr, 8.0)
            for k in (16, 256):
                cv, _, _ = measure_cv(rng, "ar1", k, snr, 8.0,
                                      TRIALS_FLOOR if k > 64 else TRIALS_LAW,
                                      ar_phi=phi)
                pred = cv_law(k, rho, inflate=1.5)
                rows.append(dict(phi=phi, snr_db=snr, k=k, rho=rho,
                                 cv_meas=cv, cv_pred=pred, q=cv / pred))
    return rows


# ----------------------------------------------------------------------------
# Experiment C3: non-stationarity floor (noise envelope modulation)
# ----------------------------------------------------------------------------
def experiment_c3(rng):
    rows = []
    # one block = K frames * hop samples; modulation period chosen inside block
    for depth in (0.0, 0.25, 0.5):
        for k in KS:
            block_samples = (k - 1) * (NFFT // 2) + NFFT
            rate = 1.5 / block_samples          # ~1.5 modulation cycles per block
            cv, _, _ = measure_cv(rng, "nonstat", k, 0.0, 8.0,
                                  TRIALS_FLOOR if k > 64 else TRIALS_LAW,
                                  mod_rate=rate, mod_depth=depth)
            rho, _ = measure_local_rho(rng, 0.0, 8.0)
            floor = nonstat_floor_analytic(depth, k)
            pred = cv_law(k, rho, inflate=1.5)
            pred_with_floor = math.sqrt(pred ** 2 + floor ** 2)
            rows.append(dict(depth=depth, k=k, cv_meas=cv,
                             cv_pred=pred, floor_pred=floor,
                             cv_pred_with_floor=pred_with_floor,
                             q=cv / pred, q_floor=cv / pred_with_floor))
    return rows


# ----------------------------------------------------------------------------
# Experiment C4: magnitude (Rice) averaging vs power averaging
# ----------------------------------------------------------------------------
def experiment_c4(rng):
    rows = []
    for snr in SNR_DB:
        rho, _ = measure_local_rho(rng, snr, 8.0)
        for k in (16, 256):
            tr = TRIALS_FLOOR if k > 64 else TRIALS_LAW
            cv_p, mean_p, _ = measure_cv(rng, "white", k, snr, 8.0, tr, "power")
            cv_m, mean_m, _ = measure_cv(rng, "white", k, snr, 8.0, tr, "mag")
            rows.append(dict(snr_db=snr, k=k, rho=rho,
                             cv_power=cv_p, cv_mag=cv_m,
                             mag_bias=(mean_m / math.sqrt(max(mean_p, 1e-12)) - 1.0),
                             q_mag=cv_m / cv_p))
    return rows


# ----------------------------------------------------------------------------
# Experiment D: design table K(J_max, SNR) and measured residual jitter
# ----------------------------------------------------------------------------
def experiment_d(rng, inflate=1.5, j_targets=(0.20, 0.10, 0.05)):
    rows = []
    for jmax in j_targets:
        for snr in (-10.0, -3.0, 0.0, 6.0, 20.0):
            rho, _ = measure_local_rho(rng, snr, 8.0)
            k = design_frames(jmax, rho, inflate)
            cv, _, _ = measure_cv(rng, "white", k, snr, 8.0,
                                  TRIALS_FLOOR if k > 64 else TRIALS_LAW)
            rows.append(dict(j_max=jmax, snr_db=snr, rho=rho,
                             k_design=k, j_pred=cv_law(k, rho, inflate),
                             j_meas=cv, margin=jmax - cv))
    return rows


# ----------------------------------------------------------------------------
# Experiment B2: per-bin overlap correlation g(b) and its two powers
# ----------------------------------------------------------------------------
def hann_overlap_g(hop_ratio=0.5):
    """Bin-independent complex-DFT lag correlation magnitude for Hann frames.

    Frame 1 uses a *local* phase exp(-j2pi b m/N); in cov(X0[b], conj X1[b])
    the per-sample phase factors cancel, leaving the constant window-segment
    inner product times a unit-magnitude phase:
        g = (1/S2) sum_m w[m] w[m + hop],   independent of b.
    For 50% Hann overlap this sums to 1/6 analytically (and in float).
    """
    n = NFFT
    hop = int(round(n * hop_ratio))
    c = sum(W[m] * W[m + hop] for m in range(n - hop))
    return c / S2


def experiment_b2(rng):
    """Noise-only per-bin lag-1 periodogram correlation vs the flat g^2 law."""
    trials = 6000
    k = 2
    g = hann_overlap_g(0.5)
    frames = gen_frames(rng, trials, k, NFFT, NFFT // 2, "white", -200.0, 8.0)
    p = periodograms(frames)
    rows = []
    for b in range(NFFT // 2 + 1):
        p0, p1 = p[:, 0, b], p[:, 1, b]
        gamma = np.cov(p0, p1)[0, 1] / p0.var(ddof=1)
        # interior bins: periodogram lag correlation = g^2 = 1/36. DC/Nyquist
        # are real single-DOF bins with a different 4th-moment factor and are
        # marked separately rather than forced onto the interior law.
        rows.append(dict(bin=b, gamma_meas=gamma, gamma_pred=g ** 2, g=g,
                         edge=(b == 0 or b == NFFT // 2)))
    return rows


# ----------------------------------------------------------------------------
# Experiment A2: refined rho- and bin-dependent inflation law
# ----------------------------------------------------------------------------
def experiment_a2(rng):
    k = 32
    g = hann_overlap_g(0.5)             # bin-independent, 1/6
    b = 8
    rows = []
    for snr in (-15.0, -10.0, -3.0, 0.0, 6.0, 20.0):
        rho, _ = measure_local_rho(rng, snr, float(b))
        cv, _, _ = measure_cv(rng, "white", k, snr, float(b), TRIALS_FLOOR)
        if_const = 1.5
        if_refined = 1.0 + 2.0 * (1.0 - 1.0 / k) * (g ** 2 + 2.0 * rho * g) \
                              / (1.0 + 2.0 * rho)
        rows.append(dict(snr_db=snr, k=k, rho=rho, cv_meas=cv,
                         j_const15=cv_law(k, rho, if_const),
                         j_refined=cv_law(k, rho, if_refined),
                         if_refined=if_refined,
                         q_vs_const=cv / cv_law(k, rho, if_const),
                         q_vs_refined=cv / cv_law(k, rho, if_refined)))
    return rows


# ----------------------------------------------------------------------------
# Experiment C1b: leakage bias vs SNR (half-bin off-grid tone), main/sidelobe
# ----------------------------------------------------------------------------
def experiment_c1b(rng):
    delta = 0.5
    tone_pos = 8.0 + delta
    k = 128
    trials = 1200
    rows = []
    for snr in (-15.0, -10.0, -5.0, 0.0, 5.0, 10.0):
        frames = gen_frames(rng, trials, k, NFFT, NFFT // 2, "leak", snr, tone_pos)
        est = avg_estimate(frames, k)
        mean = est.mean(axis=0)
        rho_leak = leakage_envelope(snr, delta, 8)
        for b, tag in ((8, "main-lobe"), (10, "main-lobe-edge"),
                       (12, "near-sidelobe"), (16, "mid-sidelobe")):
            bias = max(mean[b] - 1.0, 0.0) / mean[b]
            rows.append(dict(snr_db=snr, bin=b, region=tag,
                             bias_meas=bias, bias_pred=rho_leak[b] / (1.0 + rho_leak[b]),
                             rho_leak=rho_leak[b]))
    return rows


# ----------------------------------------------------------------------------
# Experiment C2b: coloured AR(1) noise, per-bin inflation vs the g(b) law
# ----------------------------------------------------------------------------
def experiment_c2b(rng):
    k = 64
    trials = 1500
    g = hann_overlap_g(0.5)
    rows = []
    for phi in (0.0, 0.5, 0.9):
        frames = gen_frames(rng, trials, k, NFFT, NFFT // 2, "ar1", -200.0, 8.0,
                            ar_phi=phi)
        est = avg_estimate(frames, k)
        mu = est.mean(axis=0)
        sd = est.std(axis=0, ddof=1)
        cv = sd / mu
        for b in range(NFFT // 2 + 1):
            if_pred = 1.0 + 2.0 * (1.0 - 1.0 / k) * g ** 2
            rows.append(dict(phi=phi, bin=b, cv_meas=cv[b],
                             cv_pred_white=cv_law(k, 0.0, if_pred),
                             psd=mu[b], q=cv[b] / cv_law(k, 0.0, if_pred)))
    return rows


# ----------------------------------------------------------------------------
# Experiment C3b: analytic non-stationarity floor from the envelope
# ----------------------------------------------------------------------------
def nonstat_floor_analytic(depth, k, nfft=NFFT):
    """Between-frame CV of the mean power delivered by env(t)=1+d sin(wt).

    Each frame's expected noise power is the mean of env^2 over its samples;
    the floor jitter is the CV of those frame means across the K-frame block.
    """
    hop = nfft // 2
    total = (k - 1) * hop + nfft
    rate = 1.5 / total
    t = np.arange(total)
    env2 = (1.0 + depth * np.sin(2.0 * np.pi * rate * t)) ** 2
    frame_means = np.array([env2[f * hop:f * hop + nfft].mean() for f in range(k)])
    return frame_means.std() / frame_means.mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    print("[1/6] experiment A: 1/K law across SNR ...")
    a = experiment_a(rng)
    with open(os.path.join(args.out, "expA_law.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(a[0].keys()))
        wr.writeheader()
        wr.writerows(a)

    print("[2/6] experiment B: Hann overlap inflation ...")
    b, bsum = experiment_b(rng)
    with open(os.path.join(args.out, "expB_overlap.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(b[0].keys()))
        wr.writeheader()
        wr.writerows(b)

    print("[2b] experiment B2: per-bin overlap correlation g(b) ...")
    b2 = experiment_b2(rng)
    with open(os.path.join(args.out, "expB2_perbin_g.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(b2[0].keys()))
        wr.writeheader()
        wr.writerows(b2)

    print("[2c] experiment A2: refined rho-dependent inflation law ...")
    a2 = experiment_a2(rng)
    with open(os.path.join(args.out, "expA2_refined_law.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(a2[0].keys()))
        wr.writeheader()
        wr.writerows(a2)

    print("[3/6] experiment C1: leakage floor per band ...")
    c1, c1meta = experiment_c1(rng)
    with open(os.path.join(args.out, "expC1_leakage.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(c1[0].keys()))
        wr.writeheader()
        wr.writerows(c1)

    print("[3b] experiment C1b: leakage bias vs SNR ...")
    c1b = experiment_c1b(rng)
    with open(os.path.join(args.out, "expC1b_leakage_snr.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(c1b[0].keys()))
        wr.writeheader()
        wr.writerows(c1b)

    print("[4/6] experiment C2/C3/C4: coloured, non-stationary, magnitude ...")
    c2 = experiment_c2(rng)
    with open(os.path.join(args.out, "expC2_colored.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(c2[0].keys()))
        wr.writeheader()
        wr.writerows(c2)
    c2b = experiment_c2b(rng)
    with open(os.path.join(args.out, "expC2b_colored_perbin.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(c2b[0].keys()))
        wr.writeheader()
        wr.writerows(c2b)
    c3 = experiment_c3(rng)
    with open(os.path.join(args.out, "expC3_nonstat.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(c3[0].keys()))
        wr.writeheader()
        wr.writerows(c3)
    c4 = experiment_c4(rng)
    with open(os.path.join(args.out, "expC4_magnitude.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(c4[0].keys()))
        wr.writeheader()
        wr.writerows(c4)

    print("[5/6] experiment D: frame-count design table ...")
    d = experiment_d(rng)
    with open(os.path.join(args.out, "expD_design.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(d[0].keys()))
        wr.writeheader()
        wr.writerows(d)

    print("[6/6] writing summary ...")
    with open(os.path.join(args.out, "summary.txt"), "w") as f:
        f.write("Hann inflation factor (CV^2*K fit, C(0)=1):\n")
        for hop, ift in bsum.items():
            f.write(f"  hop={hop:>4}: IF={ift:.4f}\n")
        f.write("\nExp A worst |q-1| over all (SNR,K): "
                f"{max(abs(r['q'] - 1) for r in a):.4f}\n")
        f.write("Exp A2 refined-law worst |q-1|: "
                f"{max(abs(r['q_vs_refined'] - 1) for r in a2):.4f}"
                f"  (constant-IF-1.5 law: {max(abs(r['q_vs_const'] - 1) for r in a2):.4f})\n")
        f.write("Exp B2 per-bin gamma worst |meas-pred| (interior): "
                f"{max(abs(r['gamma_meas'] - r['gamma_pred']) for r in b2 if not r['edge']):.4f}\n")
        f.write("Exp C1 worst excess ratio (meas/pred): "
                f"{max(r['excess'] for r in c1):.3f}\n")
        f.write("Exp C1b worst |bias meas-pred|: "
                f"{max(abs(r['bias_meas'] - r['bias_pred']) for r in c1b):.4f}\n")
        f.write("Exp C2 worst q (coloured): "
                f"{max(r['q'] for r in c2):.3f}\n")
        f.write("Exp C2b per-bin colored q range: "
                f"[{min(r['q'] for r in c2b):.3f},{max(r['q'] for r in c2b):.3f}]\n")
        f.write("Exp C3 worst q (nonstat, K>=128): "
                f"{max(r['q'] for r in c3 if r['k'] >= 128):.3f}"
                f"  with-floor model: "
                f"{max(r['q_floor'] for r in c3 if r['k'] >= 128):.3f}\n")
        f.write("Exp C4 worst mag/power CV ratio: "
                f"{max(r['q_mag'] for r in c4):.3f}\n")
        f.write("Exp D all designs within budget: "
                f"{all(r['margin'] > 0 for r in d)}\n")
    print("done ->", args.out)


if __name__ == "__main__":
    main()
