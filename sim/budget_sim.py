#!/usr/bin/env python3
"""
budget_sim.py — Background-budget adaptive audio pipeline simulation.

Scenario
--------
The same music segment runs through a 5-stage processing chain (one stage per
C module in this repo: mixer, resampler, spatial, limiter, fft). When the app
moves to the background the processing budget shrinks, so each stage can run at
one of three precision levels (high / mid / low). Three budget tiers
(foreground / constrained / background) are simulated, one full pipeline run
per tier, and every run records the actual precision of each stage plus the
spectral difference of its output versus the full-precision reference.

Degradation order is decided by perceptual contribution estimated from spectral
masking: each stage is alone reduced to "low" while the others stay "high"; the
audible part of its error (error energy above the masking threshold of the
reference) is its masking contribution. The least-contributing stage is
downgraded first.

A second estimator with the masking sign flipped is also evaluated. At the
lowest budget tier, if that estimator downgrades a stage the correct estimator
protects (highest contribution), the inversion check fires and both
contribution values and both orders are listed side by side.

Run:  python -I budget_sim.py [--out <output_dir>]
Outputs JSON / CSV / WAV / PNG into the given directory (default
../artifacts/budget_sim relative to this script).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import wave
from dataclasses import dataclass, field

import numpy as np

# --------------------------------------------------------------------------
# Global audio constants
# --------------------------------------------------------------------------
SR = 48_000
DUR = 2.0
N = int(SR * DUR)
SEED = 20261124

STAGES = ["mixer", "resampler", "spatial", "limiter", "fft"]
STAGE_CN = {
    "mixer": "混音器",
    "resampler": "重采样器",
    "spatial": "空间声像",
    "limiter": "限幅器",
    "fft": "FFT 频谱处理",
}
LEVELS = ["high", "mid", "low"]
LEVEL_CN = {"high": "高精度", "mid": "中精度", "low": "低精度"}

# Compute cost per stage per precision level (abstract processing units).
COST = {
    "mixer":     {"high": 3,  "mid": 2, "low": 1},
    "resampler": {"high": 8,  "mid": 3, "low": 1},
    "spatial":   {"high": 4,  "mid": 2, "low": 1},
    "limiter":   {"high": 5,  "mid": 2, "low": 1},
    "fft":       {"high": 12, "mid": 5, "low": 2},
}
# Three background budget tiers. All-high baseline costs 32; minimum is 6.
BUDGETS = [("foreground", 32), ("constrained", 20), ("background", 10)]

# Analysis STFT used both for the masking model and the spectral-diff metric.
AN_N = 2048
AN_HOP = 512
AN_WIN = np.hanning(AN_N)

# Categorical palette (dataviz reference, light mode), fixed order per stage.
STAGE_COLOR = {
    "mixer": "#2a78d6",
    "resampler": "#eb6834",
    "spatial": "#1baf7a",
    "limiter": "#eda100",
    "fft": "#e87ba4",
}
LEVEL_COLOR = {"high": "#1c5cab", "mid": "#86b6ef", "low": "#cde2fb"}
LEVEL_EN = {"high": "high", "mid": "mid", "low": "low"}


# --------------------------------------------------------------------------
# Test signal: deterministic, music-like stereo segment, identical every run
# --------------------------------------------------------------------------
def make_signal() -> np.ndarray:
    rng = np.random.default_rng(SEED)
    t = np.arange(N) / SR

    def sine(f, a, ph=0.0):
        return a * np.sin(2 * np.pi * f * t + ph)

    # Bass: 110 Hz saw-like tone (harmonics 1..6, 1/k) — strong low masker.
    bass = sum(sine(110 * k, 0.30 / k) for k in range(1, 7))
    # Mid chord with a slow amplitude LFO.
    lfo = 0.75 + 0.25 * np.sin(2 * np.pi * 0.7 * t)
    chord = lfo * (sine(440, 0.10) + sine(554.37, 0.09) + sine(659.25, 0.08))
    # High partials, deliberately quiet (below the bass masker level so the
    # high-frequency region keeps a low masking threshold).
    high = sine(3520, 0.030) + sine(5280, 0.020) + sine(7040, 0.014) \
        + sine(9600, 0.010)
    # Pink-ish noise: white noise through a simple 1/f approximation.
    white = rng.standard_normal(N)
    pink = np.zeros(N)
    b = np.zeros(7)
    for i in range(N):
        w = white[i]
        b[0] = 0.99765 * b[0] + w * 0.0990460
        b[1] = 0.96300 * b[1] + w * 0.2965164
        b[2] = 0.57000 * b[2] + w * 1.0526913
        b[3] = -0.18480 * b[3] + w * 0.9604890
        b[4] = -0.46200 * b[4] + w * 0.6787290
        b[5] = -0.73800 * b[5] + w * 0.2490200
        b[6] = -0.89100 * b[6] + w * 0.0538120
        pink[i] = sum(b) * 0.11
    pink *= 0.03 / (np.max(np.abs(pink)) + 1e-12)
    # Periodic drum hits: each one crosses the limiter ceiling, so limiter
    # envelope timing (block-hold + coefficient quantization) leaves audible
    # gain-pump artifacts in the otherwise quiet tails between hits.
    beat = 0.0
    for b0 in np.arange(0.25, DUR - 0.05, 0.5):
        env_b = np.exp(-np.maximum(t - b0, 0.0) / 0.045) * (t >= b0)
        beat += 1.15 * env_b * (sine(180, 1.0) + 0.6 * sine(320, 1.0))
    # Single decaying transient near t = 1.0 s as an extra limiter trigger.
    env = np.exp(-np.maximum(t - 1.0, 0.0) / 0.030) * (t >= 1.0)
    transient = 0.95 * env * (sine(880, 1.0) + 0.5 * sine(1760, 1.0))
    # Global gentle fade envelope.
    shape = np.minimum(1.0, np.minimum(t / 0.05, (DUR - t) / 0.05))

    bus_a = (bass + pink + 0.6 * beat) * shape
    bus_b = (chord + high + transient + 0.4 * beat) * shape
    return np.stack([bus_a, bus_b], axis=1)


# --------------------------------------------------------------------------
# Stage 1: mixer — ramp-gain summing of two buses (mixer.c)
# Precision = gain-ramp quantization depth and arithmetic width.
# --------------------------------------------------------------------------
def stage_mixer(x: np.ndarray, level: str) -> np.ndarray:
    # low precision emulates a coarse 4-step hardware gain ladder plus f16-ish
    # accumulation. Both buses are quantized independently, so their rounding
    # noise does not cancel and spreads broadband across the spectrum.
    dt = {"high": np.float64, "mid": np.float32, "low": np.float16}[level]
    bits = {"high": 0, "mid": 8, "low": 2}[level]  # 0 = unquantized
    i = np.arange(N)
    # Ramp from 0.7 to 1.0 over the first 0.5 s, then hold (mirrors mix_ramp).
    g = 0.7 + 0.3 * np.minimum(i / (0.5 * SR), 1.0)
    if bits:
        q = 1 << bits
        g = np.round(g * q) / q
    g = g.astype(dt)
    a = x[:, 0].astype(dt)
    b = x[:, 1].astype(dt)
    y = (a * g + b * g).astype(np.float64)
    return np.stack([y, y], axis=1)


# --------------------------------------------------------------------------
# Stage 2: resampler — upsample x2 then downsample x2, length preserved
# (resampler.c). Precision = interpolation kernel.
# --------------------------------------------------------------------------
def _resample_chain(sig: np.ndarray, kernel: str) -> np.ndarray:
    def interp(s, positions):
        if kernel == "nearest":
            idx = np.clip(np.round(positions).astype(int), 0, len(s) - 1)
            return s[idx]
        if kernel == "linear":
            i0 = np.floor(positions).astype(int)
            i1 = np.minimum(i0 + 1, len(s) - 1)
            frac = positions - i0
            return (1 - frac) * s[i0] + frac * s[i1]
        # windowed sinc, 16 taps (8 each side)
        half = 8
        base = np.floor(positions).astype(int)
        frac = positions - base
        out = np.zeros_like(positions)
        for k in range(-half + 1, half + 1):
            n = k - frac
            h = np.sinc(n) * np.blackman(2 * half)[k + half - 1]
            idx = np.clip(base + k, 0, len(s) - 1)
            out += h * s[idx]
        return out

    up_pos = np.arange(2 * len(sig)) / 2.0
    up = interp(sig, up_pos)
    down_pos = np.arange(len(sig)) * 2.0 / 2.0  # back to original grid
    out = interp(up, down_pos)
    return out[: len(sig)]


def stage_resampler(x: np.ndarray, level: str) -> np.ndarray:
    kernel = {"high": "sinc16", "mid": "linear", "low": "nearest"}[level]
    out = np.empty_like(x)
    for ch in range(2):
        out[:, ch] = _resample_chain(x[:, ch], kernel)
    return out


# --------------------------------------------------------------------------
# Stage 3: spatial — constant-power pan + inverse-distance attenuation
# (spatial.c). Precision = trig lookup table size.
# --------------------------------------------------------------------------
def stage_spatial(x: np.ndarray, level: str) -> np.ndarray:
    lut_n = {"high": 0, "mid": 128, "low": 3}[level]  # 0 = direct trig, f64
    t = np.arange(N) / SR
    pan = 0.6 * np.sin(2 * np.pi * 0.15 * t)
    angle = (pan + 1.0) * 0.25 * np.pi
    if lut_n:
        idx = np.round((pan + 1.0) * 0.5 * (lut_n - 1)).astype(int)
        grid = np.linspace(-1.0, 1.0, lut_n)
        ang_grid = (grid + 1.0) * 0.25 * np.pi
        left_lut = np.cos(ang_grid).astype(np.float32)
        right_lut = np.sin(ang_grid).astype(np.float32)
        left, right = left_lut[idx], right_lut[idx]
    else:
        left, right = np.cos(angle), np.sin(angle)
    # attenuate(): dist 2.5, ref 2.0, rolloff 1.0 -> 2 / (2 + 0.5)
    att = 2.0 / (2.0 + 1.0 * (2.5 - 2.0))
    mono = x[:, 0] + x[:, 1]
    return np.stack([mono * left * att, mono * right * att], axis=1)


# --------------------------------------------------------------------------
# Stage 4: limiter — one-pole attack/release envelope follower (limiter.c)
# Precision = coefficient quantization and envelope update block size.
# --------------------------------------------------------------------------
def stage_limiter(x: np.ndarray, level: str) -> np.ndarray:
    block = {"high": 1, "mid": 16, "low": 64}[level]
    coef_bits = {"high": 0, "mid": 10, "low": 3}[level]
    dtype = {"high": np.float64, "mid": np.float32, "low": np.float32}[level]
    ceiling = 0.9
    ca = np.exp(-1.0 / (0.005 * SR))
    cr = np.exp(-1.0 / (0.050 * SR))
    if coef_bits:
        q = 1 << coef_bits
        ca = float(np.round(ca * q) / q)
        cr = float(np.round(cr * q) / q)

    out = np.empty_like(x)
    for ch in range(2):
        sig = x[:, ch].astype(dtype)
        env = dtype(1.0)
        held_target = 1.0
        y = np.empty(N, dtype=np.float64)
        for i in range(N):
            if i % block == 0:
                window = sig[i : i + block]
                peak = float(np.max(np.abs(window))) if len(window) else 0.0
                held_target = ceiling / peak if peak > ceiling else 1.0
            c = ca if held_target < env else cr
            env = c * env + (1.0 - c) * held_target
            v = float(sig[i]) * float(env)
            y[i] = np.clip(v, -ceiling, ceiling)
        out[:, ch] = y
    return out


# --------------------------------------------------------------------------
# Stage 5: fft — STFT domain gentle tone shaping + overlap-add (fft.c)
# Precision = frame size, arithmetic width and retained bandwidth.
# --------------------------------------------------------------------------
def _stft(sig, n, hop):
    frames = 1 + (len(sig) - n) // hop
    spec = np.empty((frames, n // 2 + 1), dtype=np.complex128)
    for k in range(frames):
        seg = sig[k * hop : k * hop + n] * (
            AN_WIN if n == AN_N else np.hanning(n))
        spec[k] = np.fft.rfft(seg, n)
    return spec


def _istft(spec, n, hop, length):
    win = np.hanning(n)
    out = np.zeros(length)
    wsum = np.zeros(length)
    for k in range(spec.shape[0]):
        seg = np.fft.irfft(spec[k], n) * win
        s = k * hop
        e = min(s + n, length)
        out[s:e] += seg[: e - s]
        wsum[s:e] += win[: e - s] ** 2
    return out / np.maximum(wsum, 1e-12)


def stage_fft(x: np.ndarray, level: str) -> np.ndarray:
    n = {"high": 1024, "mid": 512, "low": 128}[level]
    hop = n // 4
    band_limit = {"high": SR / 2, "mid": 10_000.0, "low": 4_000.0}[level]
    dtype = {"high": np.complex128, "mid": np.complex64, "low": np.complex64}[level]
    freqs = np.fft.rfftfreq(n, 1.0 / SR)
    # Gentle tone control: -0.2*(f/nyq)^2 shelf, applied in the frequency domain.
    shelf = 1.0 - 0.2 * (freqs / (SR / 2.0)) ** 2
    shelf[freqs > band_limit] = 0.0

    out = np.empty_like(x)
    for ch in range(2):
        spec = _stft(x[:, ch], n, hop).astype(dtype).astype(np.complex128)
        spec *= shelf
        out[:, ch] = _istft(spec, n, hop, N)
    return out


STAGE_FN = {
    "mixer": stage_mixer,
    "resampler": stage_resampler,
    "spatial": stage_spatial,
    "limiter": stage_limiter,
    "fft": stage_fft,
}


def pipeline(x: np.ndarray, levels: dict) -> np.ndarray:
    y = x
    for s in STAGES:
        y = STAGE_FN[s](y, levels[s])
    return y


# --------------------------------------------------------------------------
# Analysis helpers: STFT magnitudes, Bark scale, masking threshold
# --------------------------------------------------------------------------
def mag_frames(sig: np.ndarray):
    """Magnitude spectrum per analysis frame for one mono channel."""
    frames = 1 + (len(sig) - AN_N) // AN_HOP
    mag = np.empty((frames, AN_N // 2 + 1))
    for k in range(frames):
        seg = sig[k * AN_HOP : k * AN_HOP + AN_N] * AN_WIN
        mag[k] = np.abs(np.fft.rfft(seg))
    return mag


def bark(f: np.ndarray) -> np.ndarray:
    return 13.0 * np.arctan(0.00076 * f) + 3.5 * np.arctan((f / 7500.0) ** 2)


def ath(f: np.ndarray) -> np.ndarray:
    """Absolute threshold of hearing (dB SPL), Terhardt approximation.

    The fit is only meaningful roughly between 20 Hz and 16 kHz: below 20 Hz
    the kHz**-0.8 term diverges and above 16 kHz the kHz**4 term does. Bins
    outside that range are clamped to the threshold at the nearest edge.
    """
    lo, hi = 20.0, 16_000.0
    edge = np.clip(f, lo, hi)
    kHz = edge / 1000.0
    return (
        3.64 * np.power(kHz, -0.8)
        - 6.5 * np.exp(-0.6 * (kHz - 3.3) ** 2)
        + 1e-3 * kHz ** 4
    )


def masking_threshold_db(mag: np.ndarray, freqs: np.ndarray):
    """
    Per-frame masking threshold in dB SPL.
      * input is the reference signal's magnitude spectrum, normalized so its
        peak maps to 90 dB SPL,
      * Johnston spreading function across Bark,
      * maskers summed in the power domain with a fixed 9 dB masking offset,
      * floor at the absolute threshold of hearing.
    Returns threshold array with the same frame/bin shape as mag.
    """
    peak = max(float(mag.max()), 1e-12)
    sig_db = 90.0 + 20.0 * np.log10(mag / peak + 1e-12)
    z = bark(freqs)
    t_ath = ath(freqs)
    spread_power = np.zeros_like(sig_db)
    # Spreading matrix: entry (i, j) = spread of masker j onto bin i.
    dz = z[:, None] - z[None, :]
    spread_db = (
        15.81 + 7.5 * (dz + 0.474) - 17.5 * np.sqrt(1.0 + (dz + 0.474) ** 2)
    )
    for k in range(sig_db.shape[0]):
        masked_db = sig_db[k][None, :] + spread_db - 9.0  # 9 dB masking offset
        spread_power[k] = np.sum(10.0 ** (masked_db / 10.0), axis=1)
    thr_spread = 10.0 * np.log10(spread_power + 1e-12)
    return np.maximum(np.broadcast_to(t_ath, sig_db.shape), thr_spread)


def error_spectrum_db(ref: np.ndarray, deg: np.ndarray):
    """
    Frame-wise error magnitude on the SAME SPL scale as the masking threshold:
    mean absolute magnitude difference over the two channels, normalized so the
    reference peak maps to 90 dB SPL.
    """
    lin_rows = []
    for ch in range(2):
        mr = mag_frames(ref[:, ch])
        md = mag_frames(deg[:, ch])
        nf = min(mr.shape[0], md.shape[0])
        lin_rows.append(np.abs(mr[:nf] - md[:nf]))
    lin = np.mean(lin_rows, axis=0)
    ref_peak = max(float(np.mean(
        [mag_frames(ref[:, ch]) for ch in range(2)], axis=0).max()), 1e-12)
    return 90.0 + 20.0 * np.log10(lin / ref_peak + 1e-12)


def signal_spectrum_db(sig: np.ndarray):
    specs = [mag_frames(sig[:, ch]) for ch in range(2)]
    mag = np.mean(specs, axis=0)
    peak = max(float(mag.max()), 1e-12)
    return 90.0 + 20.0 * np.log10(mag / peak + 1e-12)


# --------------------------------------------------------------------------
# Contribution estimators
# --------------------------------------------------------------------------
def contribution_correct(err_db: np.ndarray, thr_db: np.ndarray) -> float:
    """Audible error: mean excess of error over the masking threshold (dB)."""
    return float(np.mean(np.maximum(err_db - thr_db, 0.0)))


def contribution_flipped(err_db: np.ndarray, thr_db: np.ndarray) -> float:
    """
    Masking sign flipped: the estimator credits error energy that IS hidden
    below the threshold (mean threshold headroom consumed by inaudible error),
    i.e. well-masked regions are treated as perceptually important.
    """
    return float(np.mean(np.maximum(thr_db - err_db, 0.0)))


def rms_diff_db(ref: np.ndarray, deg: np.ndarray) -> float:
    e = ref - deg
    r = float(np.sqrt(np.mean(ref ** 2)))
    d = float(np.sqrt(np.mean(e ** 2)))
    return 20.0 * np.log10(d / r + 1e-12)


def band_diff_curve(ref: np.ndarray, deg: np.ndarray, freqs: np.ndarray):
    """Mean absolute magnitude error in 8 log-spaced bands, returned in dB."""
    edges = np.logspace(np.log10(30.0), np.log10(SR / 2.0), 9)
    rows = []
    for ch in range(2):
        mr = mag_frames(ref[:, ch])
        md = mag_frames(deg[:, ch])
        nf = min(mr.shape[0], md.shape[0])
        rows.append(np.abs(mr[:nf] - md[:nf]))
    lin = np.mean(rows, axis=0).mean(axis=0)
    ref_peak = max(float(np.mean(
        [mag_frames(ref[:, ch]) for ch in range(2)], axis=0).max()), 1e-12)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (freqs >= lo) & (freqs < hi)
        out.append(90.0 + 20.0 * np.log10(
            float(np.mean(lin[sel])) / ref_peak + 1e-12))
    centers = np.sqrt(edges[:-1] * edges[1:])
    return centers, np.array(out)


# --------------------------------------------------------------------------
# Budget allocation: greedy downgrade in contribution order
# --------------------------------------------------------------------------
def allocate_budget(budget: int, order: list) -> dict:
    levels = {s: "high" for s in STAGES}
    rank = {s: i for i, s in enumerate(order)}
    while sum(COST[s][levels[s]] for s in STAGES) > budget:
        candidates = [s for s in STAGES if levels[s] != "low"]
        if not candidates:
            break
        # Lowest contribution first; on a tie take the move that saves most.
        s = min(
            candidates,
            key=lambda c: (
                rank[c],
                -COST[c][levels[c]]
                + COST[c][LEVELS[LEVELS.index(levels[c]) + 1]],
            ),
        )
        levels[s] = LEVELS[LEVELS.index(levels[s]) + 1]
    return levels


# --------------------------------------------------------------------------
# WAV / reporting helpers
# --------------------------------------------------------------------------
def write_wav(path: str, sig: np.ndarray):
    pcm = np.clip(sig, -1.0, 1.0)
    pcm = (pcm * 32767.0).astype("<i2")
    with wave.open(path, "wb") as w:
        w.setnchannels(2)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes(pcm.tobytes())


def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


# --------------------------------------------------------------------------
# Figures (matplotlib, Agg; fixed categorical palette, direct labels, grids)
# --------------------------------------------------------------------------
def make_figures(out_dir, ref, thr_db, sig_db, freqs, probes, results,
                 correct_order, flip_order, allocations, band_curves):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 10,
        "axes.edgecolor": "#c3c2b7",
        "axes.labelcolor": "#52514e",
        "axes.titlesize": 12,
        "xtick.color": "#898781",
        "ytick.color": "#898781",
        "grid.color": "#e1e0d9",
        "grid.linewidth": 0.8,
        "figure.facecolor": "white",
        "axes.facecolor": "#fcfcfb",
    })
    fmax = 16_000.0
    sel = freqs <= fmax

    # Fig 1: masking model + per-stage probe error spectra.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    frame = sig_db.shape[0] // 2
    ax1.plot(freqs[sel], sig_db[frame][sel], color="#0b0b0b", lw=1.6,
             label="signal spectrum (mid frame)")
    ax1.plot(freqs[sel], thr_db[frame][sel], color="#d03b3b", lw=1.6,
             ls="--", label="masking threshold")
    ax1.plot(freqs[sel], ath(freqs)[sel], color="#898781", lw=1.2, ls=":",
             label="absolute threshold of hearing")
    ax1.set_ylabel("dB SPL")
    ax1.set_title("Spectral masking model on the reference output")
    ax1.grid(True, alpha=0.6)
    ax1.legend(loc="upper right", frameon=False)
    for s in STAGES:
        ax2.plot(freqs[sel], probes[s].err_db.mean(axis=0)[sel],
                 color=STAGE_COLOR[s], lw=1.8, label=s)
    ax2.plot(freqs[sel], thr_db.mean(axis=0)[sel], color="#d03b3b", lw=1.6,
             ls="--", label="mean masking threshold")
    ax2.set_xlabel("frequency (Hz)")
    ax2.set_ylabel("error magnitude (dB SPL scale)")
    ax2.set_title("Single-stage low-precision error spectra vs. masking threshold")
    ax2.set_xscale("log")
    ax2.grid(True, alpha=0.6)
    ax2.legend(loc="upper right", frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig1_masking_and_errors.png"), dpi=150)
    plt.close(fig)

    # Fig 2: contribution ranking, correct vs flipped, with degrade order.
    # Log scale keeps the tiny-but-nonzero audible contributions visible.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 5), sharey=False)
    order_c = sorted(STAGES, key=lambda s: probes[s].c_correct)
    order_f = sorted(STAGES, key=lambda s: probes[s].c_flip)
    pos_c = {s: i for i, s in enumerate(order_c)}
    pos_f = {s: i for i, s in enumerate(order_f)}
    for s in STAGES:
        ax1.barh(pos_c[s], probes[s].c_correct, color=STAGE_COLOR[s],
                 height=0.62, edgecolor="white", linewidth=1.5)
        ax1.text(probes[s].c_correct * 1.12, pos_c[s],
                 f"degrade #{correct_order.index(s) + 1}  "
                 f"{probes[s].c_correct:.5f} dB  {s}",
                 va="center", ha="left", fontsize=9, color="#0b0b0b")
        vf = probes[s].c_flip
        ax2.barh(pos_f[s], vf, color=STAGE_COLOR[s], height=0.62,
                 edgecolor="white", linewidth=1.5)
        ax2.text(vf * 1.02, pos_f[s],
                 f"flipped #{flip_order.index(s) + 1}  {vf:.2f} dB  {s}",
                 va="center", ha="left", fontsize=9, color="#0b0b0b")
    ax1.set_xscale("log")
    ax1.set_xlim(min(probes[s].c_correct for s in STAGES) * 0.7,
                 max(probes[s].c_correct for s in STAGES) * 6)
    ax2.set_xlim(0, max(probes[s].c_flip for s in STAGES) * 1.55)
    for ax, title in (
        (ax1, "Correct masking model: audible error (higher = protect)"),
        (ax2, "Flipped masking model: masked error counted as contribution"),
    ):
        ax.set_yticks(range(5))
        ax.set_yticklabels([""] * 5)
        ax.set_title(title, fontsize=11)
        ax.grid(True, axis="x", alpha=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig2_contribution_ranking.png"), dpi=150)
    plt.close(fig)

    # Fig 3: precision actually used at each budget tier.
    fig, ax = plt.subplots(figsize=(9, 4.6))
    tier_names = [b[0] for b in BUDGETS]
    lvl_val = {"high": 2, "mid": 1, "low": 0}
    for i, s in enumerate(STAGES):
        for j, tname in enumerate(tier_names):
            lv = allocations[tname][s]
            ax.add_patch(plt.Rectangle((j, i), 0.86, 0.78,
                                       facecolor=LEVEL_COLOR[lv],
                                       edgecolor="white", linewidth=2))
            ax.text(j + 0.43, i + 0.39, LEVEL_EN[lv],
                    ha="center", va="center", fontsize=10,
                    color="#0b0b0b" if lv != "high" else "white")
    ax.set_xlim(-0.1, len(tier_names))
    ax.set_ylim(-0.1, len(STAGES))
    ax.set_xticks(np.arange(len(tier_names)) + 0.43)
    ax.set_xticklabels([f"{t}\nbudget {b}" for t, b in BUDGETS])
    ax.set_yticks(np.arange(len(STAGES)) + 0.39)
    ax.set_yticklabels(STAGES)
    ax.invert_yaxis()
    ax.set_title("Actual per-stage precision under each budget tier")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig3_tier_precision.png"), dpi=150)
    plt.close(fig)

    # Fig 4: spectral difference curves of the three tiers vs reference.
    # The all-high reference differs only at the numerical epsilon floor
    # (~-150 dB); curves are clipped at -80 dB for display so the two
    # degraded tiers stay legible, with the clip noted on the plot.
    fig, ax = plt.subplots(figsize=(10, 5.2))
    tier_colors = {"foreground": "#1c5cab", "constrained": "#eb6834",
                   "background": "#d03b3b"}
    display_floor = -80.0
    for tname in tier_names:
        c, d = band_curves[tname]
        ax.plot(c, np.maximum(d, display_floor), lw=2.0, marker="o", ms=6,
                color=tier_colors[tname],
                label=f"{tname} (budget {dict(BUDGETS)[tname]})")
    ax.axhline(display_floor, color="#898781", lw=1.0, ls=":")
    ax.text(30, display_floor + 1.5,
            "display floor -80 dB; foreground residual is ~-150 dB epsilon",
            fontsize=8.5, color="#898781")
    ax.set_ylim(-85, 75)
    ax.set_xscale("log")
    ax.set_xlabel("frequency (Hz), band centers")
    ax.set_ylabel("spectral magnitude difference (dB SPL scale)")
    ax.set_title("Output spectrum difference vs full-precision reference")
    ax.grid(True, alpha=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig4_spectral_diff_tiers.png"), dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
@dataclass
class Probe:
    err_db: object = None
    c_correct: float = 0.0
    c_flip: float = 0.0
    rms_db: float = 0.0


def main():
    ap = argparse.ArgumentParser()
    default_out = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "artifacts", "budget_sim")
    ap.add_argument("--out", default=default_out)
    args = ap.parse_args()
    out = args.out
    os.makedirs(os.path.join(out, "wav"), exist_ok=True)

    x = make_signal()
    freqs = np.fft.rfftfreq(AN_N, 1.0 / SR)

    # Full-precision reference: every stage at high.
    ref_levels = {s: "high" for s in STAGES}
    ref = pipeline(x, ref_levels)
    write_wav(os.path.join(out, "wav", "reference_all_high.wav"), ref)

    # Masking threshold of the reference (shared by every contribution calc).
    sig_db = signal_spectrum_db(ref)
    thr_db = masking_threshold_db(
        np.mean([mag_frames(ref[:, ch]) for ch in range(2)], axis=0), freqs)

    # Single-stage probes: only that stage drops to low precision.
    probes = {}
    for s in STAGES:
        lv = dict(ref_levels)
        lv[s] = "low"
        ys = pipeline(x, lv)
        err = error_spectrum_db(ref, ys)
        nf = min(err.shape[0], thr_db.shape[0])
        probes[s] = Probe(
            err_db=err,
            c_correct=contribution_correct(err[:nf], thr_db[:nf]),
            c_flip=contribution_flipped(err[:nf], thr_db[:nf]),
            rms_db=rms_diff_db(ref, ys),
        )
        write_wav(os.path.join(out, "wav", f"probe_only_{s}_low.wav"), ys)

    correct_order = sorted(STAGES, key=lambda s: probes[s].c_correct)
    flip_order = sorted(STAGES, key=lambda s: probes[s].c_flip)

    # Three budget tiers with the correct ordering.
    allocations, tier_outputs, band_curves, tier_rows = {}, {}, {}, []
    for tname, budget in BUDGETS:
        lv = allocate_budget(budget, correct_order)
        y = pipeline(x, lv)
        allocations[tname] = lv
        tier_outputs[tname] = y
        write_wav(os.path.join(out, "wav", f"tier_{tname}.wav"), y)
        centers, bands = band_diff_curve(ref, y, freqs)
        band_curves[tname] = (centers, bands)
        err = error_spectrum_db(ref, y)
        nf2 = min(err.shape[0], thr_db.shape[0])
        tier_rows.append({
            "tier": tname,
            "budget": budget,
            "cost": sum(COST[s][lv[s]] for s in STAGES),
            "levels": lv,
            "rms_diff_db": rms_diff_db(ref, y),
            "audible_error_db": contribution_correct(err[:nf2], thr_db[:nf2]),
            "band_centers": centers.tolist(),
            "band_diff_db": bands.tolist(),
        })

    # Same lowest budget, but allocated with the flipped masking order.
    flip_low_levels = allocate_budget(BUDGETS[-1][1], flip_order)
    flip_low_out = pipeline(x, flip_low_levels)
    write_wav(os.path.join(out, "wav", "tier_background_flipped_masking.wav"),
              flip_low_out)

    # ---- inversion check -------------------------------------------------
    most_protected = correct_order[-1]          # highest correct contribution
    flipped_downgraded = [
        s for s in STAGES if flip_low_levels[s] != "high"
    ]
    inversion = most_protected in flipped_downgraded

    # ---- CSV / JSON ------------------------------------------------------
    contrib_rows = []
    for s in STAGES:
        contrib_rows.append([
            s, STAGE_CN[s],
            f"{probes[s].c_correct:.6f}",
            correct_order.index(s) + 1,
            round(probes[s].c_flip, 4),
            flip_order.index(s) + 1,
            round(probes[s].rms_db, 3),
        ])
    write_csv(os.path.join(out, "contributions.csv"),
              ["stage", "stage_cn", "masking_contribution_db",
               "correct_degrade_order", "flipped_contribution_db",
               "flipped_degrade_order", "single_stage_rms_diff_db"],
              contrib_rows)

    prec_rows = []
    for tname, budget in BUDGETS:
        for s in STAGES:
            prec_rows.append([tname, budget, s, STAGE_CN[s],
                              allocations[tname][s],
                              LEVEL_CN[allocations[tname][s]],
                              COST[s][allocations[tname][s]]])
    write_csv(os.path.join(out, "tier_precision.csv"),
              ["tier", "budget", "stage", "stage_cn",
               "level", "level_cn", "cost_units"], prec_rows)

    spec_header = ["tier", "budget", "rms_diff_db", "audible_error_db"]
    centers = tier_rows[0]["band_centers"]
    spec_header += [f"band_{c:.0f}hz_db" for c in centers]
    spec_rows = []
    for r in tier_rows:
        spec_rows.append([r["tier"], r["budget"], round(r["rms_diff_db"], 3),
                          round(r["audible_error_db"], 4)]
                         + [round(v, 3) for v in r["band_diff_db"]])
    write_csv(os.path.join(out, "spectral_diff.csv"), spec_header, spec_rows)

    report = {
        "sample_rate": SR, "duration_s": DUR, "seed": SEED,
        "cost_table": COST, "budgets": dict(BUDGETS),
        "correct_degrade_order": correct_order,
        "flipped_degrade_order": flip_order,
        "contributions": {
            s: {
                "masking_contribution_db": probes[s].c_correct,
                "correct_order": correct_order.index(s) + 1,
                "flipped_contribution_db": probes[s].c_flip,
                "flipped_order": flip_order.index(s) + 1,
                "single_stage_rms_diff_db": probes[s].rms_db,
            } for s in STAGES
        },
        "tiers": [
            {k: v for k, v in r.items() if k != "levels"}
            | {"levels": r["levels"]} for r in tier_rows
        ],
        "flipped_masking_low_budget": {
            "levels": flip_low_levels,
            "downgraded_stages": flipped_downgraded,
            "rms_diff_db": rms_diff_db(ref, flip_low_out),
        },
        "inversion_check": {
            "most_protected_stage": most_protected,
            "flipped_model_downgraded_it": inversion,
            "verdict": "INVERTED_MASKING_DETECTED" if inversion else "ok",
        },
    }
    with open(os.path.join(out, "report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    make_figures(out, ref, thr_db, sig_db, freqs, probes, tier_rows,
                 correct_order, flip_order, allocations, band_curves)

    # ---- console report --------------------------------------------------
    line = "=" * 78
    print(line)
    print("自适应处理预算模拟：三档预算 × 五环节，掩蔽贡献决定降级顺序")
    print(line)
    print(f"{'环节':<12}{'掩蔽贡献(正确)':>16}{'正确降级序':>12}"
          f"{'翻转贡献':>12}{'翻转序':>10}{'单环节RMS差dB':>16}")
    for s in STAGES:
        print(f"{s:<12}{probes[s].c_correct:>16.5f}"
              f"{correct_order.index(s) + 1:>12d}{probes[s].c_flip:>12.3f}"
              f"{flip_order.index(s) + 1:>10d}{probes[s].rms_db:>16.2f}")
    print()
    print("正确模型降级顺序（先降 -> 最后降）:",
          " -> ".join(f"{s}({STAGE_CN[s]})" for s in correct_order))
    print("翻转模型降级顺序（先降 -> 最后降）:",
          " -> ".join(f"{s}({STAGE_CN[s]})" for s in flip_order))
    print()
    for tname, budget in BUDGETS:
        lv = allocations[tname]
        r = next(r for r in tier_rows if r["tier"] == tname)
        prec = " | ".join(f"{s}:{LEVEL_CN[lv[s]]}" for s in STAGES)
        print(f"[{tname} budget={budget} 实际成本={r['cost']}]")
        print(f"  {prec}")
        print(f"  输出相对参考: RMS谱差={r['rms_diff_db']:.2f} dB, "
              f"可闻误差={r['audible_error_db']:.3f} dB")
    print()
    print("[最低预算 + 翻转掩蔽排序]")
    print("  精度:", " | ".join(
        f"{s}:{LEVEL_CN[flip_low_levels[s]]}" for s in STAGES))
    print(f"  被降级环节: {', '.join(flipped_downgraded)}")
    print()
    if inversion:
        print(f"!! 检测到掩蔽效应算反: 最不该降级的环节 "
              f"'{most_protected}({STAGE_CN[most_protected]})' "
              f"在翻转模型的最低档被降级。")
        print("   上表已将各环节掩蔽贡献值与两套降级顺序并列输出；")
        print("   详见 contributions.csv / fig2_contribution_ranking.png")
    else:
        print("翻转模型未在最低档降级最高贡献环节，未触发反转判定。")
    print(line)
    print(f"输出目录: {out}")


if __name__ == "__main__":
    main()
