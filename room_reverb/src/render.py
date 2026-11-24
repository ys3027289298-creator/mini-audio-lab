"""Render band-limited impulse responses and measure their decay.

Pipeline per receiver and per octave band:

  1. Early/specular field: every image source becomes a unit-energy raised
     cosine tone pulse, scaled so its energy equals the image's band energy
     (reflection gain plus 1/d^2 geometric spreading).
  2. Late field: image sources of finite order stop at a crossover time.  A
     band-limited Gaussian noise tail takes over, with an exponential envelope
     whose level and decay rate are calibrated to the analytic image-set EDC.
  3. T30 measurement: Schroeder backward integration of the rendered IR,
     linear fit between the -5 dB and -35 dB levels, extrapolated to T60.

The analytic EDC from images.py is run through the same T30 fitter so rendered
and exact specular estimates can be cross-checked.
"""

import numpy as np

from scipy.signal import butter, sosfiltfilt

from images import (analytic_edc, band_energy_envelope,
                    receiver_arrivals)
from room import BAND_FREQS

FS = 16000
DURATION = 1.6


def _tone_pulse(fc, fs):
    """Unit-energy raised-cosine tone pulse, two cycles long."""
    n_periods = 2
    n_samp = max(8, int(round(n_periods * fs / fc)))
    t = np.arange(n_samp) / fs
    p = np.sin(2 * np.pi * fc * t) * 0.5 * (1 - np.cos(2 * np.pi * t / (n_samp / fs)))
    return p / np.sqrt(np.sum(p ** 2))


def _bandpass_noise(n, fc, fs, rng):
    """White Gaussian noise band-limited to roughly the octave around fc."""
    low = max(20.0, fc / np.sqrt(2))
    high = min(fs / 2 * 0.95, fc * np.sqrt(2))
    sos = butter(4, [low / (fs / 2), high / (fs / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, rng.standard_normal(n))


def _late_tail(images, receiver, band_idx, t0, fs, duration, rng):
    """Exponentially decaying band-limited late field, energy-calibrated."""
    n = int(round((duration - t0) * fs))
    if n <= 1:
        return np.zeros(0)

    # Decay rate from the analytic EDC of the specular image set.
    t_a, edc_a = analytic_edc(images, receiver)
    slope_t60, _ = _fit_decay(t_a, edc_a[:, band_idx])
    if not np.isfinite(slope_t60) or slope_t60 <= 0:
        slope_t60 = 1.0  # fallback: 1 s, tail will still be negligible-ish

    # Absolute energy level at the crossover (linear scale).
    energy = band_energy_envelope(images, receiver)[:, band_idx]
    total_e = energy.sum()
    edc_at_t0 = np.interp(t0, t_a, edc_a[:, band_idx])
    target_energy = total_e * 10.0 ** (edc_at_t0 / 10.0)

    g = _bandpass_noise(n, BAND_FREQS[band_idx], fs, rng)
    t = np.arange(n) / fs
    env = 10.0 ** (-3.0 * t / slope_t60)  # energy envelope, 60 dB per T60
    tail = g * np.sqrt(env)
    # Calibrate total late energy to the analytic EDC level.
    tail *= np.sqrt(target_energy / max(np.sum(tail ** 2), 1e-30))
    return tail


def render_band_ir(images, receiver, band_idx, fs=FS, duration=DURATION,
                   add_late=True, seed=0):
    """Render one octave-band impulse response."""
    rng = np.random.default_rng(seed + band_idx)
    delays, _, spreading = receiver_arrivals(images, receiver)
    n_samp = int(round(duration * fs))
    ir = np.zeros(n_samp)
    pulse = _tone_pulse(BAND_FREQS[band_idx], fs)

    amp = images["gain"][:, band_idx] * spreading
    idx = np.round(delays * fs).astype(int)
    valid = (idx >= 0) & (idx + len(pulse) <= n_samp)
    for i in np.where(valid)[0]:
        j = idx[i]
        ir[j:j + len(pulse)] += amp[i] * pulse

    if add_late:
        t_cross = delays[valid].max()
        n_cross = int(round(t_cross * fs))
        tail = _late_tail(images, receiver, band_idx, t_cross,
                          fs, duration, rng)
        ir[n_cross:n_cross + len(tail)] += tail
    return ir


def render_ir(images, receiver, fs=FS, duration=DURATION, seed=0):
    """All band IRs plus a broadband sum used for waveform plots."""
    bands = [render_band_ir(images, receiver, b, fs, duration, seed=seed)
             for b in range(len(BAND_FREQS))]
    broadband = np.sum(bands, axis=0)
    return np.array(bands), broadband


def schroeder_edc(ir):
    """Schroeder backward-integrated decay curve in dB."""
    power = ir ** 2
    tail_energy = np.sum(power) - np.concatenate(([0.0], np.cumsum(power)[:-1]))
    total = max(np.sum(power), 1e-30)
    return 10.0 * np.log10(np.maximum(tail_energy, 1e-20) / total)


def _fit_decay(t, edc_db, upper=-5.0, lower=-35.0):
    """Linear fit between two dB levels; returns (T60, fit quality dict)."""
    finite = np.isfinite(edc_db)
    t, edc_db = t[finite], edc_db[finite]
    if len(t) < 10:
        return np.nan, {"r2": np.nan, "range_db": 0.0}

    def cross(level):
        below = np.where(edc_db <= level)[0]
        return t[below[0]] if len(below) else np.nan

    t_upper, t_lower = cross(upper), cross(lower)
    if not (np.isfinite(t_upper) and np.isfinite(t_lower)):
        return np.nan, {"r2": np.nan,
                        "range_db": float(edc_db.min() - edc_db.max())}

    mask = (t >= t_upper) & (t <= t_lower)
    slope, intercept = np.polyfit(t[mask], edc_db[mask], 1)
    pred = slope * t[mask] + intercept
    ss_res = np.sum((edc_db[mask] - pred) ** 2)
    ss_tot = np.sum((edc_db[mask] - edc_db[mask].mean()) ** 2)
    r2 = 1 - ss_res / max(ss_tot, 1e-30)
    t60 = -60.0 / slope
    return float(t60), {"r2": float(r2),
                        "range_db": float(upper - lower)}


def measure_t30(ir, fs=FS):
    """Estimate T60 from a rendered IR via Schroeder integration (T30 fit)."""
    t = np.arange(len(ir)) / fs
    return _fit_decay(t, schroeder_edc(ir))


def measure_edt(ir, fs=FS):
    """Early Decay Time: fit the 0 to -10 dB segment, extrapolated to T60.

    EDT reacts to the local early field (source/receiver geometry and nearby
    materials) much more strongly than T30, which is a late, diffuse-field
    quantity.
    """
    t = np.arange(len(ir)) / fs
    return _fit_decay(t, schroeder_edc(ir), upper=0.0, lower=-10.0)


def analytic_t30(images, receiver):
    """T30 fit applied to the exact specular EDC, per band."""
    t, edc = analytic_edc(images, receiver)
    out = []
    for b in range(edc.shape[1]):
        t60, quality = _fit_decay(t, edc[:, b])
        out.append((t60, quality))
    return out
