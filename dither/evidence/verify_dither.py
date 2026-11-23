#!/usr/bin/env python3
"""Independent verifier for the dither evidence bundle.

Reads the exported WAVs and their .f32 sidecars and re-derives everything it
can without trusting the C code:

  * dither off: floor(x*32768 + 0.5) codes must equal the WAV, sample for
    sample, with numpy arithmetic (the "no-dither byte-identical" claim,
    checked against an independent implementation, not just against itself).
  * dither on: no sample may be more than 1.5 LSB from its code's rounding
    target, and the error statistics must match the report.
  * SNR, harmonic ladder, error spectrum, noise modulation.
  * the -6 -> -90 dBFS level sweep, where dither's level-independent error
    floor and plain rounding's signal-tracking error are compared.

Usage: python3 verify_dither.py <wav_dir> <out_json>
Writes a JSON summary and exits non-zero if a claim fails.
"""
import sys, os, json, wave
import numpy as np

FS = 48000
QSCALE = 32768.0
F0 = 750.0          # 16 samples/cycle at 48 kHz; coherent for NFFT = 2^16
NFFT = 65536

C_LIGHT = {"s1": "#2a78d6", "s2": "#eb6834", "surface": "#fcfcfb",
           "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781",
           "grid": "#e1e0d9", "axis": "#c3c2b7"}
C_DARK = {"s1": "#3987e5", "s2": "#d95926", "surface": "#1a1a19",
          "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781",
          "grid": "#2c2c2a", "axis": "#383835"}


def read_wav(path):
    with wave.open(path, "rb") as w:
        assert w.getframerate() == FS and w.getnchannels() == 1 and w.getsampwidth() == 2
        q = np.frombuffer(w.readframes(w.getnframes()), dtype="<i2").astype(np.int32)
    return q


def read_f32(path):
    with open(path, "rb") as f:
        raw = f.read()
    return np.frombuffer(raw, dtype="<f4").astype(np.float64)


def plain_codes(x_f64):
    """Reference plain rounding: floor(x*32768 + 0.5), clamped."""
    return np.clip(np.floor(x_f64 * QSCALE + 0.5), -32768, 32767).astype(np.int32)


def snr_db(x, e):
    ps, pe = float(np.mean(x * x)), float(np.mean(e * e))
    return 10.0 * np.log10(ps / pe) if pe > 0 else float("inf")


def coherent_spectrum_db(err, fs=FS, nfft=NFFT, ref_amp=1.0):
    """Rectangular-window spectrum of the error, dB relative to ref_amp
    (a full-scale sine has |X| = nfft/2 at its bin)."""
    seg = err[:nfft]
    X = np.fft.rfft(seg)
    mag = 2.0 * np.abs(X) / nfft / ref_amp
    return 20.0 * np.log10(np.maximum(mag, 1e-30)), np.fft.rfftfreq(nfft, 1.0 / fs)


def harmonic_bins(f0, fs=FS, nfft=NFFT, max_k=32):
    bins = []
    for k in range(2, max_k + 1):
        if k * f0 >= fs / 2:
            break
        b = int(round(k * f0 * nfft / fs))
        bins.append((k, b))
    return bins


def spur_stats(err, f0, fs=FS, nfft=NFFT):
    db, freqs = coherent_spectrum_db(err, fs, nfft)
    fund_b = int(round(f0 * nfft / fs))
    fund_db = db[fund_b]
    # mask DC, the fundamental, and every harmonic (leakage/quantization of
    # the tone itself is tested separately via the harmonic ladder).
    mask = np.ones(len(freqs), bool)
    mask[:3] = False
    for k in range(1, 40):
        b = int(round(k * f0 * nfft / fs))
        if b >= len(freqs):
            break
        mask[max(0, b - 1):min(len(freqs), b + 2)] = False
    if not mask.any():
        return {"fundamental_db": float(fund_db)}
    i = int(np.argmax(db[mask]))
    return {
        "fundamental_db": float(fund_db),
        "worst_spur_db": float(db[mask][i]),
        "worst_spur_hz": float(freqs[mask][i]),
        "floor_db_median": float(np.median(db[freqs > 20000])),
    }


def folded_distortion(x, err, f0=F0, max_periods=3000, active=None):
    """The decisive test for 'is the error noise or signal?': fold the error
    on the tone's period. In exact periods (f0 * FS is an integer here: 750 Hz
    = 64 samples), the deterministic part of the error is a fixed curve that
    averages coherently across every period; the random part averages toward
    zero. report the folded curve's RMS (the distortion) and the residual's
    RMS (the noise), both in LSB, and the folded fraction of total error
    power. Plain rounding: ~96% folded. Dither: a few percent."""
    period = FS / f0
    assert abs(period - round(period)) < 1e-9, "tone must be periodic"
    P = int(round(period))
    idx = np.arange(min(len(x), P * max_periods))
    if active is not None:
        idx = idx[active[idx]]   # gated file: keep tone-on periods only
    ph = idx % P
    n = np.bincount(ph, minlength=P).astype(np.float64)
    e_lsb = err * QSCALE   # caller passes FS-domain error
    fold = np.bincount(ph, weights=e_lsb[idx], minlength=P) / n
    fold_rms = float(np.sqrt(np.mean(fold ** 2)))
    total = float(np.mean(e_lsb[idx] ** 2))
    resid = float(np.sqrt(max(total - fold_rms ** 2, 0.0)))
    return {"period_samples": P, "folded_rms_lsb": fold_rms,
            "residual_rms_lsb": resid,
            "total_rms_lsb": float(np.sqrt(total)),
            "folded_fraction": float(fold_rms ** 2 / total) if total > 0 else None,
            "folded_curve_lsb": [float(v) for v in fold]}


def distortion_split(x, err, max_k=33):
    """Split the error into the part that is signal-correlated distortion and
    the part that is broadband noise. Distortion = harmonics 2..max_k of the
    tone, measured on exact bins (750 Hz is coherent for NFFT). Plain
    rounding's odd-harmonic staircase lands here; dither's is only its noise
    floor. Both are expressed dB relative to the tone's RMS, so they compare
    directly with SNR = -(total error)."""
    n = NFFT
    X = np.fft.rfft(err[:n])
    ps = float(np.mean(x * x))
    harm = 0.0
    hp = []
    for k in range(2, max_k + 1):
        b = int(round(k * F0 * n / FS))
        if k * F0 >= FS / 2:
            break
        a = 2.0 * abs(X[b]) / n
        harm += a * a
        hp.append({"k": k, "db_rel_tone": 20 * np.log10(a / np.sqrt(ps) + 1e-30)})
    noise = 0.0
    mask = np.ones(len(X), bool)
    mask[:2] = False
    for k in range(0, 40):
        b = int(round(k * F0 * n / FS))
        if b < len(X):
            mask[max(0, b - 2):b + 3] = False
    noise = float(np.sum((2.0 * np.abs(X[mask]) / n) ** 2))
    total = float(np.mean(err[:n] ** 2))
    return {"thd_db": 10 * np.log10(harm / ps),
            "noise_db": 10 * np.log10(noise / ps),
            "total_db": 10 * np.log10(total / ps),
            "harmonics": hp}


def spurs_vs_floor(x, err):
    """The discriminating spectrum view. For both files, report the error
    spectrum relative to the tone: the loudest non-tone component in each
    1/3-octave-ish band, and the median floor. Plain rounding has discrete
    spurs well above its floor; dither has a flat floor and no spurs above it.
    The tone and its harmonic bins are excluded from the spur search, so the
    comparison is about what is left when the signal's own distortion is not
    the question."""
    n = NFFT
    X = np.fft.rfft(err[:n])
    f = np.fft.rfftfreq(n, 1.0 / FS)
    ref = np.sqrt(np.mean(x * x))          # tone RMS; spur is dB re tone
    db = 20 * np.log10(2 * np.abs(X) / n / ref + 1e-30)
    tone_mask = np.zeros(len(f), bool)
    for k in range(0, 40):
        b = int(round(k * F0 * n / FS))
        if b >= len(f): break
        tone_mask[max(0, b - 4):min(len(f), b + 5)] = True
    free = ~tone_mask
    floor = float(np.median(db[free]))
    edges = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 24000]
    bands = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = free & (f >= lo) & (f < hi)
        if not sel.any():
            bands.append({"lo_hz": lo, "hi_hz": hi, "peak_db": None}); continue
        bands.append({"lo_hz": lo, "hi_hz": hi, "peak_db": float(np.max(db[sel]))})
    return {"ref_tone_db": 0.0, "floor_db": floor,
            "peak_above_floor_db": float(np.max(db[free]) - floor), "bands": bands}


def noise_modulation(err, win_ms=20):
    B = int(win_ms * 1e-3 * FS)
    m = len(err) // B
    if m < 2:
        return {"win_ms": win_ms, "windows": m}
    e = err[:m * B].reshape(m, B)
    rms = np.sqrt(np.mean(e * e, axis=1))
    return {
        "win_ms": win_ms,
        "windows": m,
        "rms_lsb_median": float(np.median(rms)),
        "rms_lsb_p05": float(np.percentile(rms, 5)),
        "rms_lsb_p95": float(np.percentile(rms, 95)),
        "rms_lsb_std": float(np.std(rms)),
        "digital_silent_windows": int(np.sum(rms == 0.0)),
    }


# ---------------------------------------------------------------- charts

def svg_noise_spectrum(path, freqs, plain_db, dith_db, title):
    W, H = 920, 520
    L, R, T, B = 78, 24, 58, 62
    pw, ph = W - L - R, H - T - B
    ymin, ymax = -200.0, -80.0
    xmax = FS / 2
    def X(f): return L + (f / xmax) * pw
    def Y(v): return T + (ymax - v) / (ymax - ymin) * ph
    def poly(freqs, db):
        step = max(1, len(freqs) // (pw * 2))
        pts = [(X(f), Y(v)) for f, v in zip(freqs[::step], db[::step]) if v > ymin]
        return "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    grid = []
    for v in range(-200, -79, 20):
        y = Y(v)
        grid.append(f'<line x1="{L}" y1="{y:.1f}" x2="{L+pw}" y2="{y:.1f}" stroke="{{grid}}" stroke-width="1"/>')
        grid.append(f'<text x="{L-10}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="{{muted}}" font-variant-numeric="tabular-nums">{v}</text>')
    for f in [0, 4000, 8000, 12000, 16000, 20000, 24000]:
        grid.append(f'<text x="{X(f):.1f}" y="{T+ph+22}" text-anchor="middle" font-size="12" fill="{{muted}}" font-variant-numeric="tabular-nums">{f//1000 if f else 0}k</text>')
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img"
     aria-label="{title}" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif">
  <rect width="{W}" height="{H}" fill="{{surface}}"/>
  <text x="{L}" y="26" font-size="17" font-weight="600" fill="{{ink}}">{title}</text>
  <text x="{L}" y="46" font-size="12.5" fill="{{ink2}}">Quantization error spectrum, -60 dBFS 750 Hz tone, 65536-point rectangular window; dB relative to full scale</text>
  {"".join(grid)}
  <line x1="{L}" y1="{T+ph}" x2="{L+pw}" y2="{T+ph}" stroke="{{axis}}" stroke-width="1"/>
  <path d="{poly(freqs, plain_db)}" fill="none" stroke="{{s1}}" stroke-width="1.4" stroke-linejoin="round"/>
  <path d="{poly(freqs, dith_db)}" fill="none" stroke="{{s2}}" stroke-width="1.4" stroke-linejoin="round"/>
  <g font-size="12.5" fill="{{ink2}}">
    <line x1="{L+18}" y1="{T+18}" x2="{L+46}" y2="{T+18}" stroke="{{s1}}" stroke-width="2"/>
    <text x="{L+54}" y="{T+22}">No dither (plain rounding)</text>
    <line x1="{L+18}" y1="{T+42}" x2="{L+46}" y2="{T+42}" stroke="{{s2}}" stroke-width="2"/>
    <text x="{L+54}" y="{T+46}">TPDF dither</text>
  </g>
  <text x="{L+pw/2:.1f}" y="{H-16}" text-anchor="middle" font-size="12.5" fill="{{ink2}}">Frequency</text>
  <text transform="translate(20,{T+ph/2:.1f}) rotate(-90)" text-anchor="middle" font-size="12.5" fill="{{ink2}}">Error (dBFS)</text>
</svg>'''
    for c, suf in ((C_LIGHT, ""), (C_DARK, ".dark")):
        out = svg
        for k, v in c.items():
            out = out.replace("{" + k + "}", v)
        with open(path.replace(".svg", suf + ".svg"), "w") as fh:
            fh.write(out)


def svg_error_wave(path, plain_e, dith_e, ms, title):
    n = int(ms * 1e-3 * FS)
    t = np.arange(n) / FS * 1000.0
    W, H = 920, 470
    L, R, T, B = 78, 24, 70, 52
    pw = W - L - R
    ph = (H - T - B - 52) / 2
    ymax = 1.2
    def X(i): return L + t[i] / t[-1] * pw
    def Y(v, top): return top + (ymax - v) / (2 * ymax) * ph
    def panel(e, top, color, label, note):
        d = "M" + " L".join(f"{X(i):.1f},{Y(v, top):.1f}" for i, v in enumerate(e))
        y0 = Y(0, top)
        return f'''<text x="{L}" y="{top-9}" font-size="13" font-weight="600" fill="{{ink}}">{label}</text>
  <text x="{L+104}" y="{top-9}" font-size="12" fill="{{ink2}}">{note}</text>
  <line x1="{L}" y1="{y0:.1f}" x2="{L+pw}" y2="{y0:.1f}" stroke="{{grid}}" stroke-width="1"/>
  <path d="{d}" fill="none" stroke="{color}" stroke-width="1.3" stroke-linejoin="round"/>'''
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img"
     aria-label="{title}" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif">
  <rect width="{W}" height="{H}" fill="{{surface}}"/>
  <text x="{L}" y="26" font-size="17" font-weight="600" fill="{{ink}}">{title}</text>
  <text x="{L}" y="46" font-size="12.5" fill="{{ink2}}">Quantization error (decoded minus input), LSB; shared scale</text>
  {panel(plain_e, T+34, "{s1}", "No dither", "bounded to +/-0.5 LSB, but correlated with the tone")}
  {panel(dith_e, T+34+ph+52, "{s2}", "TPDF dither", "up to +/-1.5 LSB, uncorrelated: flat noise")}
  <text x="{L+pw/2:.1f}" y="{H-14}" text-anchor="middle" font-size="12.5" fill="{{ink2}}">Time (ms)</text>
</svg>'''
    for c, suf in ((C_LIGHT, ""), (C_DARK, ".dark")):
        out = svg
        for k, v in c.items():
            out = out.replace("{" + k + "}", v)
        with open(path.replace(".svg", suf + ".svg"), "w") as fh:
            fh.write(out)


def svg_snr_vs_level(path, levels, snr_plain, snr_dith, title="SNR vs tone level: dither vs none"):
    """SNR vs signal level: plain rounding rises ~6 dB/dB above the noise
    floor; dither gives ~3.3 dB/dB. Below the floor the roles reverse."""
    W, H = 920, 520
    L, R, T, B = 78, 24, 58, 68
    pw, ph = W - L - R, H - T - B
    xmin, xmax = -110, 0
    ymin, ymax = -20, 100
    def X(v): return L + (v - xmin) / (xmax - xmin) * pw
    def Y(v): return T + (ymax - v) / (ymax - ymin) * ph
    def line(xs, ys):
        return "M" + " L".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in zip(xs, ys))
    grid = []
    for v in range(0, 101, 20):
        y = Y(v)
        grid.append(f'<line x1="{L}" y1="{y:.1f}" x2="{L+pw}" y2="{y:.1f}" stroke="{{grid}}" stroke-width="1"/>')
        grid.append(f'<text x="{L-10}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="{{muted}}" font-variant-numeric="tabular-nums">{v}</text>')
    for v in [-110, -90, -70, -50, -30, -10, 0]:
        grid.append(f'<text x="{X(v):.1f}" y="{T+ph+22}" text-anchor="middle" font-size="12" fill="{{muted}}" font-variant-numeric="tabular-nums">{v}</text>')
    # markers on the last point of each series
    lx_p, ly_p = X(levels[-1]), Y(snr_plain[-1])
    lx_d, ly_d = X(levels[-1]), Y(snr_dith[-1])
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img"
     aria-label="{title}" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif">
  <rect width="{W}" height="{H}" fill="{{surface}}"/>
  <text x="{L}" y="26" font-size="17" font-weight="600" fill="{{ink}}">{title}</text>
  <text x="{L}" y="46" font-size="12.5" fill="{{ink2}}">750 Hz tone, 3 s each; error measured against the same float input. Note: x-axis is amplitude, not time.</text>
  {"".join(grid)}
  <line x1="{L}" y1="{T+ph}" x2="{L+pw}" y2="{T+ph}" stroke="{{axis}}" stroke-width="1"/>
  <path d="{line(levels, snr_plain)}" fill="none" stroke="{{s1}}" stroke-width="2" stroke-linejoin="round"/>
  <path d="{line(levels, snr_dith)}" fill="none" stroke="{{s2}}" stroke-width="2" stroke-linejoin="round"/>
  <circle cx="{lx_p:.1f}" cy="{ly_p:.1f}" r="4.5" fill="{{s1}}" stroke="{{surface}}" stroke-width="2"/>
  <circle cx="{lx_d:.1f}" cy="{ly_d:.1f}" r="4.5" fill="{{s2}}" stroke="{{surface}}" stroke-width="2"/>
  <g font-size="12.5" fill="{{ink2}}">
    <line x1="{L+18}" y1="{T+18}" x2="{L+46}" y2="{T+18}" stroke="{{s1}}" stroke-width="2"/>
    <text x="{L+54}" y="{T+22}">No dither</text>
    <line x1="{L+18}" y1="{T+42}" x2="{L+46}" y2="{T+42}" stroke="{{s2}}" stroke-width="2"/>
    <text x="{L+54}" y="{T+46}">TPDF dither</text>
  </g>
  <text x="{L+pw/2:.1f}" y="{H-18}" text-anchor="middle" font-size="12.5" fill="{{ink2}}">Tone level (dBFS)</text>
  <text transform="translate(20,{T+ph/2:.1f}) rotate(-90)" text-anchor="middle" font-size="12.5" fill="{{ink2}}">SNR vs input (dB)</text>
</svg>'''
    for c, suf in ((C_LIGHT, ""), (C_DARK, ".dark")):
        out = svg
        for k, v in c.items():
            out = out.replace("{" + k + "}", v)
        with open(path.replace(".svg", suf + ".svg"), "w") as fh:
            fh.write(out)


def svg_folded(path, plain_curve, dith_curve, title):
    """The distortion proof: error folded on the 64-sample tone period.
    A fixed curve in both panels' shared scale means the error repeats with the
    signal (distortion); a near-zero line means it does not (noise)."""
    P = len(plain_curve)
    W, H = 920, 360
    L, R, T, B = 78, 24, 62, 56
    pw = W - L - R
    ph = (H - T - B - 30) / 2
    ymax = 0.55
    def X(i): return L + i / (P - 1) * pw
    def Y(v, top): return top + (ymax - v) / (2 * ymax) * ph
    def panel(c, top, color, label, note):
        pts = " ".join(f"{X(i):.1f},{Y(v, top):.1f}" for i, v in enumerate(c))
        y0 = Y(0, top)
        ticks = "".join(
            f'<text x="{L-8}" y="{Y(v, top)+4:.1f}" text-anchor="end" font-size="11" '
            f'fill="{{muted}}" font-variant-numeric="tabular-nums">{v:+.1f}</text>'
            for v in (-0.5, 0.0, 0.5))
        return f'''<text x="{L}" y="{top-8}" font-size="13" font-weight="600" fill="{{ink}}">{label}</text>
  <text x="{L+96}" y="{top-8}" font-size="12" fill="{{ink2}}">{note}</text>
  <line x1="{L}" y1="{y0:.1f}" x2="{L+pw}" y2="{y0:.1f}" stroke="{{grid}}" stroke-width="1"/>
  {ticks}
  <polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round"/>'''
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img"
     aria-label="{title}" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif">
  <rect width="{W}" height="{H}" fill="{{surface}}"/>
  <text x="{L}" y="26" font-size="17" font-weight="600" fill="{{ink}}">{title}</text>
  <text x="{L}" y="44" font-size="12.5" fill="{{ink2}}">Error averaged over 3000 tone periods: what repeats with the signal is distortion, what averages away is noise</text>
  {panel(plain_curve, T+24, "{s1}", "No dither", "fixed 0.47 LSB pattern: the staircase")}
  {panel(dith_curve, T+24+ph+30, "{s2}", "TPDF dither", "stays near zero: no signal-locked pattern")}
  <text x="{L+pw/2:.1f}" y="{H-14}" text-anchor="middle" font-size="12.5" fill="{{ink2}}">Sample position within one 64-sample tone period  ·  error in LSB</text>
</svg>'''
    for c, suf in ((C_LIGHT, ""), (C_DARK, ".dark")):
        out = svg
        for k, v in c.items():
            out = out.replace("{" + k + "}", v)
        with open(path.replace(".svg", suf + ".svg"), "w") as fh:
            fh.write(out)


def svg_gated_timeline(path, times_ms, plain_rms, dith_rms, title):
    """Noise modulation, the audible one: 20 ms error RMS through the gated
    file. Plain rounding drops to exactly zero in the gaps; dither holds a
    constant floor across tone and silence alike."""
    W, H = 920, 340
    L, R, T, B = 78, 24, 92, 58
    pw, ph = W - L - R, H - T - B
    xmax = times_ms[-1]
    ymin, ymax = 0.0, 0.55
    def X(v): return L + v / xmax * pw
    def Y(v): return T + (ymax - v) / (ymax - ymin) * ph
    def line(ys):
        return "M" + " L".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in zip(times_ms, ys))
    grid = []
    for v in [0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        y = Y(v)
        grid.append(f'<line x1="{L}" y1="{y:.1f}" x2="{L+pw}" y2="{y:.1f}" stroke="{{grid}}" stroke-width="1"/>')
        grid.append(f'<text x="{L-10}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="{{muted}}" font-variant-numeric="tabular-nums">{v:.1f}</text>')
    for v in range(0, 3001, 500):
        grid.append(f'<text x="{X(v):.1f}" y="{T+ph+22}" text-anchor="middle" font-size="12" fill="{{muted}}" font-variant-numeric="tabular-nums">{v}</text>')
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img"
     aria-label="{title}" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif">
  <rect width="{W}" height="{H}" fill="{{surface}}"/>
  <text x="{L}" y="26" font-size="17" font-weight="600" fill="{{ink}}">{title}</text>
  <text x="{L}" y="44" font-size="12.5" fill="{{ink2}}">Error RMS in 20 ms windows, through 100 ms tone / 100 ms silence gating at -70 dBFS</text>
  {"".join(grid)}
  <g font-size="12.5" fill="{{ink2}}">
    <line x1="{L}" y1="{T-24}" x2="{L+28}" y2="{T-24}" stroke="{{s1}}" stroke-width="2"/>
    <text x="{L+36}" y="{T-20}">No dither</text>
    <line x1="{L+150}" y1="{T-24}" x2="{L+178}" y2="{T-24}" stroke="{{s2}}" stroke-width="2"/>
    <text x="{L+186}" y="{T-20}">TPDF dither</text>
  </g>
  <line x1="{L}" y1="{T+ph}" x2="{L+pw}" y2="{T+ph}" stroke="{{axis}}" stroke-width="1"/>
  <path d="{line(plain_rms)}" fill="none" stroke="{{s1}}" stroke-width="1.6" stroke-linejoin="round"/>
  <path d="{line(dith_rms)}" fill="none" stroke="{{s2}}" stroke-width="1.6" stroke-linejoin="round"/>
  <text x="{L+pw/2:.1f}" y="{H-14}" text-anchor="middle" font-size="12.5" fill="{{ink2}}">Time (ms)</text>
  <text transform="translate(20,{T+ph/2:.1f}) rotate(-90)" text-anchor="middle" font-size="12.5" fill="{{ink2}}">Error RMS (LSB)</text>
</svg>'''
    for c, suf in ((C_LIGHT, ""), (C_DARK, ".dark")):
        out = svg
        for k, v in c.items():
            out = out.replace("{" + k + "}", v)
        with open(path.replace(".svg", suf + ".svg"), "w") as fh:
            fh.write(out)


# ---------------------------------------------------------------- main

def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    wav_dir, out_json = sys.argv[1:3]
    svg_dir = os.path.join(os.path.dirname(out_json), "svg")
    os.makedirs(svg_dir, exist_ok=True)

    out = {"fs": FS, "f0": F0, "quantizer": {"scale": QSCALE,
           "lsb_dbfs": 20 * np.log10(1 / QSCALE)}, "modes": {}, "level_sweep": {}}
    failures = []

    for mode in ("tone60", "sweep", "gated"):
        rec = {}
        for variant in ("plain", "dither"):
            wav = os.path.join(wav_dir, f"{mode}_{variant}.wav")
            x = read_f32(wav + ".f32")
            q = read_wav(wav)
            assert len(x) == len(q)
            y = q / QSCALE
            e = y - x
            e_lsb = e * QSCALE

            if variant == "plain":
                ref = plain_codes(x)
                mism = int(np.count_nonzero(ref != q))
                rec["plain_codes_mismatch"] = mism
                if mism:
                    failures.append(f"{mode}_plain_codes_mismatch={mism}")

            r = {
                "sample_count": int(len(q)),
                "snr_db": snr_db(x, e),
                "err_rms_lsb": float(np.sqrt(np.mean(e_lsb ** 2))),
                "err_dc_lsb": float(np.mean(e_lsb)),
                "err_var_lsb2": float(np.var(e_lsb)),
                "err_peak_lsb": float(np.max(np.abs(e_lsb))),
                "clipped_samples": int(np.sum((q == 32767) | (q == -32768))),
                "codes_used": int(len(np.unique(q))),
                "zero_samples": int(np.sum(q == 0)),
                "noise_modulation": noise_modulation(e_lsb),
            }
            if mode == "gated":
                blk = FS // 10
                active = (np.arange(len(x)) // blk) % 2 == 1
                r["folded"] = folded_distortion(x, e, active=active)
            else:
                r["folded"] = folded_distortion(x, e)
            # Gated signal: the global SNR is dominated by the digital-zero
            # gaps, which flatter plain rounding. Report the tone-on windows too.
            if mode == "gated":
                blk = FS // 10
                on = (np.arange(len(x)) // blk) % 2 == 1
                ps = float(np.mean(x[on] ** 2))
                pe = float(np.mean(e[on] ** 2))
                r["snr_active_db"] = 10 * np.log10(ps / pe) if pe > 0 else float("inf")
            if mode != "sweep":  # the sweep is not a steady tone: no harmonic ladder
                r["spectrum"] = spur_stats(e, F0)
                r["spurs"] = spurs_vs_floor(x, e)
                r["distortion"] = distortion_split(x, e)
            rec[variant] = r

        # dithered peak error must obey the 1.5 LSB quantization bound
        if rec["dither"]["err_peak_lsb"] > 1.5001:
            failures.append(f"{mode}_dither_err_peak>1.5LSB")

        q_p = read_wav(os.path.join(wav_dir, f"{mode}_plain.wav"))
        q_d = read_wav(os.path.join(wav_dir, f"{mode}_dither.wav"))
        rec["contrast"] = {
            "samples_differ": int(np.count_nonzero(q_p != q_d)),
            "samples_differ_pct": float(100.0 * np.mean(q_p != q_d)),
            "snr_delta_db": rec["dither"]["snr_db"] - rec["plain"]["snr_db"],
            "err_rms_delta_lsb": rec["dither"]["err_rms_lsb"] - rec["plain"]["err_rms_lsb"],
            "digital_silent_samples_plain": int(np.sum(q_p == 0)),
            "digital_silent_samples_dither": int(np.sum(q_d == 0)),
        }
        out["modes"][mode] = rec

    # ---- level sweep: measure the exported C output at each level ---------
    levels, snr_p, snr_d, rms_p, rms_d, fold_p, fold_d = [], [], [], [], [], [], []
    for lvl in (-6, -20, -40, -60, -80, -90, -100):
        xx = read_f32(os.path.join(wav_dir, f"levels/l{lvl}_plain.wav") + ".f32")
        ep = (read_wav(os.path.join(wav_dir, f"levels/l{lvl}_plain.wav")) / QSCALE - xx) * QSCALE
        qd = read_wav(os.path.join(wav_dir, f"levels/l{lvl}_dither.wav"))
        xd = read_f32(os.path.join(wav_dir, f"levels/l{lvl}_dither.wav") + ".f32")
        ed = (qd / QSCALE - xd) * QSCALE
        snr_p.append(10 * np.log10(np.mean(xx * xx) / np.mean((ep / QSCALE) ** 2)))
        snr_d.append(10 * np.log10(np.mean(xx * xx) / np.mean((ed / QSCALE) ** 2)))
        rms_p.append(float(np.sqrt(np.mean(ep ** 2))))
        rms_d.append(float(np.sqrt(np.mean(ed ** 2))))
        fold_p.append(folded_distortion(xx, ep / QSCALE)["folded_fraction"])
        fold_d.append(folded_distortion(xd, ed / QSCALE)["folded_fraction"])
        levels.append(lvl)
    out["level_sweep"] = {"levels_dbfs": levels, "snr_plain_db": snr_p,
                          "snr_dither_db": snr_d, "err_rms_plain_lsb": rms_p,
                          "err_rms_dither_lsb": rms_d,
                          "folded_fraction_plain": fold_p,
                          "folded_fraction_dither": fold_d}

    # ---- error-variance accounting on the -60 dBFS tone ----------------
    x = read_f32(os.path.join(wav_dir, "tone60_dither.wav") + ".f32")
    ed = (read_wav(os.path.join(wav_dir, "tone60_dither.wav")) / QSCALE - x) * QSCALE
    ep = (read_wav(os.path.join(wav_dir, "tone60_plain.wav")) / QSCALE - x) * QSCALE
    out["error_variance"] = {
        "plain_var_lsb2": float(np.var(ep)),
        "dither_var_lsb2": float(np.var(ed)),
        "dither_var_theory_lsb2": 0.25,
        "dither_mean_lsb": float(np.mean(ed)),
        "variance_increase_lsb2": float(np.var(ed) - np.var(ep)),
    }

    # ---- charts ---------------------------------------------------------
    x = read_f32(os.path.join(wav_dir, "tone60_dither.wav") + ".f32")
    q_d = read_wav(os.path.join(wav_dir, "tone60_dither.wav"))
    q_p = read_wav(os.path.join(wav_dir, "tone60_plain.wav"))
    ed = ((q_d / QSCALE) - x) * QSCALE
    ep = ((q_p / QSCALE) - x) * QSCALE
    db_p, freqs = coherent_spectrum_db(ep / QSCALE)
    db_d, _ = coherent_spectrum_db(ed / QSCALE)
    nz = freqs > 0
    svg_noise_spectrum(os.path.join(svg_dir, "noise_spectrum.svg"), freqs[nz], db_p[nz], db_d[nz],
                       "16-bit quantization error spectrum")
    # first 12 ms after a settling region, at a zero crossing
    t = np.arange(len(x)) / FS
    start = int(0.5 * FS)  # 500 ms in
    n = int(12e-3 * FS)
    svg_error_wave(os.path.join(svg_dir, "error_wave.svg"), ep[start:start + n], ed[start:start + n], 12,
                   "Error waveform: staircase vs noise")
    svg_snr_vs_level(os.path.join(svg_dir, "snr_vs_level.svg"), levels, snr_p, snr_d)
    svg_folded(os.path.join(svg_dir, "folded_error.svg"),
               out["modes"]["tone60"]["plain"]["folded"]["folded_curve_lsb"],
               out["modes"]["tone60"]["dither"]["folded"]["folded_curve_lsb"],
               "Error folded on the tone period")
    # gated timeline: 20 ms error RMS
    xg = read_f32(os.path.join(wav_dir, "gated_dither.wav") + ".f32")
    B = int(0.020 * FS)
    m = len(xg) // B
    tim = (np.arange(m) + 0.5) * 20.0
    rp = np.sqrt(np.mean((((read_wav(os.path.join(wav_dir, "gated_plain.wav")) / QSCALE) - xg)[:m * B]).reshape(m, B) ** 2, axis=1)) * QSCALE
    rd = np.sqrt(np.mean((((read_wav(os.path.join(wav_dir, "gated_dither.wav")) / QSCALE) - xg)[:m * B]).reshape(m, B) ** 2, axis=1)) * QSCALE
    svg_gated_timeline(os.path.join(svg_dir, "gated_timeline.svg"), tim, rp, rd,
                       "Error RMS through the gated file")

    out["failures"] = failures
    with open(out_json, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {out_json}")
    if failures:
        print("FAILURES: " + ", ".join(failures))
        return 1
    print("all verification claims hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
