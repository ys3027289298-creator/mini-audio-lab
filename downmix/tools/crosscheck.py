#!/usr/bin/env python3
"""Independent cross-check for the C downmix module (stdlib only).

Two layers:

  A. On the float32 frames dumped by the C self-test: recompute the matrix
     path in Python float64 and check the same identities. The C test tones
     are quantized to float32, which breaks their exact DFT-bin orthogonality
     at the ~1e-7 relative level, so tolerances there reflect float32 reality.

  B. On float64 tones generated inside this script: the same identities hold
     to machine precision, proving the matrix itself is exactly
     energy-conserving and that layer-A deviations come from input
     quantization, not the math.

Run from anywhere with: python3 -I tools/crosscheck.py
"""
import math
import os
import struct
import sys

N = 256
CH = 6
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
failures = 0
LF_BINS = [0, 1, 2, 3, 4, N - 4, N - 3, N - 2, N - 1]
NAMES = ("L", "R", "C", "LFE", "Ls", "Rs")


def check(name, ok):
    global failures
    print(f"{name}={1 if ok else 0}")
    if not ok:
        failures += 1


def read_f32(path):
    with open(os.path.join(HERE, path), "rb") as f:
        raw = f.read()
    return list(struct.unpack(f"<{len(raw) // 4}f", raw))


def dft_bin(x, k):
    re = im = 0.0
    w = -2.0 * math.pi * k / N
    for n, v in enumerate(x):
        a = w * n
        re += v * math.cos(a)
        im += v * math.sin(a)
    return re, im


def lf_energy(x):
    """Energy in DFT bins 0..4 (+ conjugate mirrors); Parseval: E = sum|X|^2/N."""
    return sum(re * re + im * im
               for re, im in (dft_bin(x, k) for k in LF_BINS)) / N


def matrix():
    q = 1.0 / math.sqrt(2.0)
    c, s = math.cos(math.pi / 6.0), 0.5
    wl = [1.0, 0.0, q, q, c, s]          # Lt row
    wr = [0.0, 1.0, q, q, s, c]          # Rt row
    return wl, wr


def apply_mx(wl, wr, x):
    lt = [sum(wl[j] * x[j][i] for j in range(CH)) for i in range(N)]
    rt = [sum(wr[j] * x[j][i] for j in range(CH)) for i in range(N)]
    return lt, rt


def tones(quantize):
    """Six orthogonal DFT-bin tones; quantize=True rounds each to float32."""
    bins = (20, 24, 12, 2, 18, 16)
    amps = (0.5, 0.5, 0.5, 0.8, 0.5, 0.5)
    out = []
    for k, a in zip(bins, amps):
        ch = []
        for i in range(N):
            v = a * math.cos(2.0 * math.pi * k * i / N)
            ch.append(f32(v) if quantize else v)
        out.append(ch)
    return out


def f32(v):
    return struct.unpack("<f", struct.pack("<f", v))[0]


def layer_a_float32_dump():
    print("== layer A: C float32 dump, recomputed in float64 ==")
    flat = read_f32("downmix_in.f32")
    x = [[float(v) for v in flat[c * N:(c + 1) * N]] for c in range(CH)]
    lt_c = [float(v) for v in read_f32("downmix_lt.f32")]
    rt_c = [float(v) for v in read_f32("downmix_rt.f32")]

    wl, wr = matrix()
    lt_py, rt_py = apply_mx(wl, wr, x)
    maxd = max(max(abs(lt_py[i] - lt_c[i]), abs(rt_py[i] - rt_c[i]))
               for i in range(N))
    print(f"max |C float32 output - python float64| = {maxd:.3e}")
    check("a_matches_c_samples", maxd < 1e-5)

    e_lfe = lf_energy(x[3])
    e_stereo = lf_energy(lt_py) + lf_energy(rt_py)
    mono = [(lt_py[i] + rt_py[i]) / math.sqrt(2.0) for i in range(N)]
    e_mono = lf_energy(mono)
    print(f"LF energy: LFE in={e_lfe:.6f} stereo={e_stereo:.6f} mono={e_mono:.6f}")
    # float32-quantized inputs: orthogonality leaks at ~1e-7 relative
    check("a_lf_stereo_conserved", abs(e_stereo - e_lfe) / e_lfe < 1e-6)
    check("a_lf_mono_preserved", abs(e_mono - e_lfe) / e_lfe < 1e-6)

    wr_bad = wr.copy()
    wr_bad[3] = -wr_bad[3]
    blt, brt = apply_mx(wl, wr_bad, x)
    bad_mono = [(blt[i] + brt[i]) / math.sqrt(2.0) for i in range(N)]
    e_bad_mono = lf_energy(bad_mono)
    e_bad_stereo = lf_energy(blt) + lf_energy(brt)
    print(f"opposite-polarity LFE: stereo={e_bad_stereo:.6f} "
          f"mono={e_bad_mono:.3e}")
    check("a_buggy_mono_cancels", e_bad_mono < 1e-6 * e_lfe)
    check("a_buggy_stereo_hides", abs(e_bad_stereo - e_lfe) / e_lfe < 1e-6)

    e_in = sum(v * v for ch in x for v in ch)
    e_out = sum(lt_py[i] ** 2 + rt_py[i] ** 2 for i in range(N))
    print(f"total energy: in={e_in:.6f} out={e_out:.6f} "
          f"relerr={abs(e_out - e_in) / e_in:.3e}")
    check("a_total_energy_float32_signal", abs(e_out - e_in) / e_in < 1e-6)

    with open(os.path.join(HERE, "bypass_in.f32"), "rb") as f:
        a_bytes = f.read()
    with open(os.path.join(HERE, "bypass_out.f32"), "rb") as f:
        b_bytes = f.read()
    check("a_bypass_byte_identical", a_bytes == b_bytes and len(a_bytes) > 0)


def layer_b_exact_double():
    print("== layer B: float64 tones, identities at machine precision ==")
    wl, wr = matrix()
    check("b_column_norms_unit", all(
        abs(wl[j] ** 2 + wr[j] ** 2 - 1.0) <= 1e-15 for j in range(CH)))
    check("b_lfe_in_phase", wl[3] > 0.0 and wr[3] > 0.0)

    x = tones(False)
    lt, rt = apply_mx(wl, wr, x)

    # exact orthogonality of distinct DFT bins => instantaneous identity
    e_in = sum(v * v for ch in x for v in ch)
    e_out = sum(lt[i] ** 2 + rt[i] ** 2 for i in range(N))
    print(f"total energy: in={e_in:.12f} out={e_out:.12f} "
          f"relerr={abs(e_out - e_in) / e_in:.3e}")
    check("b_total_energy_identity", abs(e_out - e_in) / e_in < 1e-12)

    e_lfe = lf_energy(x[3])
    e_stereo = lf_energy(lt) + lf_energy(rt)
    mono = [(lt[i] + rt[i]) / math.sqrt(2.0) for i in range(N)]
    e_mono = lf_energy(mono)
    print(f"LF energy: LFE={e_lfe:.12f} stereo={e_stereo:.12f} mono={e_mono:.12f}")
    check("b_lf_stereo_conserved", abs(e_stereo - e_lfe) / e_lfe < 1e-12)
    check("b_lf_mono_preserved", abs(e_mono - e_lfe) / e_lfe < 1e-12)

    # mono-fold gain check for every channel: (wL+wR)/sqrt(2) must be nonzero
    folds = [(wl[j] + wr[j]) / math.sqrt(2.0) for j in range(CH)]
    print("mono-fold gains: " +
          " ".join(f"{NAMES[j]}={folds[j]:+.6f}" for j in range(CH)))
    check("b_no_channel_cancels_in_mono", all(abs(g) > 0.5 for g in folds))

    wr_bad = wr.copy()
    wr_bad[3] = -wr_bad[3]
    blt, brt = apply_mx(wl, wr_bad, x)
    bad_mono = [(blt[i] + brt[i]) / math.sqrt(2.0) for i in range(N)]
    e_bad_mono = lf_energy(bad_mono)
    e_bad_stereo = lf_energy(blt) + lf_energy(brt)
    print(f"opposite-polarity LFE: stereo={e_bad_stereo:.12f} "
          f"mono={e_bad_mono:.3e}")
    check("b_buggy_mono_cancels", e_bad_mono < 1e-12 * e_lfe)
    check("b_buggy_stereo_hides", abs(e_bad_stereo - e_lfe) / e_lfe < 1e-12)

    # same float64 tones through f32 quantization reproduces layer-A error size,
    # pinpointing quantization as the only error source
    xq = tones(True)
    qlt, qrt = apply_mx(wl, wr, [[float(v) for v in ch] for ch in xq])
    eq_in = sum(v * v for ch in xq for v in ch)
    eq_out = sum(qlt[i] ** 2 + qrt[i] ** 2 for i in range(N))
    rel = abs(eq_out - eq_in) / eq_in
    print(f"float32-quantized tones identity relerr={rel:.3e} (expected ~1e-7)")
    check("b_quantization_explains_layer_a", 1e-9 < rel < 1e-5)


def main():
    layer_a_float32_dump()
    layer_b_exact_double()
    print(f"py_all_ok={1 if failures == 0 else 0}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
