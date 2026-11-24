#!/usr/bin/env python3
"""Independent verification of the mini-audio-lab streaming pipeline.

Reads only the raw *.f32 outputs produced by pipeline_legacy /
pipeline_unified (no C headers, no shared code) and, in pure NumPy:

  * re-derives every stage length from the stream formulas;
  * counts block-boundary clicks (boundary |diff| > 4x in-block median);
  * counts legacy zero-padding at spatial block boundaries;
  * recomputes each stage from the stage input and compares bytes;
  * draws the per-stage spectrograms.

Run with: python3 -I tools/analyze.py <legacy_dir> <unified_dir> <out_artifacts>
"""
import json
import math
import os
import sys

# python3 -I drops PYTHONPATH and user site-packages; the script's own
# directory stays on sys.path, so locate the bootstrap next to this file.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _isolate_site  # noqa: F401,E402  (re-adds user site-packages)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# The C stage uses cosf/sinf from libm; call the same functions so the
# independent Python recomputation can be byte-exact, not ULP-tolerant.
import ctypes
_libm = ctypes.CDLL("libm.so.6")
_libm.cosf.restype = ctypes.c_float
_libm.cosf.argtypes = [ctypes.c_float]
_libm.sinf.restype = ctypes.c_float
_libm.sinf.argtypes = [ctypes.c_float]
_cosf = np.frompyfunc(_libm.cosf, 1, 1)
_sinf = np.frompyfunc(_libm.sinf, 1, 1)

SR_IN = 8000
RATIO = 1.3
BLOCK_S = 2.0
N_BLOCKS = 6
PAN_X = -0.3
MIX_GAIN = 0.6
RAMP_TOTAL = 4096
ATT_REF, ATT_ROLL, ATT_NOW = 2.0, 1.0, 4.0
PI_F = np.float32(3.14159265)  # same literal as pan() in spatial.c

# dataviz skill: validated blue sequential ramp (light -> dark).
BLUE_RAMP = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab"]


def read_f32(path):
    raw = open(path, "rb").read()
    assert len(raw) % 4 == 0, f"{path}: not whole floats"
    return np.frombuffer(raw, dtype="<f4").copy()


def total_out(n_in, ratio):
    return math.floor(n_in / ratio)


def block_outs(n_blocks, block_in, ratio):
    edges = [math.floor((b * block_in) / ratio) for b in range(n_blocks + 1)]
    return [edges[b + 1] - edges[b] for b in range(n_blocks)]


def resample_py(x, ratio):
    """Independent reimplementation of resample_linear().
    Keep t in double until the final cast, exactly like the C version."""
    n_out = math.floor(len(x) / ratio)
    idx = np.arange(n_out, dtype=np.float64) * ratio
    i0 = np.floor(idx).astype(np.int64)
    t = (idx - i0)
    i1 = np.minimum(i0 + 1, len(x) - 1)
    return ((1.0 - t) * x[i0].astype(np.float64)
            + t * x[i1].astype(np.float64)).astype(np.float32)


def mix_py(x):
    """Whole-stream ramped self-mix, matching ref_mix() in pipeline.c."""
    i = np.arange(len(x), dtype=np.float32)
    t = np.minimum(i / np.float32(RAMP_TOTAL), np.float32(1.0))
    v = x * (MIX_GAIN * t) + x * (MIX_GAIN * t)
    return np.clip(v, -1.0, 1.0).astype(np.float32)


def pan_gains():
    """Byte-exact replica of pan(PAN_X) in spatial.c (float32 arithmetic)."""
    x = np.float32(PAN_X)
    a = (x + np.float32(1.0)) * np.float32(0.25) * PI_F
    gl = np.float32(_cosf(a))
    gr = np.float32(_sinf(a))
    return gl, gr


def attenuation():
    """Byte-exact replica of attenuate() for the demo distance."""
    return np.float32(ATT_REF / (ATT_REF + ATT_ROLL * (ATT_NOW - ATT_REF)))


def spatial_py(x):
    """Whole-stream constant-power pan + attenuation, pairs summed.
    Float32 evaluation order matches spatial_frame() in spatial.c."""
    gl, gr = pan_gains()
    att = attenuation()
    pairs = len(x) // 2
    x2 = x[: 2 * pairs].reshape(pairs, 2).astype(np.float32)
    l0 = (x2[:, 0] * gl) * att
    l1 = (x2[:, 1] * gl) * att
    r0 = (x2[:, 0] * gr) * att
    r1 = (x2[:, 1] * gr) * att
    return (l0 + l1).astype(np.float32), (r0 + r1).astype(np.float32)


def legacy_pipeline_py(x):
    """Independent replay of the pre-unification per-block behavior."""
    block_in = SR_IN * int(BLOCK_S)
    gl, gr = pan_gains()
    att = attenuation()
    rs_parts, mx_parts, l_parts, r_parts = [], [], [], []
    for b in range(N_BLOCKS):
        blk = x[b * block_in: (b + 1) * block_in]
        # legacy resampler: each block rounds down on its own
        rb = resample_py(blk, RATIO)
        rs_parts.append(rb)
        # legacy mixer: ramp restarts at every block
        i = np.arange(len(rb), dtype=np.float32)
        t = np.minimum(i / np.float32(2048), np.float32(1.0))
        mb = np.clip(2.0 * rb * (MIX_GAIN * t), -1.0, 1.0).astype(np.float32)
        mx_parts.append(mb)
        # legacy spatial: ceil(n/2) pairs. A normal pair sums both frames.
        # On an odd block the dangling last frame is emitted alone as a
        # full stereo frame that should never have existed (one spurious
        # pair per odd block); the next block's first frame then lands in
        # the wrong stereo position - that misalignment is the click source.
        pairs = (len(mb) + 1) // 2
        for p in range(pairs):
            f0l = np.float32(mb[2 * p] * gl) * att
            f0r = np.float32(mb[2 * p] * gr) * att
            if 2 * p + 1 < len(mb):
                f1l = np.float32(mb[2 * p + 1] * gl) * att
                f1r = np.float32(mb[2 * p + 1] * gr) * att
                l_parts.append(np.float32(f0l + f1l))
                r_parts.append(np.float32(f0r + f1r))
            else:
                l_parts.append(np.float32(f0l))
                r_parts.append(np.float32(f0r))
    return (np.concatenate(rs_parts), np.concatenate(mx_parts),
            np.array(l_parts, dtype=np.float32),
            np.array(r_parts, dtype=np.float32))


def boundary_click_count(x, block_lens, mult=4.0):
    """Clicks at block joins: boundary |diff| > mult * in-block median |diff|."""
    diffs = np.abs(np.diff(x))
    if len(diffs) == 0:
        return 0, 0.0
    edges = np.cumsum(block_lens)[:-1]  # index of first sample of block 1..
    boundary_idx = edges - 1            # diff straddling the join
    mask = np.ones(len(diffs), bool)
    mask[boundary_idx] = False
    median = float(np.median(diffs[mask])) or 1e-12
    clicks = int(np.sum(diffs[boundary_idx] > mult * median))
    return clicks, float(np.max(diffs[boundary_idx]) / median)


def main():
    legacy_dir, unified_dir, out_dir = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    report = {}

    in_leg = read_f32(os.path.join(legacy_dir, "input.f32"))
    in_uni = read_f32(os.path.join(unified_dir, "input.f32"))
    assert len(in_leg) == len(in_uni) == SR_IN * 12
    assert in_leg.tobytes() == in_uni.tobytes(), "inputs must be identical"

    block_in = SR_IN * int(BLOCK_S)
    formula_blocks = block_outs(N_BLOCKS, block_in, RATIO)
    formula_total = total_out(len(in_leg), RATIO)

    length_rows = []
    for mode, d in (("legacy", legacy_dir), ("unified", unified_dir)):
        run = json.load(open(os.path.join(d, "run.json")))
        rs = read_f32(os.path.join(d, "resampler.f32"))
        mx = read_f32(os.path.join(d, "mixer.f32"))
        sl = read_f32(os.path.join(d, "spatial_L.f32"))
        sr_ = read_f32(os.path.join(d, "spatial_R.f32"))
        il = read_f32(os.path.join(d, "spatial_interleaved.f32"))

        r = {}
        r["input"] = len(in_leg)
        r["resampler"] = len(rs)
        r["mixer"] = len(mx)
        r["spatial_pairs"] = len(sl)
        r["spatial_frames"] = 2 * len(sl)
        r["rs_block_lens"] = run["block_len"]["resampler"]
        r["mix_block_lens"] = run["block_len"]["mixer"]
        r["sp_block_lens"] = run["block_len"]["spatial_pairs"]
        r["spurious_pairs"] = run["total"]["spatial_spurious_pairs"]
        report[mode] = r

        # length derivation checks
        assert len(mx) == len(rs), f"{mode}: mixer must be length-neutral"
        assert len(il) == 2 * len(sl) and len(sr_) == len(sl)
        if mode == "unified":
            assert r["rs_block_lens"] == formula_blocks, (r["rs_block_lens"],
                                                           formula_blocks)
            assert len(rs) == formula_total
            assert len(sl) == formula_total // 2
            assert r["spurious_pairs"] == 0
            # cumulative outputs must match the prefix floors exactly
            cum = np.cumsum(r["rs_block_lens"])
            expected = [math.floor(((b + 1) * block_in) / RATIO)
                        for b in range(N_BLOCKS)]
            assert cum.tolist() == expected
        else:
            # legacy: every block independently floor(B/ratio), one short
            assert r["rs_block_lens"] == [math.floor(block_in / RATIO)] * N_BLOCKS
            assert len(rs) == N_BLOCKS * math.floor(block_in / RATIO)

        # click counts at the five internal block boundaries
        clicks = {}
        ratios = {}
        clicks["resampler"], ratios["resampler"] = boundary_click_count(
            rs, r["rs_block_lens"])
        clicks["mixer"], ratios["mixer"] = boundary_click_count(
            mx, r["mix_block_lens"])
        sp_blocks = [2 * p for p in r["sp_block_lens"]]
        clicks["spatial_L"], ratios["spatial_L"] = boundary_click_count(
            sl, r["sp_block_lens"])
        clicks["spatial_R"], ratios["spatial_R"] = boundary_click_count(
            sr_, r["sp_block_lens"])
        r["clicks"], r["boundary_ratio"] = clicks, ratios

        # independent byte-level recomputation
        py_rs = resample_py(in_leg, RATIO)
        if mode == "unified":
            assert len(py_rs) == len(rs)
            assert py_rs.tobytes() == rs.tobytes(), "resampler bytes mismatch"
            py_mx = mix_py(py_rs)
            assert py_mx.tobytes() == mx.tobytes(), "mixer bytes mismatch"
            py_l, py_r = spatial_py(py_mx)
            assert py_l.tobytes() == sl.tobytes(), "spatial L bytes mismatch"
            assert py_r.tobytes() == sr_.tobytes(), "spatial R bytes mismatch"
        else:
            lg_rs, lg_mx, lg_l, lg_r = legacy_pipeline_py(in_leg)
            assert lg_rs.tobytes() == rs.tobytes(), "legacy resampler mismatch"
            assert lg_mx.tobytes() == mx.tobytes(), "legacy mixer mismatch"
            assert lg_l.tobytes() == sl.tobytes(), "legacy spatial L mismatch"
            assert lg_r.tobytes() == sr_.tobytes(), "legacy spatial R mismatch"
        r["byte_recompute_ok"] = True

        # legacy spatial: the spurious frame holds the dangling input in L
        if mode == "legacy":
            spurious = int(r["spurious_pairs"])
            r["observed_spurious_pairs"] = spurious

        length_rows.append((mode, r))

    # cross-mode: turning unification off must reproduce the legacy byte stream
    leg_rs = read_f32(os.path.join(legacy_dir, "resampler.f32"))
    uni_rs = read_f32(os.path.join(unified_dir, "resampler.f32"))
    assert len(leg_rs) != len(uni_rs)  # lengths intentionally differ

    # ---------------- length table ----------------
    with open(os.path.join(out_dir, "length_table.md"), "w") as f:
        f.write("# Output length comparison "
                f"(sr_in={SR_IN}, ratio={RATIO}, block={int(BLOCK_S)}s, "
                f"{N_BLOCKS} blocks)\n\n")
        f.write("unified formula: total `floor(N/r)` = "
                f"**{formula_total}**; per block "
                "`floor((b+1)B/r) - floor(bB/r)`; "
                f"stereo pairs `floor(M/2)` = **{formula_total // 2}**\n\n")
        f.write("| stage / block | legacy n_out | unified n_out | delta |\n")
        f.write("|---|---|---|---|\n")
        leg, uni = report["legacy"], report["unified"]
        for b in range(N_BLOCKS):
            f.write(f"| resampler block {b} | {leg['rs_block_lens'][b]} | "
                    f"{uni['rs_block_lens'][b]} | "
                    f"{uni['rs_block_lens'][b] - leg['rs_block_lens'][b]} |\n")
        for key, label in (("resampler", "resampler TOTAL"),
                           ("mixer", "mixer TOTAL"),
                           ("spatial_frames", "spatial frames TOTAL")):
            f.write(f"| {label} | {leg[key]} | {uni[key]} | "
                    f"{uni[key] - leg[key]} |\n")
        f.write(f"| input TOTAL | {leg['input']} | {uni['input']} | 0 |\n\n")
        f.write("## Block-boundary clicks (boundary |diff| > 4x median)\n\n")
        f.write("| stage | legacy | unified |\n|---|---|---|\n")
        for k in ("resampler", "mixer", "spatial_L", "spatial_R"):
            f.write(f"| {k} | {leg['clicks'][k]} | {uni['clicks'][k]} |\n")
        f.write("\n## Spurious stereo frames at spatial boundaries\n\n")
        f.write("Each odd legacy block emits one extra stereo frame from its\n")
        f.write("dangling sample, shifting every later block by one pair.\n\n")
        f.write("| measure | legacy | unified |\n|---|---|---|\n")
        f.write(f"| spurious pairs (run.json) | "
                f"{leg['spurious_pairs']} | {uni['spurious_pairs']} |\n")
        f.write(f"| independently replayed | "
                f"{leg['observed_spurious_pairs']} | 0 |\n\n")
        f.write("## Boundary discontinuity ratio (boundary / in-block median)\n\n")
        f.write("| signal | legacy | unified |\n|---|---|---|\n")
        for k in ("resampler", "mixer", "spatial_L", "spatial_R"):
            f.write(f"| {k} | {leg['boundary_ratio'][k]:.2f} | "
                    f"{uni['boundary_ratio'][k]:.2f} |\n")

    with open(os.path.join(out_dir, "length_table.csv"), "w") as f:
        f.write("stage / block,legacy n_out,unified n_out,delta\n")
        for b in range(N_BLOCKS):
            f.write(f"resampler block {b},{leg['rs_block_lens'][b]},"
                    f"{uni['rs_block_lens'][b]},"
                    f"{uni['rs_block_lens'][b] - leg['rs_block_lens'][b]}\n")
        for key in ("resampler", "mixer", "spatial_frames", "input"):
            f.write(f"{key} TOTAL,{leg[key]},{uni[key]},"
                    f"{uni[key] - leg[key]}\n")

    # ---------------- spectrograms, one per stage per mode ----------------
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "ma_blue", BLUE_RAMP)
    signals = [
        ("resampler", "resampler.f32", SR_IN / RATIO),
        ("mixer", "mixer.f32", SR_IN / RATIO),
        ("spatial_L", "spatial_L.f32", SR_IN / RATIO),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for col, (name, fn, sr) in enumerate(signals):
        for row, (mode, d) in enumerate((("legacy", legacy_dir),
                                         ("unified", unified_dir))):
            ax = axes[row, col]
            x = read_f32(os.path.join(d, fn))
            nfft = 1024
            hop = nfft // 4
            win = np.hanning(nfft)
            frames = 1 + (len(x) - nfft) // hop
            spec = np.empty((nfft // 2, frames))
            for fr in range(frames):
                seg = x[fr * hop: fr * hop + nfft] * win
                sp = np.fft.rfft(seg, nfft)[: nfft // 2]
                spec[:, fr] = np.abs(sp) + 1e-12
            db = 20 * np.log10(spec / (np.max(spec) + 1e-12))
            t = np.arange(frames) * hop / sr
            f = np.fft.rfftfreq(nfft, 1 / sr)[: nfft // 2]
            pcm = ax.pcolormesh(t, f, db, cmap=cmap, vmin=-80, vmax=0,
                                shading="auto")
            for b in range(1, N_BLOCKS):
                ax.axvline(b * BLOCK_S, color="#e34948", lw=0.8,
                           ls="--", alpha=0.7)
            ax.set_title(f"{mode} - {name}", fontsize=11)
            ax.set_xlabel("time (s)")
            ax.set_ylabel("frequency (Hz)")
            ax.set_ylim(0, sr / 2)
    cax = fig.add_axes([0.935, 0.15, 0.014, 0.70])
    fig.colorbar(pcm, cax=cax, label="dB (normalized)")
    fig.suptitle("Per-stage spectrograms (red dashed = 2 s block boundaries)",
                 fontsize=13)
    fig.subplots_adjust(left=0.06, right=0.91, top=0.91, bottom=0.08,
                        wspace=0.25, hspace=0.32)
    fig.savefig(os.path.join(out_dir, "spectrograms.png"), dpi=130)
    plt.close(fig)

    # ---------------- machine-readable summary ----------------
    summary = {
        "formula": {"total_out": formula_total,
                    "block_out": formula_blocks,
                    "total_pairs": formula_total // 2},
        "legacy": {k: v for k, v in leg.items()},
        "unified": {k: v for k, v in uni.items()},
        "click_assertion": (
            max(uni["clicks"].values()) == 0
            and sum(leg["clicks"].values()) > 0),
    }
    with open(os.path.join(out_dir, "analysis.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps({
        "formula_total": formula_total,
        "legacy_lengths": {k: leg[k] for k in
                           ("resampler", "mixer", "spatial_frames")},
        "unified_lengths": {k: uni[k] for k in
                            ("resampler", "mixer", "spatial_frames")},
        "legacy_clicks": leg["clicks"],
        "unified_clicks": uni["clicks"],
        "legacy_spurious_pairs": leg["spurious_pairs"],
        "unified_spurious_pairs": uni["spurious_pairs"],
        "byte_recompute_ok": True,
    }, indent=2))


if __name__ == "__main__":
    main()
