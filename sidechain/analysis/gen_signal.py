#!/usr/bin/env python3
"""Generate a labeled sidechain test signal.

Background: Gaussian white noise at a fixed RMS.
Events: separated sine pings with random start, duration, level and frequency,
each with 5 ms raised-cosine fade in/out. Ground-truth labels are written per
frame (10 ms), 1 for any frame touched by an event.

Outputs (little-endian, in the directory given as argv[1]):
  key.f32          mono float32 key signal
  labels.u8        one label byte per 10 ms frame
  meta.json        generation parameters and the exact event list
"""
import json
import os
import sys

import numpy as np

SR = 48_000
DURATION_S = 60.0
FRAME_MS = 10.0
SEED = 42

NOISE_RMS_DBFS = -40.0
N_EVENTS = 60
DUR_MS = (80.0, 400.0)
LEVEL_DBFS = (-24.0, -12.0)          # uniform event level range
FREQ_HZ = (180.0, 2_000.0)
FADE_MS = 5.0
MIN_GAP_MS = 250.0


def main(out_dir: str) -> None:
    rng = np.random.default_rng(SEED)
    n = int(SR * DURATION_S)
    frame = int(SR * FRAME_MS / 1000.0)
    n_frames = n // frame

    noise_rms = 10.0 ** (NOISE_RMS_DBFS / 20.0)
    x = rng.normal(0.0, noise_rms, n).astype(np.float32)

    fade_n = int(SR * FADE_MS / 1000.0)
    min_gap = int(SR * MIN_GAP_MS / 1000.0)
    events = []
    # rejection sampling keeps events non-overlapping with a minimum gap
    attempts = 0
    while len(events) < N_EVENTS and attempts < 100_000:
        attempts += 1
        dur = int(SR * rng.uniform(*DUR_MS) / 1000.0)
        start = int(rng.integers(0, n - dur))
        end = start + dur
        if any(not (end + min_gap <= e["start"] or start >= e["end"] + min_gap)
               for e in events):
            continue
        level = 10.0 ** (rng.uniform(*LEVEL_DBFS) / 20.0)
        freq = float(rng.uniform(*FREQ_HZ))
        t = np.arange(dur) / SR
        tone = np.sin(2.0 * np.pi * freq * t)
        env = np.ones(dur)
        ramp = 0.5 * (1.0 - np.cos(np.pi * np.arange(fade_n) / fade_n))
        env[:fade_n] = ramp
        env[-fade_n:] = ramp[::-1]
        x[start:end] += (level * env * tone).astype(np.float32)
        events.append({"start": start, "end": end, "level_dbfs":
                       round(20.0 * np.log10(level), 3), "freq_hz": round(freq, 1)})

    events.sort(key=lambda e: e["start"])
    labels = np.zeros(n_frames, dtype=np.uint8)
    for e in events:
        f0, f1 = e["start"] // frame, (e["end"] - 1) // frame + 1
        labels[f0:f1] = 1

    x = np.clip(x, -1.0, 1.0)
    os.makedirs(out_dir, exist_ok=True)
    x.astype("<f4").tofile(os.path.join(out_dir, "key.f32"))
    labels.tofile(os.path.join(out_dir, "labels.u8"))
    meta = {
        "sample_rate": SR, "duration_s": DURATION_S, "frame_ms": FRAME_MS,
        "seed": SEED, "noise_rms_dbfs": NOISE_RMS_DBFS, "n_events": len(events),
        "events": events,
    }
    with open(os.path.join(out_dir, "meta.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print(f"samples={n} frames={n_frames} events={len(events)} "
          f"event_frames={int(labels.sum())} noise_frames={int((labels == 0).sum())}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
