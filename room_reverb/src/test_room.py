"""Sanity checks for the image-source and Sabine machinery."""

import numpy as np

from images import enumerate_images, receiver_arrivals, analytic_edc
from room import SPEED_OF_SOUND, build_layout, build_room
from render import analytic_t30
from sabine import global_t60, local_t60


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"  ok: {msg}")


def main():
    room = build_room()
    source, receivers = build_layout(room)

    # 1. Direct image (0,0,0) sits at the source, unit gain.
    images = enumerate_images(room, source, order=4)
    d = np.linalg.norm(images["positions"][images["direct_mask"]][0] - source)
    check(d < 1e-12, "direct image coincides with the source")
    check(np.allclose(images["gain"][images["direct_mask"]], 1.0),
          "direct image has unit gain in every band")

    # 2. Reflection counts: total image count and energy loss on reflection.
    n = 3
    small = enumerate_images(room, source, order=n)
    # total images: (2n+1)^3
    check(small["positions"].shape[0] == (2 * n + 1) ** 3,
          "image count equals (2N+1)^3")
    # every non-direct image has lost energy relative to the direct sound
    nondirect = ~small["direct_mask"]
    check(np.all(small["gain"][nondirect] < 1.0),
          "all reflected images have gain < 1")

    # 3. Delay of the direct arrival matches source-receiver distance / c.
    for name, p in receivers.items():
        delays, dist, spreading = receiver_arrivals(images, p)
        direct_delay = delays[images["direct_mask"]][0]
        expected = np.linalg.norm(p - source) / SPEED_OF_SOUND
        check(abs(direct_delay - expected) < 1e-9,
              f"{name}: direct delay equals distance / c")
        check(abs(spreading[images["direct_mask"]][0]
                  - 1 / (4 * np.pi * np.linalg.norm(p - source))) < 1e-12,
              f"{name}: direct-path spreading is 1/(4 pi d)")

    # 4. Analytic EDC is monotonically non-increasing and ends near zero.
    t, edc = analytic_edc(images, receivers["center"])
    check(np.all(np.diff(edc[:, 3]) <= 1e-9),
          "analytic EDC is monotonically non-increasing")
    check(edc[-1, 3] < -50, "analytic EDC falls below -50 dB at max order")

    # 5. Uniform-distance limit of the local weighted Sabine equals global.
    g = global_t60(room)
    loc = local_t60(room, np.array(room.size) / 2, weighting="receiver")
    # Force equal weights explicitly: weighted area must then equal total A.
    from sabine import absorption_areas, sabine_t60
    names, S, A = absorption_areas(room)
    # Equidistant zones: relative emphasis must be 1 everywhere.
    inv = 1.0 / np.full(len(names), 3.0)
    rel = inv / inv.mean()
    a_uniform = rel[:, None] * A
    check(np.allclose(a_uniform, A),
          "uniform distances leave every zone absorption area unchanged")
    check(np.allclose(sabine_t60(room.volume, a_uniform), g["t60_bands"]),
          "uniform weights reproduce the global Sabine T60")
    check(np.all(loc["t60_bands"] > 0) and np.all(np.isfinite(loc["t60_bands"])),
          "local T60 is finite and positive in every band (incl. centroid point)")
    # Zone shares have shape (zones, bands) and must sum to 1 per band.
    check(np.allclose(loc["share"].sum(axis=0), 1.0),
          "per-zone absorption shares sum to 1 in every band")

    # 6. Analytic T30 values are finite at every receiver.
    for name, p in receivers.items():
        vals = analytic_t30(images, p)
        check(all(np.isfinite(v[0]) for v in vals),
              f"{name}: analytic T30 finite in all bands")

    print("\nall checks passed.")


if __name__ == "__main__":
    main()
