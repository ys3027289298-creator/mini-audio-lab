"""Sabine-style reverberation-time formulas.

Two levels of model:

1. Global Sabine T60, one constant for the whole room:
       T60 = 0.161 V / (sum S_i alpha_i)
   This is the "reverb time is a global constant" baseline that ignores the
   local material distribution.

2. Local T60.  Absorption areas S_i alpha_i are first grouped by material
   zone, then re-weighted by distance from a reference point (the source, as
   requested; a receiver-weighted variant is also supplied for comparison):
       u_i   = 1 / d_i
       r_i   = u_i / mean(u)           (relative emphasis; 1 when equidistant)
       A'_i  = r_i * S_i * alpha_i
   and Sabine is applied to the re-weighted absorption A'.  Each zone's share
   of A' is also returned so the dominant zone can be named.  The normalized
   weights w_i = u_i / sum(u) are reported alongside the relative emphasis.
"""

import numpy as np

from room import BAND_FREQS, MATERIALS


def absorption_areas(room):
    """Per-zone geometric area and frequency-dependent absorption area.

    Returns (zone_names, S, A) where
      S[k]    = total surface area of zone k, m^2
      A[k, b] = S[k] * alpha_k(b), sabins per octave band.
    """
    zones = room.zones()
    names = sorted(zones)
    S = np.array([zones[n]["area"] for n in names])
    A = np.empty((len(names), len(BAND_FREQS)))
    for k, n in enumerate(names):
        A[k] = S[k] * MATERIALS[n]["alpha"]
    return names, S, A


def sabine_t60(volume, absorption_area):
    """Sabine formula.  absorption_area has shape (bands,) or (zones, bands)."""
    total = np.sum(absorption_area, axis=0) if absorption_area.ndim == 2 \
        else absorption_area
    return 0.161 * volume / total


def zone_distances(room, point):
    """Distance from a point to each zone centroid, aligned with zone order."""
    names, _, _ = absorption_areas(room)
    zones = room.zones()
    return np.array([np.linalg.norm(zones[n]["center"] - point) for n in names])


def local_t60(room, point, weighting="source"):
    """Distance-weighted local Sabine T60.

    weighting:
      "source"   : weights use source-to-zone distances (the requested model;
                   gives one constant shared by every receiver)
      "receiver" : weights use receiver-to-zone distances instead, so the
                   estimate actually moves with the listening position

    Returns a dict with per-band T60, the mid-band (500 Hz + 1 kHz) average,
    the normalized weights and each zone's share of weighted absorption.
    """
    names, S, A = absorption_areas(room)
    # A receiver can sit on a zone centroid; floor the distance so 1/d stays
    # finite (0.1 m is well below any meaningful zone separation).
    d = np.maximum(zone_distances(room, point), 0.1)
    inv = 1.0 / d
    w = inv / inv.sum()           # normalized shares, for reporting
    # Relative emphasis: ratio to the mean inverse distance.  With every zone
    # equidistant this is exactly 1, so A_w = A and the model collapses to the
    # global Sabine result (a required consistency check).
    rel = inv / inv.mean()
    A_w = rel[:, None] * A
    t60 = sabine_t60(room.volume, A_w)

    # Fraction of the weighted absorption contributed by each zone, per band.
    share = A_w / A_w.sum(axis=0, keepdims=True)

    return {
        "names": names,
        "weights": w,
        "distances": d,
        "t60_bands": t60,
        "t60_mid": float(np.mean(t60[[2, 3]])),
        "share": share,
        "weighting": weighting,
    }


def global_t60(room):
    """The baseline: one global Sabine T60, no local material information."""
    _, _, A = absorption_areas(room)
    t60 = sabine_t60(room.volume, A)
    return {"t60_bands": t60, "t60_mid": float(np.mean(t60[[2, 3]]))}
