"""Specular image-source model for a rectangular room.

Every image source is identified by an integer triple (nx, ny, nz).  Its
distance to a receiver gives the reflection delay, and the number of times the
ray hits each face determines the reflection gain in each octave band:

    r_zone(band) = sqrt(1 - alpha_zone(band))

so an image's band gain is the product of one reflection coefficient per
wall hit.  The same image list serves two purposes:

  * analytic energy decay curves (exact for the specular field, no rendering
    noise), used to cross-check the T30 estimates taken from rendered IRs;
  * early-field rendering of band-limited impulse responses (render.py).
"""

import numpy as np

from room import BAND_FREQS, MATERIALS, SPEED_OF_SOUND


def _image_positions(n, L, s):
    """1-D image positions for image indices n along an axis of length L."""
    n = np.asarray(n)
    pos = np.empty(n.shape)
    even = n % 2 == 0
    pos[even] = n[even] * L + s
    odd = ~even
    k = (n[odd] - 1) // 2
    pos[odd] = 2 * (k + 1) * L - s
    return pos


def enumerate_images(room, source, order=30):
    """All image sources up to reflection `order` per axis.

    Returns arrays:
      delays    : arrival time at the source itself is not meaningful, so this
                  is computed per receiver; here we return positions instead.
      positions : (M, 3) image positions, m
      hit_counts: (M, 6) reflection counts, ordered axis0-, axis0+, axis1-, ...
      gain      : (M, B) amplitude gain per octave band (incl. no spreading)
    """
    size = room.size
    axis_rng = np.arange(-order, order + 1)
    nx, ny, nz = np.meshgrid(axis_rng, axis_rng, axis_rng, indexing="ij")
    n = np.stack([nx.ravel(), ny.ravel(), nz.ravel()], axis=1)  # (M, 3)
    M = n.shape[0]

    pos = np.empty((M, 3))
    for ax in range(3):
        pos[:, ax] = _image_positions(n[:, ax], size[ax], source[ax])

    # Reflection counts on the low (-1) and high (+1) face of each axis.
    counts = np.empty((M, 6), dtype=np.int64)
    for ax in range(3):
        counts[:, 2 * ax] = np.maximum(0, -n[:, ax])  # low face
        counts[:, 2 * ax + 1] = np.maximum(0, n[:, ax])  # high face

    # Map each of the six faces to its material zone.
    face_zone = [None] * 6
    for f in room.faces:
        face_zone[2 * f.axis + (1 if f.sign > 0 else 0)] = f.zone

    # log amplitude gain per band: 0.5 * log(1-alpha) accumulated per hit.
    log_gain = np.zeros((M, len(BAND_FREQS)))
    for face_idx, zone in enumerate(face_zone):
        c = counts[:, face_idx]
        if not np.any(c):
            continue
        log_r = 0.5 * np.log(1.0 - MATERIALS[zone]["alpha"])  # per band
        log_gain += c[:, None] * log_r[None, :]
    gain = np.exp(log_gain)

    direct = np.all(n == 0, axis=1)
    return {
        "positions": pos,
        "counts": counts,
        "gain": gain,
        "face_zone": face_zone,
        "direct_mask": direct,
        "order": order,
    }


def receiver_arrivals(images, receiver):
    """Delay (s) and geometric-spreading amplitude for one receiver."""
    delta = images["positions"] - receiver[None, :]
    dist = np.linalg.norm(delta, axis=1)
    delays = dist / SPEED_OF_SOUND
    spreading = 1.0 / (4 * np.pi * dist)
    return delays, dist, spreading


def band_energy_envelope(images, receiver):
    """Per-image energy (B bands) at the receiver, including 1/d^2 spreading."""
    _, _, spreading = receiver_arrivals(images, receiver)
    return (images["gain"] ** 2) * (spreading ** 2)[:, None]


def analytic_edc(images, receiver, t=None):
    """Schroeder energy decay curve from the discrete image set.

    EDC(t) in dB = 10 log10( sum of energies arriving after t / total ).
    """
    delays, _, _ = receiver_arrivals(images, receiver)
    energy = band_energy_envelope(images, receiver)  # (M, B)

    if t is None:
        # One sample just after every distinct arrival time.
        t = np.sort(delays)
    order = np.argsort(delays)
    sorted_e = energy[order]
    total = sorted_e.sum(axis=0)
    cum_after = total - np.cumsum(sorted_e, axis=0)

    edc = np.empty((len(t), energy.shape[1]))
    idx = np.searchsorted(delays[order], t, side="right")
    for b in range(energy.shape[1]):
        tail = np.where(idx < len(order),
                        cum_after[np.minimum(idx, len(order) - 1), b], 0.0)
        edc[:, b] = 10.0 * np.log10(np.maximum(tail, 1e-20) / total[b])
    return t, edc
