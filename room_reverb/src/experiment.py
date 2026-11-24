"""End-to-end room-reverberation experiment.

Builds a three-material room, measures impulse responses at three receivers,
estimates T60 (T30 fit on Schroeder curves) in six octave bands, and compares
them against:

  * the global Sabine constant (one T60 for the whole room);
  * the requested local formula: zone absorption weighted by source-to-zone
    distance (a single value shared by all receivers);
  * a receiver-weighted variant of the same formula, which can actually vary
    with listening position.

Outputs (printed tables + CSV + PNG) go to room_reverb/artifacts/.
Run:  make run   (or  python3 experiment.py from this directory)
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from images import enumerate_images, analytic_edc
from render import (FS, DURATION, measure_t30, measure_edt, analytic_t30,
                    render_ir, schroeder_edc)
from room import BAND_FREQS, MATERIALS, build_layout, build_room
from sabine import absorption_areas, global_t60, local_t60

ART = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                   "artifacts")

# Validated dataviz categorical palette (light surface); aqua gets direct
# labels/table relief for its sub-3:1 contrast.
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a"}
ZONE_COLOR = {"concrete": "#52514e", "wood": "#b07a3a", "foam": "#1baf7a"}
RCV_COLOR = {"center": C["blue"], "hard_corner": C["orange"],
             "foam_wall": C["aqua"]}
RCV_LABEL = {"center": "room centre", "hard_corner": "hard corner",
             "foam_wall": "in front of foam wall"}
IMAGE_ORDER = 30


def fmt(x):
    return "  n/a " if not np.isfinite(x) else f"{x:6.3f}"


def main():
    os.makedirs(ART, exist_ok=True)
    room = build_room()
    source, receivers = build_layout(room)
    names, S, A = absorption_areas(room)

    print("=" * 72)
    print("ROOM")
    print(f"  size {tuple(room.size)} m, volume {room.volume:.1f} m^3, "
          f"surface {room.surface_area:.1f} m^2")
    zones = room.zones()
    print("  zones:")
    for k, n in enumerate(names):
        print(f"    {n:9s} S={S[k]:5.1f} m^2  centroid="
              f"{np.round(zones[n]['center'], 2)}  "
              f"alpha={MATERIALS[n]['alpha']}")

    # ---- Formula side -----------------------------------------------------
    g = global_t60(room)
    src_local = local_t60(room, source, weighting="source")
    rcv_local = {r: local_t60(room, p, weighting="receiver")
                 for r, p in receivers.items()}

    print("\nSOURCE-DISTANCE WEIGHTS (the requested local model)")
    for k, n in enumerate(names):
        print(f"  {n:9s} d={src_local['distances'][k]:5.2f} m  "
              f"weight={src_local['weights'][k]:.3f}")
    print(f"  -> local Sabine T60 (mid-band) = {src_local['t60_mid']:.3f} s; "
          f"global Sabine = {g['t60_mid']:.3f} s")
    print("  note: source-weighted value is one constant for every receiver.")

    # ---- Measurement side --------------------------------------------------
    print("\nBuilding image sources ...")
    images = enumerate_images(room, source, order=IMAGE_ORDER)
    print(f"  {images['positions'].shape[0]} images (order {IMAGE_ORDER}/axis)")

    rows = []
    rendered = {}
    for rcv_name, rcv_pos in receivers.items():
        print(f"\n--- receiver: {rcv_name} {np.round(rcv_pos, 2)} ---")
        band_irs, broadband = render_ir(images, rcv_pos, seed=42)
        rendered[rcv_name] = (band_irs, broadband)

        t30_render, edt_render, r2 = [], [], []
        for b, fc in enumerate(BAND_FREQS):
            t60, q = measure_t30(band_irs[b], FS)
            edt, _ = measure_edt(band_irs[b], FS)
            t30_render.append(t60)
            edt_render.append(edt)
            r2.append(q["r2"])
        t30_exact = [v[0] for v in analytic_t30(images, rcv_pos)]

        for b, fc in enumerate(BAND_FREQS):
            rows.append({
                "receiver": rcv_name, "band_hz": int(fc),
                "t60_measured_t30": t30_render[b],
                "edt_measured": edt_render[b],
                "t60_analytic_ism": t30_exact[b],
                "t60_sabine_global": g["t60_bands"][b],
                "t60_sabine_source_weighted": src_local["t60_bands"][b],
                "t60_sabine_receiver_weighted":
                    rcv_local[rcv_name]["t60_bands"][b],
                "fit_r2": r2[b],
            })

        mid = [2, 3]
        measured_mid = float(np.mean([t30_render[i] for i in mid]))
        exact_mid = float(np.mean([t30_exact[i] for i in mid]))
        print("  band Hz   T30(IR)   EDT(IR)  T30(ISM EDC)  Sabine-global  "
              "Sabine-src-w  Sabine-rcv-w    fit R2")
        for b, fc in enumerate(BAND_FREQS):
            print(f"  {int(fc):5d}  {fmt(t30_render[b])}  "
                  f"{fmt(edt_render[b])}  {fmt(t30_exact[b])}  "
                  f"{fmt(g['t60_bands'][b])}  "
                  f"{fmt(src_local['t60_bands'][b])}  "
                  f"{fmt(rcv_local[rcv_name]['t60_bands'][b])}  "
                  f"{r2[b]:.4f}")
        print(f"  mid-band(500/1k) measured={measured_mid:.3f} s, "
              f"exact ISM={exact_mid:.3f} s, "
              f"src-weighted formula={src_local['t60_mid']:.3f} s, "
              f"rcv-weighted formula={rcv_local[rcv_name]['t60_mid']:.3f} s")

        # "Muffled" indicator: low-frequency energy relative to high frequency.
        e_lo = np.sum(band_irs[1] ** 2)  # 250 Hz band
        e_hi = np.sum(band_irs[4] ** 2)  # 2 kHz band
        print(f"  E(250 Hz)/E(2 kHz) = {10*np.log10(e_lo/e_hi):+.1f} dB "
              f"(higher = darker/muffler)")

    # ---- CSV ---------------------------------------------------------------
    csv_path = os.path.join(ART, "t60_comparison.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # Zone shares (dominant-zone attribution), receiver-weighted model.
    share_path = os.path.join(ART, "zone_shares.csv")
    with open(share_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["receiver", "zone", "distance_m", "weight",
                    *[f"share_{int(fc)}Hz" for fc in BAND_FREQS]])
        for rcv_name in receivers:
            loc = rcv_local[rcv_name]
            for k, n in enumerate(names):
                w.writerow([rcv_name, n, f"{loc['distances'][k]:.3f}",
                            f"{loc['weights'][k]:.4f}",
                            *[f"{v:.4f}" for v in loc["share"][k, :]]])

    # ---- Worst discrepancy vs requested formula ---------------------------
    # Ranking uses the noise-free analytic ISM measurement: rendered T30 is
    # shifted a little by pulse-width smearing (deterministic, verified stable
    # across six noise seeds), and the analytic EDC is the ground truth of the
    # modeled physics.  Rendered figures appear in the per-receiver tables.
    print("\n" + "=" * 72)
    print("DISCREPANCY: measured T30 vs source-weighted Sabine")
    summary = []
    for rcv_name in receivers:
        r = [x for x in rows if x["receiver"] == rcv_name]
        meas_render = np.array([x["t60_measured_t30"] for x in r])
        meas_exact = np.array([x["t60_analytic_ism"] for x in r])
        formula_all = np.array([x["t60_sabine_source_weighted"] for x in r])
        err_exact = float(np.mean(np.abs(meas_exact - formula_all)))
        err_render = float(np.mean(np.abs(meas_render - formula_all)))
        mid_exact = np.mean(meas_exact[[2, 3]])
        summary.append((rcv_name, mid_exact, err_exact, err_render,
                        float(mid_exact - src_local["t60_mid"])))
        print(f"  {rcv_name:12s} all-band mean |delta|: analytic={err_exact:.3f} s, "
              f"rendered={err_render:.3f} s; mid-band analytic delta "
              f"{mid_exact - src_local['t60_mid']:+.3f} s")
    worst = max(summary, key=lambda t: t[2])
    loc = rcv_local[worst[0]]
    mid_share = loc["share"][:, [2, 3]].mean(axis=1)  # (zones,)
    dominant = names[int(np.argmax(mid_share))]
    print(f"\n  WORST MATCH (analytic, all bands): {worst[0]} "
          f"(mean |delta| {worst[2]:.3f} s; rendered cross-check {worst[3]:.3f} s, "
          f"mid-band delta {worst[4]:+.3f} s)")
    print(f"  dominant zone in its receiver-weighted absorption: {dominant} "
          f"({mid_share.max()*100:.1f}% of weighted mid-band absorption)")
    for k, n in enumerate(names):
        print(f"    {n:9s} share={mid_share[k]*100:5.1f}%  "
              f"d={loc['distances'][k]:.2f} m  w={loc['weights'][k]:.3f}")
    print("\n  Physical reading: T30 is a late, near-diffuse quantity and")
    print("  barely moves with the receiver; the position effect shows up in")
    print("  EDT and in the spectral balance (E250/E2k). The foam wall sets")
    print("  the fast initial decay, while slow paths that avoid that wall")
    print("  dominate the tail - hence the flat ~1.4 s T30 across bands that")
    print("  neither Sabine variant predicts for the high frequencies.")

    # ---- Plots -------------------------------------------------------------
    plot_geometry(room, source, receivers)
    plot_ir(rendered)
    plot_edc(images, receivers, rendered)
    plot_t60(rows, g, src_local, rcv_local)
    print(f"\nArtifacts written to {os.path.abspath(ART)}")


def plot_geometry(room, source, receivers):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    Lx, Ly, Lz = room.size
    ax.add_patch(plt.Rectangle((0, 0), Lx, Ly, fill=False, lw=1.5,
                               edgecolor="#52514e"))

    def wall(x, color, label, lw=8):
        ax.plot([x, x], [0, Ly], color=color, lw=lw, solid_capstyle="butt",
                label=label, alpha=0.85)

    wall(0, ZONE_COLOR["foam"], "foam wall (x=0)")
    wall(Lx, ZONE_COLOR["concrete"], "concrete wall (x=Lx)")
    ax.plot([0, Lx], [0, 0], color=ZONE_COLOR["wood"], lw=8,
            solid_capstyle="butt", label="wood walls (y faces)")
    ax.plot([0, Lx], [Ly, Ly], color=ZONE_COLOR["wood"], lw=8,
            solid_capstyle="butt")

    ax.scatter(*source[:2], marker="*", s=260, color="#0b0b0b",
               zorder=5, label="source")
    annot_kw = {
        "center": ((12, 10), "left"),
        "hard_corner": ((0, -18), "center"),
        "foam_wall": ((0, 16), "center"),
    }
    for name, p in receivers.items():
        ax.scatter(p[0], p[1], s=110, color=RCV_COLOR[name],
                   edgecolor="white", linewidth=1.2, zorder=5,
                   label=f"receiver: {RCV_LABEL[name]}")
        offset, ha = annot_kw[name]
        ax.annotate(RCV_LABEL[name], (p[0], p[1]),
                    textcoords="offset points", xytext=offset, fontsize=9,
                    ha=ha)
    ax.set_xlim(-0.4, Lx + 0.4)
    ax.set_ylim(-0.7, Ly + 1.0)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Room plan view: material zones, source and receivers "
                 "(floor/ceiling are concrete)")
    ax.set_aspect("equal")
    ax.grid(True, color="#e1e0d9", lw=0.7)
    ax.legend(fontsize=8, loc="lower right", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(os.path.join(ART, "fig1_geometry.png"), dpi=150)
    plt.close(fig)


def plot_ir(rendered):
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for ax, name in zip(axes, rendered):
        _, ir = rendered[name]
        t = np.arange(len(ir)) / FS * 1000
        ax.plot(t, ir, color=RCV_COLOR[name], lw=0.6)
        ax.set_ylabel("amplitude")
        ax.set_title(RCV_LABEL[name], fontsize=10, loc="left")
        ax.grid(True, color="#e1e0d9", lw=0.6)
    axes[-1].set_xlabel("time (ms)")
    axes[0].set_xlim(0, DURATION * 1000)
    fig.suptitle("Broadband impulse responses (sum of six rendered octave "
                 "bands + calibrated late tail)")
    fig.tight_layout()
    fig.savefig(os.path.join(ART, "fig2_impulse_responses.png"), dpi=150)
    plt.close(fig)


def plot_edc(images, receivers, rendered):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    # Left: rendered Schroeder curves, 250 Hz vs 2 kHz, all receivers.
    for name in receivers:
        band_irs, _ = rendered[name]
        for b, ls in ((1, "-"), (4, "--")):
            t = np.arange(len(band_irs[b])) / FS
            axes[0].plot(t * 1000, schroeder_edc(band_irs[b]),
                         color=RCV_COLOR[name], ls=ls, lw=1.4,
                         label=f"{RCV_LABEL[name]}, {int(BAND_FREQS[b])} Hz")
    axes[0].axhspan(-35, -5, color="#eda100", alpha=0.12,
                    label="T30 fit window (-5 to -35 dB)")
    axes[0].set_title("Rendered Schroeder decay curves")
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("EDC (dB)")
    axes[0].set_ylim(-70, 2)
    axes[0].grid(True, color="#e1e0d9", lw=0.6)
    axes[0].legend(fontsize=7, loc="lower left", framealpha=0.95)

    # Right: exact ISM analytic EDC per receiver, mid bands.
    for name, p in receivers.items():
        t, edc = analytic_edc(images, p)
        mid = edc[:, [2, 3]].mean(axis=1)
        axes[1].plot(t * 1000, mid, color=RCV_COLOR[name], lw=1.4,
                     label=RCV_LABEL[name])
    axes[1].axhspan(-35, -5, color="#eda100", alpha=0.12)
    axes[1].set_title("Analytic image-source EDC (500 Hz + 1 kHz average)")
    axes[1].set_xlabel("time (ms)")
    axes[1].grid(True, color="#e1e0d9", lw=0.6)
    axes[1].legend(fontsize=8, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(os.path.join(ART, "fig3_schroeder_edc.png"), dpi=150)
    plt.close(fig)


def plot_t60(rows, g, src_local, rcv_local):
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 9), sharex=True)
    width = 0.16
    x = np.arange(len(BAND_FREQS))
    rcv_names = list(rcv_local)

    for ax, key, title in (
        (axes[0], "t60_measured_t30",
         "T60 (T30 fit, -5 to -35 dB): late decay, nearly position-independent"),
        (axes[1], "edt_measured",
         "EDT (0 to -10 dB fit): early decay, moves with listening position"),
    ):
        for i, name in enumerate(rcv_names):
            r = [row for row in rows if row["receiver"] == name]
            ax.bar(x + (i - 1.5) * width, [v[key] for v in r], width,
                   color=RCV_COLOR[name],
                   label=f"measured, {RCV_LABEL[name]}",
                   edgecolor="white", linewidth=0.8)
        ax.plot(x, g["t60_bands"], "k--", lw=1.6, marker="o", ms=5,
                label="Sabine, global (one constant)")
        ax.plot(x, src_local["t60_bands"], color="#6250d6", lw=1.6,
                marker="s", ms=5,
                label="Sabine, source-distance weighted (requested)")
        for name in rcv_names:
            ax.plot(x, rcv_local[name]["t60_bands"], ":",
                    color=RCV_COLOR[name], lw=1.4,
                    label=f"Sabine, receiver-weighted, {RCV_LABEL[name]}")
        ax.set_ylabel("time (s)")
        ax.set_title(title, fontsize=10, loc="left")
        ax.grid(True, axis="y", color="#e1e0d9", lw=0.7)
    axes[0].legend(fontsize=7.5, ncol=2, framealpha=0.95)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f"{int(f)}" for f in BAND_FREQS])
    axes[1].set_xlabel("octave band centre frequency (Hz)")
    fig.suptitle("Measured vs formula decay times per octave band")
    fig.tight_layout()
    fig.savefig(os.path.join(ART, "fig4_t60_comparison.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
