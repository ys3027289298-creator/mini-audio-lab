#!/usr/bin/env python3
"""Window amplitude compensation analysis.

Derives and validates:
  (1) amplitude deviation g(delta)-1 as a function of window length N and
      tone frequency f (via delta = N*f - round(N*f));
  (2) an error bound of the closed-form g near the bin boundary (delta=0.5);
and uses both to compute the usable frequency band and compensation
coefficients c(delta)=1/g(delta) per window. Residual > 1e-3 => incomplete.
"""
import csv
import os
import numpy as np

FS = 48000.0
TOL = 1e-3
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# Cosine-sum windows: w[n] = sum_m a[m] cos(2 pi m n / N), periodic form.
WINDOWS = {
    "rect":     [1.0],
    "hann":     [0.5, -0.5],
    "hamming":  [0.54, -0.46],
    "blackman": [0.42, -0.5, 0.08],
}
NS = [64, 256, 1024, 4096]


def window_samples(a, N):
    n = np.arange(N)
    w = np.zeros(N)
    for m, am in enumerate(a):
        w += am * np.cos(2.0 * np.pi * m * n / N)
    return w


def dirichlet(nu, N):
    """D_N(2 pi nu / N) = sum_{n=0}^{N-1} exp(-j 2 pi nu n / N), complex."""
    nu = np.asarray(nu, float)
    out = np.full(nu.shape, complex(N, 0.0))
    nz = np.abs(np.sin(np.pi * nu / N)) > 1e-15
    v = nu[nz]
    out[nz] = (np.sin(np.pi * v) / np.sin(np.pi * v / N)
               * np.exp(-1j * np.pi * v * (N - 1.0) / N))
    return out


def g_exact(a, delta, N):
    """|W(2 pi delta / N)| / W(0): exact discrete-window gain."""
    delta = np.asarray(delta, float)
    W = np.zeros(delta.shape, complex)
    for m, am in enumerate(a):
        W += 0.5 * am * (dirichlet(delta - m, N) + dirichlet(delta + m, N))
    return np.abs(W) / (N * a[0])


def g_approx(a, delta):
    """Closed-form (large-N sinc) gain; the compensation formula."""
    delta = np.asarray(delta, float)
    W = np.zeros(delta.shape)
    for m, am in enumerate(a):
        W += 0.5 * am * ((-1.0) ** m) * (np.sinc(delta - m) + np.sinc(delta + m))
    return np.abs(W) / a[0]


def err_bound(a, delta, N):
    """Relative error bound of g_approx near bin offset delta:
    eps <= pi (delta + sum_{m>=1} m|a_m|/a_0) / (6 N^2 g(delta))."""
    delta = np.asarray(delta, float)
    s = sum(m * abs(am) for m, am in enumerate(a)) / a[0]
    return np.pi * (delta + s) / (6.0 * N * N * g_approx(a, delta))


def g_measured(a, N, k0, delta):
    """FFT of a windowed complex tone at (k0+delta)/N cycles/sample."""
    n = np.arange(N)
    w = window_samples(a, N)
    x = np.exp(2j * np.pi * (k0 + delta) * n / N)
    X = np.fft.fft(x * w)
    return abs(X[k0]) / w.sum()


def residual_real_measured(a, N, k0, delta):
    """Compensated amplitude residual for a real cosine (worst case:
    includes the negative-frequency image term)."""
    n = np.arange(N)
    w = window_samples(a, N)
    x = np.cos(2.0 * np.pi * (k0 + delta) * n / N)
    X = np.fft.fft(x * w)
    amp = 2.0 * abs(X[k0]) / w.sum() / g_approx(a, np.array([delta]))[0]
    return abs(amp - 1.0)


def main():
    os.makedirs(OUT, exist_ok=True)
    deltas = np.linspace(0.0, 0.5, 51)

    # ---- (1) deviation expression: exact vs closed form vs measured ----
    print("== deviation expression validation (N=256, k0=32, complex tone) ==")
    N1, K1 = 256, 32
    rows = []
    for name, a in WINDOWS.items():
        ge = g_exact(a, deltas, N1)
        ga = g_approx(a, deltas)
        gm = np.array([g_measured(a, N1, K1, d) for d in deltas])
        e_meas = np.max(np.abs(gm - ge) / ge)
        e_form = np.max(np.abs(ga - ge) / ge)
        print(f"  {name:8s} max|meas-exact|={e_meas:.2e}  "
              f"max|closedform-exact|={e_form:.2e}")
        for d, ve, va, vm in zip(deltas, ge, ga, gm):
            rows.append([name, f"{d:.3f}", f"{ve:.8f}", f"{va:.8f}",
                         f"{vm:.8f}", f"{20*np.log10(ve):.4f}"])
    with open(os.path.join(OUT, "deviation_curve.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["window", "delta_bins", "g_exact", "g_closedform",
                    "g_measured", "deviation_dB"])
        w.writerows(rows)

    # ---- (2) error bound near bin boundary ----
    print("\n== closed-form error vs bound, max over delta in [0,0.5] ==")
    print(f"  {'window':8s} {'N':>5s} {'max_rel_err':>12s} {'bound':>12s}  ok")
    bound_ok = True
    for name, a in WINDOWS.items():
        for N in NS:
            ge = g_exact(a, deltas, N)
            err = np.max(np.abs(g_approx(a, deltas) - ge) / ge)
            bnd = np.max(err_bound(a, deltas, N))
            ok = err <= bnd
            bound_ok &= ok
            print(f"  {name:8s} {N:5d} {err:12.3e} {bnd:12.3e}  {ok}")

    # ---- usable band: closed-form error + real-signal image leakage ----
    print("\n== usable band (compensated residual < 1e-3, real signal) ==")
    dg = np.linspace(0.0, 0.5, 21)
    band_rows = []
    worst_in_band = 0.0
    for name, a in WINDOWS.items():
        for N in NS:
            eps = np.abs(g_approx(a, dg) - g_exact(a, dg, N)) / g_exact(a, dg, N)
            k_lo = k_hi = None
            for k0 in range(1, N // 2):
                img = g_exact(a, 2.0 * k0 + dg, N) / g_exact(a, dg, N)
                if np.max(eps + img) < TOL:
                    if k_lo is None:
                        k_lo = k0
                    k_hi = k0
            if k_lo is None:
                print(f"  {name:8s} N={N:5d}: no usable band")
                band_rows.append([name, N, "", "", "", "", ""])
                continue
            # verify with real FFTs at band edges and interior, worst delta
            checks = [(k_lo, 0.5), (k_lo, 0.0), (k_hi, 0.5),
                      ((k_lo + k_hi) // 2, 0.31), ((k_lo + k_hi) // 2, 0.5)]
            rmax = max(residual_real_measured(a, N, k, d) for k, d in checks)
            worst_in_band = max(worst_in_band, rmax)
            flo, fhi = k_lo / N * FS, k_hi / N * FS
            print(f"  {name:8s} N={N:5d}: bins [{k_lo:5d},{k_hi:5d}]  "
                  f"[{flo:9.1f},{fhi:9.1f}] Hz  max_residual={rmax:.2e}")
            band_rows.append([name, N, k_lo, k_hi, f"{flo:.1f}",
                              f"{fhi:.1f}", f"{rmax:.3e}"])
    with open(os.path.join(OUT, "usable_band.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["window", "N", "k_lo", "k_hi", "f_lo_Hz@48k",
                    "f_hi_Hz@48k", "max_residual"])
        w.writerows(band_rows)

    # ---- failure point: formula breaks near DC for real signals ----
    print("\n== failure point (formula invalid) ==")
    a = WINDOWS["hann"]
    Nf, kf, df = 1024, 1, 0.5
    r = residual_real_measured(a, Nf, kf, df)
    pred = g_exact(a, np.array([2 * kf + df]), Nf)[0] / g_exact(a, np.array([df]), Nf)[0]
    print(f"  hann N={Nf} f0={kf+df} bins ({(kf+df)/Nf*FS:.1f} Hz @48k): "
          f"residual={r:.3e} (predicted image term {pred:.3e}) >> 1e-3")

    # ---- deviation curve plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for name, a in WINDOWS.items():
        ge = g_exact(a, deltas, 256)
        ax.plot(deltas, 20 * np.log10(ge), label=name)
        gm = [g_measured(a, 256, 32, d) for d in deltas[::10]]
        ax.plot(deltas[::10], 20 * np.log10(gm), ".", color=ax.lines[-1].get_color())
    ax.set_xlabel("bin offset $\\delta$ (bins)")
    ax.set_ylabel("amplitude deviation $20\\log_{10} g(\\delta)$ (dB)")
    ax.set_title("window scalloping loss, N=256 (lines=exact, dots=measured)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "deviation_curve.png"), dpi=120)

    print(f"\nbound_ok={bound_ok}  worst_in_band_residual={worst_in_band:.3e}  "
          f"tol={TOL}")
    print(f"complete={bound_ok and worst_in_band < TOL}")


if __name__ == "__main__":
    main()
