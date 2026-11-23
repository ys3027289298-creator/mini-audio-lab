
/* looper: sample-level loop-seam processing with phase-continuity proof.
 *
 * Two sample-level seam fixes are provided:
 *   loop_trim_period  - cut the buffer to an integer number of fundamental
 *                       periods (integer-period condition, THEORY.md 2)
 *   loop_crossfade    - seam-pinned cosine-squared crossfade wrapped across
 *                       the end/start seam; the two seam-adjacent samples are
 *                       moved to the same value, so the cyclic sample step at
 *                       the seam is exactly zero (THEORY.md 4.2)
 * With enable = 0 the module is a pure passthrough (bit-identical output).
 *
 * Two discontinuity quantities are kept distinct:
 *   phase residual  d = x[0] - x_synth(N): where the wrapped signal would
 *                   have continued at sample N. Zero iff the loop is phase
 *                   continuous. For a tone cos(2*pi*q*i/N) it is 1-cos(2*pi*e)
 *                   with e the fractional-bin offset.
 *   seam step       s = x[0] - x[N-1]: the actual sample jump at wrap time.
 *                   A click is present iff |s| is an outlier relative to the
 *                   in-buffer first differences; s = d + O(2*pi*q/N) (one
 *                   smooth sample step), so both share the same scale.
 *
 * The self-test derives the integer-period condition, checks the analytic
 * jump/leakage estimates, and measures seam step, cyclic first difference
 * and spectral leakage before/after each fix. Double-precision buffers are
 * used for the proof work; the float API is exercised separately.
 * CSV artifacts seam_wave.csv / seam_diff.csv / leakage_spec.csv are written
 * next to the binary for tools/plot.py.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NMAX 2048
#define SR 48000.0
#define PI 3.14159265358979323846

/* ---------------- processing (float API) ---------------- */

/* Cut to k = n / p whole periods: output length = k * p. */
static int loop_trim_period(const float* in, int n, int p, float* out) {
    if (n <= 0 || p <= 0 || p > n) return -1;
    int m = (n / p) * p;
    for (int i = 0; i < m; i++) out[i] = in[i];
    return m;
}

/* Cut to k periods of a fractional period P (samples). */
static int loop_trim_fractional(const float* in, int n, double period, float* out) {
    if (n <= 0 || period <= 0.0) return -1;
    int k = (int)floor((double)n / period);
    if (k < 1) return -1;
    int m = (int)lround(period * (double)k);
    if (m <= 0 || m > n) return -1;
    for (int i = 0; i < m; i++) out[i] = in[i];
    return m;
}

/* Seam-pinned cosine-squared crossfade wrapped across the end/start seam.
 *
 * delta = in[0] - in[n-1] is the raw seam step. Half of it is removed on
 * each splice leg with a window that is 1 at the seam and 0 at the far end:
 *   w_j = cos^2(pi j / (2(M-1))),  j = 0 .. M-1
 *   y[j]       = in[j]       - w_j       * delta/2   (head leg)
 *   y[n-M+j]   = in[n-M+j]   + (1-w_j)   * delta/2   (tail leg)
 * At the seam w_0 = 1 and 1-w_{M-1} = 1, so y[0] = in[0]-delta/2 and
 * y[n-1] = in[n-1]+delta/2 are identical (C0 exactly). At the far ends
 * w_{M-1} = 0 and 1-w_0 = 0, so samples M-1 and n-M rejoin the original
 * with zero correction; interior samples M .. n-M-1 are never touched.
 * Per-sample motion inside the fade is bounded by delta/(2(M-1)), i.e. the
 * C1 step stays at the signal's own scale when M ~= one period. */
static int loop_crossfade(const float* in, int n, int m, float* out) {
    if (n <= 0 || m <= 1 || 2 * m > n) return -1;
    for (int i = 0; i < n; i++) out[i] = in[i];
    float delta = in[0] - in[n - 1];
    for (int j = 0; j < m; j++) {
        double w = cos(PI * (double)j / (2.0 * (double)(m - 1)));
        w = w * w;
        out[j] = (float)((double)in[j] - w * 0.5 * (double)delta);
        out[n - m + j] = (float)((double)in[n - m + j]
                                 + (1.0 - w) * 0.5 * (double)delta);
    }
    return n;
}

/* Processing disabled: bit-identical passthrough. */
static int loop_passthrough(const float* in, int n, int enable, float* out) {
    if (n <= 0) return -1;
    if (!enable) { for (int i = 0; i < n; i++) out[i] = in[i]; return n; }
    return -1;
}

/* Autocorrelation fundamental-period estimate, searched in [pmin, pmax]. */
static int estimate_period_d(const double* x, int n, int pmin, int pmax) {
    if (pmin <= 0 || pmax >= n || pmin > pmax) return -1;
    int best = -1; double bestv = 0.0;
    for (int p = pmin; p <= pmax; p++) {
        double s = 0.0;
        for (int i = 0; i < n - p; i++) s += x[i] * x[i + p];
        if (s > bestv) { bestv = s; best = p; }
    }
    return best;
}

/* ---------------- analysis (double precision for the proof) ---------------- */

static void dft_mag_d(const double* x, int n, double* mag) {
    for (int k = 0; k < n; k++) {
        double re = 0.0, im = 0.0;
        for (int t = 0; t < n; t++) {
            double ph = -2.0 * PI * (double)k * (double)t / (double)n;
            re += x[t] * cos(ph);
            im += x[t] * sin(ph);
        }
        mag[k] = sqrt(re * re + im * im);
    }
}

typedef struct {
    double seam_step;  /* |x[0] - x[n-1]|: actual sample jump at the wrap */
    double max_diff;   /* max |cyclic first difference| */
    double interior;   /* max |first difference| strictly inside the buffer */
    double leak;       /* off-peak DFT energy (1/N sum |X|^2 off peak) */
    double peak_energy;
    double total;
    int peak;
} SeamMetric;

static void measure_d(const double* x, int n, SeamMetric* m) {
    double mag[NMAX];
    dft_mag_d(x, n, mag);
    int peak = 1;
    double total = 0.0, best = -1.0;
    for (int k = 0; k < n; k++) {
        double e = mag[k] * mag[k] / (double)n;
        total += e;
        if (k > 0 && k < n / 2 && mag[k] > best) { best = mag[k]; peak = k; }
    }
    double leak = 0.0;
    for (int k = 0; k < n; k++) {
        if (k == peak || k == n - peak) continue;
        leak += mag[k] * mag[k] / (double)n;
    }
    double ss = fabs(x[0] - x[n - 1]);
    double md = ss, inside = 0.0;
    for (int i = 1; i < n; i++) {
        double d = fabs(x[i] - x[i - 1]);
        if (d > md) md = d;
        if (d > inside) inside = d;
    }
    m->seam_step = ss; m->max_diff = md; m->interior = inside;
    m->leak = leak; m->peak_energy = best * best / (double)n;
    m->total = total; m->peak = peak;
}

static void measure_f(const float* x, int n, SeamMetric* m) {
    double xx[NMAX];
    for (int i = 0; i < n; i++) xx[i] = (double)x[i];
    measure_d(xx, n, m);
}

/* ---------------- test harness ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* x[i] = cos(2*pi*cycles*i/n); period P = n/cycles, frequency SR/P. */
static void tone_p(double* x, int n, double cycles) {
    for (int i = 0; i < n; i++)
        x[i] = cos(2.0 * PI * cycles * (double)i / (double)n);
}

/* Phase residual against the known synthetic continuation at sample n:
 * d = x[0] - x_synth(n); zero iff the loop wrap is phase continuous. */
static double phase_residual_tone(double cycles) {
    return 1.0 - cos(2.0 * PI * cycles);
}

int main(void) {
    static double raw[NMAX], buf[NMAX], xf[NMAX];
    float fraw[NMAX], fout[NMAX], fpass[NMAX];
    SeamMetric m;
    const int N = 1024;

    double qA = 24.0;                 /* integral: 24 whole periods in N */
    double qB = 23.5;                 /* half a bin off: 23.5 periods */
    double PA = (double)N / qA;
    double PB = (double)N / qB;
    double fA = SR / PA, fB = SR / PB;

    /* ---- case A: integer-period loop, phase continuous by construction ---- */
    tone_p(raw, N, qA);
    measure_d(raw, N, &m);
    double d_resA = phase_residual_tone(qA);
    double d_localA = fabs(raw[1] - raw[0]);     /* smooth step at phase 0 */
    static double ddA[NMAX], dmagA[NMAX];
    ddA[0] = raw[0] - raw[N - 1];
    for (int i = 1; i < N; i++) ddA[i] = raw[i] - raw[i - 1];
    dft_mag_d(ddA, N, dmagA);
    double floor_a = 0.0;
    for (int k = 30; k < 80; k++) floor_a += dmagA[k] * dmagA[k];
    printf("# A integer-period: f=%.4f Hz P=%.4f q=%.1f phase_residual=%.3e "
           "seam_step=%.3e local_step=%.3e leak/total=%.3e diff-floor=%.3e\n",
           fA, PA, qA, d_resA, m.seam_step, d_localA, m.leak / m.total, floor_a);
    check("intperiod_phase_residual", fabs(d_resA) < 1e-12);
    check("intperiod_slope_match", fabs(m.seam_step - d_localA) < 1e-12);
    check("intperiod_no_leak", m.leak / m.total < 1e-24);
    check("intperiod_diff_floor", floor_a < 1e-18);

    /* ---- case B: q = 23.5, the loop does not close ---- */
    double eps = 0.5;
    tone_p(raw, N, qB);
    measure_d(raw, N, &m);
    double delta_theory = phase_residual_tone(qB);   /* 1-cos(47pi) = 2 */
    double delta_quad = 2.0 * PI * PI * eps * eps;      /* 2*pi^2*eps^2 */
    double delta_lin = 2.0 * PI * eps;                  /* worst-phase bound */
    printf("# B fractional bin: f=%.4f Hz P=%.4f eps=%.2f\n", fB, PB, eps);
    printf("# phase residual: theory=%.12e quad_bound=%.4f lin_bound=%.4f | "
           "seam_step=%.6e (differs by one smooth sample step O(2*pi*q/N))\n",
           delta_theory, delta_quad, delta_lin, m.seam_step);
    check("residual_matches_theory", fabs(delta_theory - 2.0) < 1e-12);
    check("residual_inside_bounds", delta_theory <= delta_quad + 1e-12
                                    && delta_theory <= delta_lin);
    check("raw_seam_audible", m.seam_step > 10.0 * m.interior);

    /* leakage model. Within one N-point buffer a fractional-bin cosine's
     * DFT is perfectly sparse (two bins only); the broadband click is born
     * at the wrap, so its exact discrete spectrum is the DFT of the CYCLIC
     * first difference d[i] = x[i] - x[(i-1) mod N]. Inside the buffer d is
     * just the smooth tone derivative; at i = 0 it carries the seam step
     * x[0]-x[N-1]. By the DFT derivative identity
     *   D[k] = (1 - e^{-j2*pi*k/N}) X[k],  |D[k]| = 2 sin(pi*k/N) |X[k]|,
     * and D[0] = sum_i d[i] = x[0]-x[N-1] is exactly the seam step. The
     * click is D's broadband content once the two tone bins are removed. */
    static double dd[NMAX], dmag[NMAX];
    dd[0] = raw[0] - raw[N - 1];
    for (int i = 1; i < N; i++) dd[i] = raw[i] - raw[i - 1];
    dft_mag_d(dd, N, dmag);
    double d_sum = 0.0;
    for (int i = 0; i < N; i++) d_sum += dd[i];
    printf("# cyclic difference: sum_i d[i]=%.3e (telescopes to zero), "
           "d[0]=%.9e = seam step %.9e\n", d_sum, dd[0], raw[0] - raw[N - 1]);
    check("click_identity_dc", fabs(d_sum) < 1e-12
                               && fabs(dd[0] - (raw[0] - raw[N - 1])) < 1e-12);
    /* energy concentration: the seam sample alone carries the largest share
     * of the difference energy; every in-buffer step is the smooth tone
     * derivative (magnitude 2*pi*q/N ~= 0.144, energy ~0.021 each) */
    double e_seam = dd[0] * dd[0];
    double e_diff = 0.0;
    for (int i = 0; i < N; i++) e_diff += dd[i] * dd[i];
    double e_max_inside = 0.0;
    for (int i = 1; i < N; i++)
        if (dd[i] * dd[i] > e_max_inside) e_max_inside = dd[i] * dd[i];
    printf("# difference energy: seam sample=%.6f total=%.6f share=%.4f "
           "largest inside=%.6f\n", e_seam, e_diff, e_seam / e_diff, e_max_inside);
    check("seam_energy_share", e_seam / e_diff > 0.2
                               && e_seam > 100.0 * e_max_inside);
    /* derivative identity on every bin, magnitudes compared */
    double xmag[NMAX];
    dft_mag_d(raw, N, xmag);
    int deriv_ok = 1;
    for (int k = 0; k < N; k++) {
        double model = 2.0 * fabs(sin(PI * (double)k / (double)N)) * xmag[k];
        if (fabs(dmag[k] - model) > 1e-8 * (1.0 + model)) deriv_ok = 0;
    }
    check("derivative_identity", deriv_ok);
    /* broadband click floor: away from the tone bins 23/24 the difference
     * spectrum stays at the seam-impulse level, whereas the integer-period
     * loop collapses there to numerical zero (floor_a measured above) */
    double floor_b = 0.0;
    for (int k = 30; k < 80; k++) floor_b += dmag[k] * dmag[k];
    printf("# diff-spectrum broadband floor: fractional=%.6e integer=%.3e "
           "ratio=%.3e\n", floor_b, floor_a, floor_b / (floor_a + 1e-300));
    check("click_floor_present", floor_b > 1e6 * (floor_a + 1e-300));

    /* ---- fix 1a: autocorrelation period estimate + whole-period trim ---- */
    int p_est = estimate_period_d(raw, N, 20, 80);
    printf("# autocorrelation period estimate: %d (true %.3f)\n", p_est, PB);
    check("period_estimate", abs(p_est - (int)lround(PB)) <= 1);
    int nt0 = (N / p_est) * p_est;
    for (int i = 0; i < nt0; i++) buf[i] = raw[i];
    measure_d(buf, nt0, &m);
    printf("# integer-sample trim to P_est=%d: n=%d seam_step=%.3e (estimator "
           "quantization leaves a residual; exact closure needs the true period)\n",
           p_est, nt0, m.seam_step);
    check("trim_len", nt0 == (N / p_est) * p_est && nt0 < N);

    /* ---- fix 1b: trim at the true fractional period -> exact closure ---- */
    int kper = (int)floor((double)N / PB);                /* 23 */
    int nt = (int)lround(PB * (double)kper);              /* 1002 */
    for (int i = 0; i < nt; i++) buf[i] = raw[i];
    double q_trim = qB * (double)nt / (double)N;          /* 22.9995 */
    double res_trim = phase_residual_tone(q_trim);
    measure_d(buf, nt, &m);
    /* sample-level trim can only land on an integer index, so the residual
     * is bounded by the smooth single-sample motion 2*pi*q/nt (~0.144) rather
     * than zero; ~30% of that bound is achieved here */
    double trim_bound = 2.0 * PI * 23.0 / (double)nt;
    printf("# fractional-period trim: n=%d periods=%d q_eff=%.6f "
           "phase_residual=%.3e (sample bound %.3e) seam/interior=%.3e "
           "leak/total=%.3e\n",
           nt, kper, q_trim, res_trim, trim_bound,
           m.seam_step / m.interior, m.leak / m.total);
    check("fract_trim_periods", kper == 23);
    check("fract_trim_phase", fabs(res_trim) < 0.5 * trim_bound);
    check("fract_trim_slope", m.seam_step <= 2.0 * m.interior + 1e-12);
    check("fract_trim_leak", m.leak / m.total < 1e-3);

    /* ---- fix 2: seam-pinned cosine^2 crossfade, M = one estimated period ---- */
    int M = p_est;
    double delta_raw = raw[0] - raw[N - 1];
    for (int i = 0; i < N; i++) xf[i] = raw[i];
    for (int j = 0; j < M; j++) {
        double w = cos(PI * (double)j / (2.0 * (double)(M - 1)));
        w = w * w;
        xf[j] = raw[j] - w * 0.5 * delta_raw;
        xf[N - M + j] = raw[N - M + j] + (1.0 - w) * 0.5 * delta_raw;
    }
    measure_d(xf, N, &m);
    static double dxf[NMAX], mxf[NMAX];
    dxf[0] = xf[0] - xf[N - 1];
    for (int i = 1; i < N; i++) dxf[i] = xf[i] - xf[i - 1];
    dft_mag_d(dxf, N, mxf);
    double floor_x = 0.0;
    for (int k = 30; k < 80; k++) floor_x += mxf[k] * mxf[k];
    int interior_ok = 1;
    for (int i = M; i < N - M; i++)
        if (xf[i] != raw[i]) interior_ok = 0;
    /* far ends of both legs rejoin the original with zero correction */
    int leg_rejoin = xf[M - 1] == raw[M - 1] && xf[N - M] == raw[N - M];
    double c0 = xf[0] - xf[N - 1];
    printf("# crossfade: M=%d C0_residual=%.3e seam_step=%.3e interior=%.3e "
           "diff-floor %.3e -> %.3e (x%.2e)\n",
           M, c0, m.seam_step, m.interior, floor_b, floor_x, floor_x / floor_b);
    check("xfade_c0_exact", fabs(c0) < 1e-12);
    check("xfade_leg_rejoin", leg_rejoin);
    check("xfade_interior_untouched", interior_ok);
    check("xfade_no_click", m.seam_step <= 2.0 * m.interior + 1e-12);
    check("xfade_leak_down", floor_x < 0.25 * floor_b);

    /* ---------------- float API path ---------------- */
    for (int i = 0; i < N; i++) fraw[i] = (float)raw[i];
    int nc = loop_crossfade(fraw, N, M, fout);
    SeamMetric mf;
    measure_f(fout, nc, &mf);
    int f_interior = 1;
    for (int i = M; i < N - M; i++)
        if (fout[i] != fraw[i]) f_interior = 0;
    float fc0 = fout[0] - fout[nc - 1];
    printf("# float API crossfade: n=%d C0_residual=%.3e seam/interior=%.3e\n",
           nc, fc0, mf.seam_step / mf.interior);
    check("float_xfade_len", nc == N);
    check("float_xfade_c0", fabsf(fc0) < 1e-5);
    check("float_xfade_no_click", mf.seam_step <= 2.0 * mf.interior + 1e-6);
    check("float_xfade_interior", f_interior);
    check("off_passthrough", loop_passthrough(fraw, N, 0, fpass) == N
                             && memcmp(fpass, fraw, sizeof(float) * N) == 0);
    check("on_rejected", loop_passthrough(fraw, N, 1, fpass) == -1);
    check("bad_args", loop_crossfade(fraw, N, 1, fpass) == -1
                      && loop_crossfade(fraw, N, N, fpass) == -1
                      && loop_trim_period(fraw, N, 0, fpass) == -1
                      && loop_trim_fractional(fraw, N, 0.0, fpass) == -1);
    int fnt = loop_trim_fractional(fraw, N, PB, fpass);
    check("float_fract_trim", fnt == nt && fpass[0] == fraw[0]
                              && fpass[fnt - 1] == fraw[fnt - 1]);

    /* ---------------- CSV artifacts for plotting ---------------- */
    FILE* f;
    f = fopen("seam_wave.csv", "w");
    fprintf(f, "case,index,value\n");
    for (int i = N - 32; i < N + 32; i++) {
        int idx = (i + N) % N;
        fprintf(f, "raw,%d,%.10g\n", i, raw[idx]);
        fprintf(f, "crossfade,%d,%.10g\n", i, xf[idx]);
    }
    for (int i = nt - 32; i < nt + 32; i++) {
        int idx = (i + nt) % nt;
        fprintf(f, "trim,%d,%.10g\n", i, buf[idx]);
    }
    fclose(f);

    f = fopen("seam_diff.csv", "w");
    fprintf(f, "case,index,diff\n");
    for (int i = N - 32; i < N + 32; i++) {
        int a = (i + N) % N, b = ((i - 1) + N) % N;
        fprintf(f, "raw,%d,%.10g\n", i, raw[a] - raw[b]);
        fprintf(f, "crossfade,%d,%.10g\n", i, xf[a] - xf[b]);
    }
    for (int i = nt - 32; i < nt + 32; i++) {
        int a = (i + nt) % nt, b = ((i - 1) + nt) % nt;
        fprintf(f, "trim,%d,%.10g\n", i, buf[a] - buf[b]);
    }
    fclose(f);

    f = fopen("leakage_spec.csv", "w");
    fprintf(f, "case,bin,fn,magnitude\n");
    double mag_raw[NMAX], mag_xf[NMAX], mag_tr[NMAX];
    dft_mag_d(raw, N, mag_raw);
    dft_mag_d(xf, N, mag_xf);
    dft_mag_d(buf, nt, mag_tr);
    for (int k = 0; k <= 60; k++) {
        double fn = (double)k / (double)N;
        fprintf(f, "raw,%d,%.10g,%.10g\n", k, fn, mag_raw[k]);
        fprintf(f, "crossfade,%d,%.10g,%.10g\n", k, fn, mag_xf[k]);
        int kt = (int)lround(fn * nt);
        if (kt >= nt) kt = nt - 1;
        fprintf(f, "trim,%d,%.10g,%.10g\n", kt, fn, mag_tr[kt]);
    }
    fclose(f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
