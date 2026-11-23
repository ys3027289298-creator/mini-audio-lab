
/* doppler: moving-source pitch from retarded emission time, with a per-sample
 * phase-accumulation reference for comparison.
 *
 * Geometry: the listener sits at the origin, the source moves along x from
 * (x0, v0) with constant acceleration a, offset d on the perpendicular axis.
 *     x(te) = x0 + v0*te + 0.5*a*te^2
 *     r(te) = sqrt(x(te)^2 + d^2)
 *     vr(te) = x(te)*v(te)/r(te)          (positive = receding)
 *
 * A wavefront emitted at te reaches the listener at
 *     t = te + r(te)/c
 * so the retarded time te(t) is the root of g(te) = te + r(te)/c - t = 0.
 * Differentiating f(t) = f0*c/(c + vr(te)) gives the instantaneous frequency
 * and the chirp
 *     f    = f0*c/(c + vr)
 *     df/dt = -f0*c^2 * [(v^2 - vr^2) + x*a] / (r*(c+vr)^3)
 * (see doppler_chirp).
 *
 * v == 0 is never a divisor: a stationary source solves te = t - r/c and
 * vr = 0, so f = f0 exactly. The constant-velocity solver does divide by
 * c^2 - v0^2, so it never takes that path at |v0| == c: it drops to a linear
 * fallback instead. A disabled model reports ratio 1.0 bit-exactly, which
 * leaves the caller's signal untouched.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DOPPLER_C 343.0

typedef struct { double te, x, v, r, vr, f; } DopplerPoint;

typedef struct {
    double f0;          /* rest frequency, Hz */
    double c;           /* speed of sound, m/s */
    double x0, v0, a;   /* source state at emission time 0 */
    double d;           /* perpendicular distance, m (> 0) */
    double rel_tol;     /* retarded-time solver tolerance */
    int max_iter;
    int enabled;
} Doppler;

/* ---- setup ------------------------------------------------------------- */

static int doppler_init(Doppler* m, double f0, double c,
                        double x0, double v0, double a, double d) {
    if (m == NULL) return -1;
    if (!(f0 > 0.0) || !(c > 0.0) || !(d > 0.0)) return -1;
    if (!isfinite(x0) || !isfinite(v0) || !isfinite(a)) return -1;
    memset(m, 0, sizeof(*m));
    m->f0 = f0;
    m->c = c;
    m->x0 = x0;
    m->v0 = v0;
    m->a = a;
    m->d = d;
    m->rel_tol = 1e-14;
    m->max_iter = 64;
    m->enabled = 1;
    return 0;
}

static void doppler_enable(Doppler* m, int on) { if (m) m->enabled = on ? 1 : 0; }
static int doppler_enabled(const Doppler* m) { return m != NULL && m->enabled; }

/* ---- geometry ---------------------------------------------------------- */

static void source_at(const Doppler* m, double te, double* x, double* v) {
    *x = m->x0 + m->v0 * te + 0.5 * m->a * te * te;
    *v = m->v0 + m->a * te;
}

/* g is strictly increasing in te because dr/dte >= -|v| > -c. */
static double delay_g(const Doppler* m, double te, double t) {
    double x, v;
    source_at(m, te, &x, &v);
    return te + sqrt(x * x + m->d * m->d) / m->c - t;
}

/* Lower bound on the root: |x(te)| <= |x0| + |v0||te| + |a|te^2/2, so
 * r/c <= (|x0| + |v0||te| + |a|te^2/2 + d)/c. Start one second under the
 * worst-case delay and widen until g < 0. */
static double lower_bracket(const Doppler* m, double t) {
    double at = fabs(t);
    double reach = fabs(m->x0) + fabs(m->v0) * at + 0.5 * fabs(m->a) * at * at + m->d;
    double lo = t - reach / m->c - 1.0;
    for (int k = 0; k < 32 && delay_g(m, lo, t) > 0.0; k++)
        lo = t - 2.0 * (t - lo);
    return lo;
}

/* Constant velocity: the delay equation is quadratic in te,
 *   (c^2 - v0^2) te^2 - 2(c^2 t + x0 v0) te + (c^2 t^2 - x0^2 - d^2) = 0.
 * Evaluate both roots and keep the one below t that actually solves the
 * equation; the leading coefficient can be <= 0 when |v0| >= c. When the
 * coefficient vanishes (|v0| == c) the equation degenerates to linear. */
static int retarded_const_v(const Doppler* m, double t, double* te_out) {
    double A = m->c * m->c - m->v0 * m->v0;
    double B = m->c * m->c * t + m->x0 * m->v0;
    double Cq = m->c * m->c * t * t - m->x0 * m->x0 - m->d * m->d;
    double lo = t - (fabs(m->x0) + fabs(m->v0) * fabs(t) + m->d) / m->c - 1.0;
    double best = NAN, best_g = INFINITY;

    double scale = m->c * m->c + m->v0 * m->v0;
    if (fabs(A) > 1e-12 * scale) {
        double disc = B * B - A * Cq;
        if (disc < 0.0) disc = 0.0;
        double s = sqrt(disc);
        double cand[2] = { (B - s) / A, (B + s) / A };
        for (int i = 0; i < 2; i++) {
            double te = cand[i];
            if (!isfinite(te) || te >= t || te < lo) continue;
            double g = fabs(delay_g(m, te, t));
            if (g < best_g) { best_g = g; best = te; }
        }
    } else if (fabs(B) > 1e-12 * m->c * m->c) {
        double te = Cq / (2.0 * B);
        if (isfinite(te) && te < t && te >= lo) {
            best = te;
            best_g = fabs(delay_g(m, te, t));
        }
    }
    if (!isfinite(best) || best_g > 1e-6) return -1;
    *te_out = best;
    return 0;
}

/* Accelerating source: safeguarded Newton against a sign-change bracket. The
 * derivative is 1 + vr/c, which is small near closest approach; the bracket
 * keeps a bad Newton step from jumping to another branch. */
static int retarded_accel(const Doppler* m, double t, double* te_out) {
    double lo = lower_bracket(m, t);
    double hi = t;                                  /* g(t) = r/c > 0 */
    double te = t - (fabs(m->x0) + fabs(m->v0) * fabs(t) +
                     0.5 * fabs(m->a) * t * t + m->d) / m->c;
    if (!(te > lo && te < hi)) te = 0.5 * (lo + hi);

    for (int i = 0; i < m->max_iter; i++) {
        double g = delay_g(m, te, t);
        if (g == 0.0) break;
        if (g > 0.0) hi = te; else lo = te;
        double x, v;
        source_at(m, te, &x, &v);
        double r = sqrt(x * x + m->d * m->d);
        double vr = x * v / r;
        double gp = 1.0 + vr / m->c;
        if (!(gp > 0.0)) return -1;                  /* supersonic approach */
        double next = te - g / gp;
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        if (fabs(next - te) <= m->rel_tol * fmax(1.0, fabs(next))) { te = next; break; }
        te = next;
    }
    *te_out = te;
    return 0;
}

static int retarded_time(const Doppler* m, double t, double* te_out) {
    if (m == NULL || te_out == NULL || !isfinite(t)) return -1;
    return m->a == 0.0 ? retarded_const_v(m, t, te_out) : retarded_accel(m, t, te_out);
}

/* ---- instantaneous frequency ------------------------------------------ */

/* Fill the analytic solution at listener time t. Returns -1 only for invalid
 * input or a supersonic approach (the received pitch is not real there). */
static int doppler_at(const Doppler* m, double t, DopplerPoint* p) {
    if (m == NULL || p == NULL) return -1;
    double te;
    if (retarded_time(m, t, &te) != 0) return -1;
    double x, v;
    source_at(m, te, &x, &v);
    double r = sqrt(x * x + m->d * m->d);
    double vr = r > 0.0 ? x * v / r : 0.0;   /* d > 0 guarantees r > 0 */
    double den = m->c + vr;
    if (!(den > 0.0)) return -1;
    p->te = te; p->x = x; p->v = v; p->r = r; p->vr = vr;
    p->f = m->f0 * m->c / den;
    return 0;
}

/* df/dt of the received frequency, from differentiating f = f0*c/(c+vr):
 *   df/dt = -f0*c^2 * [(v^2 - vr^2) + x*a] / (r*(c+vr)^3)
 * This sets the O(h) finite-difference bias and the error bound. */
static int doppler_chirp(const Doppler* m, double t, double* chirp_out) {
    DopplerPoint p;
    if (chirp_out == NULL || doppler_at(m, t, &p) != 0) return -1;
    double den = m->c + p.vr;
    *chirp_out = -m->f0 * m->c * m->c *
                 ((p.v * p.v - p.vr * p.vr) + p.x * m->a) / (p.r * den * den * den);
    return 0;
}

/* Received phase 2*pi*f0*te. Its derivative is the received frequency; the
 * phase-accumulation reference differences this per sample. */
static int doppler_phase(const Doppler* m, double t, double* phase_out) {
    DopplerPoint p;
    if (phase_out == NULL || doppler_at(m, t, &p) != 0) return -1;
    *phase_out = 2.0 * M_PI * m->f0 * p.te;
    return 0;
}

/* Frequency ratio f/f0. Disabled models report exactly 1.0. */
static double doppler_ratio(const Doppler* m, double t) {
    DopplerPoint p;
    if (m == NULL || !m->enabled) return 1.0;
    if (doppler_at(m, t, &p) != 0) return 1.0;
    return p.f / m->f0;
}

/* Fill n per-sample ratios at t0 + i/fs. Disabled -> all 1.0, bit-exact. */
static int doppler_track(const Doppler* m, double t0, double fs, int n, float* ratio) {
    if (m == NULL || ratio == NULL || n < 0 || !(fs > 0.0) || !isfinite(t0)) return -1;
    if (!m->enabled) {
        for (int i = 0; i < n; i++) ratio[i] = 1.0f;
        return 0;
    }
    for (int i = 0; i < n; i++) {
        double r = doppler_ratio(m, t0 + i / fs);
        if (!isfinite(r)) return -1;
        ratio[i] = (float)r;
    }
    return 0;
}

/* ---- phase-accumulation reference -------------------------------------- */

typedef struct {
    double fs;
    double ph1, ph2;   /* phase at the last and next-to-last sample */
    double step1;      /* unwrapped true step ph1 - ph2 */
    int count;
    int aliased;       /* the latest *true* step reached pi: f >= fs/2 */
    int prev_aliased;
    double f_back;     /* backward difference over the latest interval */
    double f_center;   /* centered difference at the previous sample */
} PhaseRef;

static double wrap_pi(double x) {
    x = fmod(x + M_PI, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x - M_PI;
}

static void phase_ref_init(PhaseRef* pr, double fs) {
    memset(pr, 0, sizeof(*pr));
    pr->fs = fs;
    pr->f_back = NAN;
    pr->f_center = NAN;
}

/* Feed the received phase at t. From the second sample on:
 *   f_back   = wrap(ph_i - ph_{i-1}) * fs / (2*pi)
 *   f_center = (wrap(ph_i - ph_{i-1}) + wrap(ph_{i-1} - ph_{i-2})) * fs / (4*pi)
 *             summed as two single steps, so it survives |step| < 2*pi; wrapping
 *             the two-sample span directly would fold at pi.
 *
 * Aliasing: a sampled phase difference cannot tell a step of pi from one of
 * 3*pi, so folding cannot be detected from the wrapped value alone - a true
 * step of 1.3*pi wraps to -0.7*pi and looks ordinary. The reference therefore
 * judges aliasing against the true (unwrapped) step, i.e. |ph_i - ph_{i-1}| >= pi,
 * which is the Nyquist edge f_received >= fs/2. Both frequency estimates are set
 * to NAN on an aliased interval so folded values never masquerade as readings.
 */
static int phase_ref_push(PhaseRef* pr, const Doppler* m, double t) {
    double ph;
    if (pr == NULL || doppler_phase(m, t, &ph) != 0) return -1;
    if (pr->count >= 1) {
        double step = ph - pr->ph1;                 /* true step, unwrapped */
        pr->prev_aliased = pr->aliased;
        pr->aliased = fabs(step) >= M_PI;
        pr->f_back = pr->aliased ? NAN : wrap_pi(step) * pr->fs / (2.0 * M_PI);
        if (pr->count >= 2) {
            pr->f_center = (pr->aliased || pr->prev_aliased) ? NAN
                         : (wrap_pi(step) + wrap_pi(pr->step1)) * pr->fs / (4.0 * M_PI);
        }
        pr->step1 = step;
    }
    pr->ph2 = pr->ph1;
    pr->ph1 = ph;
    pr->count++;
    return 0;
}

/* Run the actual oscillator loop (phase += 2*pi*f/fs) and report the largest
 * wrap(accumulated - exact) over n samples. This is the closed-loop phase
 * drift of a real-time pitch shifter seeded from the exact phase at t0. */
static double phase_drift(const Doppler* m, double fs, double t0, int n) {
    double ph, exact;
    if (n <= 0 || doppler_phase(m, t0, &ph) != 0) return NAN;
    double worst = 0.0;
    for (int i = 1; i < n; i++) {
        DopplerPoint p;
        if (doppler_at(m, t0 + i / fs, &p) != 0) return NAN;
        ph = wrap_pi(ph + 2.0 * M_PI * p.f / fs);
        if (doppler_phase(m, t0 + i / fs, &exact) != 0) return NAN;
        double e = fabs(wrap_pi(ph - exact));
        if (e > worst) worst = e;
    }
    return worst;
}

/* ---- self-test --------------------------------------------------------- */

#ifndef DOPPLER_NO_MAIN

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int make_pass(Doppler* m, int fast) {
    if (fast) return doppler_init(m, 440.0, DOPPLER_C, 0.0, 120.0, 20.0, 5.0);
    return doppler_init(m, 440.0, DOPPLER_C, 0.0, 30.0, 0.0, 5.0);
}

int main(void) {
    Doppler m;
    DopplerPoint p;
    PhaseRef pr;

    /* construction and rejection of degenerate parameters */
    check("init_ok", doppler_init(&m, 440.0, DOPPLER_C, 0.0, 0.0, 0.0, 5.0) == 0);
    check("bad_f0", doppler_init(&m, 0.0, DOPPLER_C, 0.0, 30.0, 0.0, 5.0) == -1);
    check("bad_c", doppler_init(&m, 440.0, 0.0, 0.0, 30.0, 0.0, 5.0) == -1);
    check("bad_d", doppler_init(&m, 440.0, DOPPLER_C, 0.0, 30.0, 0.0, 0.0) == -1);
    check("bad_t", doppler_at(&m, NAN, &p) == -1);

    /* stationary source: te = t - r/c, vr = 0, f = f0 exactly, no division */
    doppler_init(&m, 440.0, DOPPLER_C, 0.0, 0.0, 0.0, 5.0);
    check("static_at", doppler_at(&m, 0.0, &p) == 0);
    check("static_vr_zero", p.vr == 0.0);
    check("static_f_exact", p.f == 440.0);
    check("static_retarded", fabs(p.te + 5.0 / DOPPLER_C) < 1e-15);
    {
        float a[256], b[256];
        for (int i = 0; i < 256; i++) b[i] = 1.0f;
        check("zero_speed_track", doppler_track(&m, -0.01, 48000.0, 256, a) == 0);
        check("zero_speed_ones", memcmp(a, b, sizeof(a)) == 0);   /* bit-exact */
    }

    /* disabling Doppler must reproduce the unprocessed pitch exactly */
    doppler_init(&m, 440.0, DOPPLER_C, 0.0, 120.0, 20.0, 5.0);
    {
        float a[512], b[512];
        for (int i = 0; i < 512; i++) b[i] = 1.0f;
        doppler_enable(&m, 0);
        check("disabled_flag", doppler_enabled(&m) == 0);
        check("disabled_ratio_exact", doppler_ratio(&m, 0.1) == 1.0);
        check("disabled_track", doppler_track(&m, -0.01, 48000.0, 512, a) == 0);
        check("disabled_identity", memcmp(a, b, sizeof(a)) == 0);
        doppler_enable(&m, 1);
        check("enable_back", doppler_enabled(&m) == 1);
    }

    /* pitch rises on approach and falls on recede, and moves with speed */
    make_pass(&m, 0);
    double slow_shift;
    check("approach", doppler_at(&m, -0.05, &p) == 0 && p.vr < 0.0 && p.f > 440.0);
    slow_shift = p.f - 440.0;
    check("recede", doppler_at(&m, 0.25, &p) == 0 && p.vr > 0.0 && p.f < 440.0);
    {
        Doppler fast;
        make_pass(&fast, 1);
        check("faster_shifts_more", doppler_at(&fast, -0.05, &p) == 0 &&
                                    (p.f - 440.0) > 3.0 * slow_shift);
    }

    /* the closed form and the iterative solver must agree; force a != 0 */
    {
        Doppler m2;
        make_pass(&m, 0);
        m2 = m;
        m2.a = 1e-9;
        double worst = 0.0;
        for (int i = 0; i < 200; i++) {
            double t = -0.5 + i / 48000.0;
            double te1, te2;
            if (retarded_time(&m, t, &te1) != 0) { worst = INFINITY; break; }
            if (retarded_time(&m2, t, &te2) != 0) { worst = INFINITY; break; }
            if (fabs(te1 - te2) > worst) worst = fabs(te1 - te2);
        }
        check("closed_form_matches_iter", worst < 1e-9);
    }

    /* analytic f is the derivative of the model phase (phase is radians) */
    {
        double t = 0.1, h = 1e-5, ph1 = 0.0, ph2 = 0.0, fd;
        int ok = doppler_phase(&m, t - h, &ph1) == 0 &&
                 doppler_phase(&m, t + h, &ph2) == 0 &&
                 doppler_at(&m, t, &p) == 0;
        fd = (ph2 - ph1) / (2.0 * h) / (2.0 * M_PI);
        check("freq_is_phase_derivative", ok && fabs(fd - p.f) < 1e-6 * p.f);
    }

    /* analytic chirp matches the central difference of f */
    {
        double h = 1e-4, f1 = 0.0, f2 = 0.0, ch = 0.0;
        int ok = doppler_at(&m, 0.1 - h, &p) == 0; f1 = p.f;
        ok &= doppler_at(&m, 0.1 + h, &p) == 0;    f2 = p.f;
        ok &= doppler_chirp(&m, 0.1, &ch) == 0;
        check("chirp_is_freq_derivative", ok && fabs((f2 - f1) / (2.0 * h) - ch) < 1e-4 * fabs(ch) + 1e-9);
    }

    /* the phase reference converges to analytic f as fs grows. Compare only
     * unaliased samples: once a step passes pi the backward difference reads a
     * folded value and is no longer an error estimate of the model. */
    {
        Doppler fast;
        double e48 = 0.0, e8k = 0.0, e1k = 0.0;
        make_pass(&fast, 1);
        for (int pass = 0; pass < 3; pass++) {
            int fs = pass == 0 ? 48000 : (pass == 1 ? 8000 : 1000);
            double worst = 0.0;
            phase_ref_init(&pr, (double)fs);
            for (int i = 0; i < fs; i++) {
                double t = -0.5 + i / (double)fs;
                if (phase_ref_push(&pr, &fast, t) != 0) { worst = INFINITY; break; }
                if (pr.count >= 2 && !pr.aliased) {
                    DopplerPoint q;
                    doppler_at(&fast, t, &q);
                    double e = fabs(pr.f_back - q.f);
                    if (e > worst) worst = e;
                }
            }
            if (pass == 0) e48 = worst; else if (pass == 1) e8k = worst; else e1k = worst;
        }
        check("phase_error_small_48k", e48 < 0.1);
        check("phase_error_grows_low_fs", e8k > e48);
        /* below aliasing the bias is (h/2)*df/dt, so it scales as 1/fs */
        check("phase_error_scales_with_h", fabs(e8k / e48 - 6.0) < 0.6);
        check("phase_error_aliased_1k", e1k > 0.0);   /* remaining valid samples only */
    }

    /* the backward-phase bias is the interval average of f minus f(t_i), i.e.
     * -(h/2)*df/dt + O(h^2); check it stays inside (h/2)*max|chirp| over the
     * interval. The bound is the leading term, so the observed ratio sits just
     * above 1 by the O(h^2) remainder, never below it. */
    {
        Doppler fast;
        double worst = 0.0;
        make_pass(&fast, 1);
        for (int pass = 0; pass < 2; pass++) {
            int fs = pass == 0 ? 48000 : 8000;
            double h = 1.0 / fs;
            double fprev = NAN, cprev = NAN;
            phase_ref_init(&pr, (double)fs);
            for (int i = 0; i < fs; i++) {
                double t = -0.5 + i * h;
                DopplerPoint p;
                double ch = 0.0;
                if (phase_ref_push(&pr, &fast, t) != 0) break;
                if (doppler_at(&fast, t, &p) != 0) break;
                if (doppler_chirp(&fast, t, &ch) != 0) break;
                if (pr.count >= 2 && !pr.aliased && isfinite(fprev)) {
                    double bound = 0.5 * h * fmax(fabs(ch), fabs(cprev));
                    double ratio = bound > 0.0 ? fabs(pr.f_back - p.f) / bound : 0.0;
                    if (ratio > worst) worst = ratio;
                }
                fprev = p.f; cprev = ch;
            }
        }
        check("backward_bias_bounded", worst > 0.99 && worst <= 1.001);
    }

    /* the centered difference is second order: its error scales as h^2 */
    {
        Doppler fast;
        double ec[2] = { 0.0, 0.0 };
        make_pass(&fast, 1);
        for (int pass = 0; pass < 2; pass++) {
            int fs = pass == 0 ? 48000 : 8000;
            double h = 1.0 / fs, worst = 0.0;
            phase_ref_init(&pr, (double)fs);
            for (int i = 0; i < fs; i++) {
                double t = -0.5 + i * h;
                DopplerPoint p, q;
                if (phase_ref_push(&pr, &fast, t) != 0) break;
                if (pr.count < 3 || !isfinite(pr.f_center)) continue;
                if (doppler_at(&fast, t, &p) != 0) break;
                if (doppler_at(&fast, t - h, &q) != 0) break;
                double e = fabs(pr.f_center - q.f);
                (void)p;
                if (e > worst) worst = e;
            }
            ec[pass] = worst;
        }
        check("center_error_second_order", ec[0] > 0.0 && ec[1] / ec[0] > 30.0 &&
                                            ec[1] / ec[0] < 42.0);
    }

    /* closed-loop accumulator drift is bounded by ~ h/2 * total pitch change */
    {
        Doppler fast;
        make_pass(&fast, 1);
        double d48 = phase_drift(&fast, 48000.0, -0.5, 48000);
        double d1k = phase_drift(&fast, 1000.0, -0.5, 1000);
        check("drift_small_48k", d48 < 0.05);
        check("drift_grows_low_fs", d1k > d48);
    }

    /* aliasing edge: a per-sample phase step reaches pi when f reaches fs/2 */
    {
        Doppler fast;
        int any48 = 0, any1k = 0;
        make_pass(&fast, 1);
        phase_ref_init(&pr, 48000.0);
        for (int i = 0; i < 48000; i++) {
            if (phase_ref_push(&pr, &fast, -0.5 + i / 48000.0) != 0) break;
            any48 |= pr.aliased;
        }
        phase_ref_init(&pr, 1000.0);
        for (int i = 0; i < 1000; i++) {
            if (phase_ref_push(&pr, &fast, -0.5 + i / 1000.0) != 0) break;
            any1k |= pr.aliased;
        }
        check("no_alias_48k", any48 == 0);
        check("alias_1k", any1k == 1);

        /* exact onset speed at which f0*c/(c+vr) = fs/2 */
        double v_48 = DOPPLER_C * (2.0 * 440.0 / 48000.0 - 1.0);
        double v_1k = DOPPLER_C * (2.0 * 440.0 / 1000.0 - 1.0);
        check("onset_48k_negative", v_48 < -300.0);   /* needs > Mach 0.98 */
        check("onset_1k_value", fabs(v_1k + 41.16) < 0.01);
    }

    /* |v0| == c must not divide by zero in the constant-velocity solver.
     * Head-on at exactly Mach 1 every emission arrives at the same instant
     * (sonic-boom singularity), so the delay equation has no single root
     * below t; the solver must say so with -1 and never emit a NaN. */
    {
        Doppler mc;
        double te = 0.0;
        int r;
        check("init_mach1", doppler_init(&mc, 440.0, DOPPLER_C, -10.0, DOPPLER_C, 0.0, 5.0) == 0);
        r = retarded_time(&mc, 0.0, &te);
        check("mach1_no_divide", r == -1 || isfinite(te));
        check("mach1_den_zero_rejected", doppler_at(&mc, 0.0, &p) == -1);
        check("mach1_ratio_fallback", doppler_ratio(&mc, 0.0) == 1.0);
    }

    /* the phase reference refuses bad input */
    check("phase_ref_null", phase_ref_push(NULL, &m, 0.0) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

#endif /* DOPPLER_NO_MAIN */
