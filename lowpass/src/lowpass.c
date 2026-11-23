
/* lowpass: 4th-order Butterworth lowpass (two biquads, bilinear transform).
 *
 * Sweeps the cutoff frequency and scores two objectives to minimize:
 *   A = aliasing residual (dB)  — stop-band energy above fs/4 that folds
 *                                 back under 2:1 decimation, from the measured
 *                                 frequency response;
 *   P = phase distortion (rad)  — magnitude-weighted RMS residual of the
 *                                 pass-band phase (0..125 Hz, fs=1000) from its
 *                                 best-fit straight line, i.e. deviation from
 *                                 pure delay.
 * Non-dominated candidates form the Pareto front; three non-dominated
 * points are evaluated through a linear-resampler cascade (algorithm copied
 * from ../resampler/src/resampler.c).
 */
#include <stdio.h>
#include <math.h>

#define PI       3.14159265358979323846
#define FS       1000.0          /* sample rate used by every measurement */
#define NFFT     1024            /* radix-2, one block = 1.024 s */
#define NBIN     (NFFT / 2)      /* one-sided bins, 0..512 (0..500 Hz) */
#define NSEC     2               /* biquad sections -> 4th order */
#define NCND     15              /* cutoff sweep candidates */
#define F_SIG    125.0           /* signal band edge, kept below every cutoff */
#define F_ALIAS  250.0           /* decimation Nyquist; energy above folds */

typedef struct { float b0, b1, b2, a1, a2, z1, z2; } Biquad;

/* double-precision section used for design-time response measurement, so
   the objective numbers do not pick up float32 coefficient quantization */
typedef struct { double b0, b1, b2, a1, a2, z1, z2; } BiquadD;

typedef struct {
    double fc;                   /* cutoff as fraction of fs */
    double alias_db;             /* objective A */
    double phase_rad;            /* objective P */
    double gd_spread_ms;         /* corroborating group-delay spread in band */
    int    dominated;            /* Pareto flag */
} Candidate;

/* ----- filter design: analog Butterworth poles -> digital biquads -------- */

static void design_section_d(BiquadD* q, double fc_frac, int k) {
    /* analog 4th-order Butterworth pole on the left half-plane,
       theta = pi/2 + pi*(2k+1)/(2*4), k = 0..1 for the two conjugate pairs */
    double theta = PI / 2.0 + PI * (2.0 * k + 1.0) / 8.0;
    double re = cos(theta);
    double wd = 2.0 * PI * fc_frac * FS;           /* cutoff, rad/s */
    double W = 2.0 * FS * tan(wd / (2.0 * FS));    /* frequency prewarp */
    double W2 = W * W;
    /* bilinear transform of (s^2 - 2*re*W*s + W^2) denominator,
       s = 2*fs*(z-1)/(z+1) */
    double c = 2.0 * FS;
    double den = W2 - 2.0 * re * W * c + c * c;
    double a1 = 2.0 * (W2 - c * c) / den;
    double a2 = (W2 + 2.0 * re * W * c + c * c) / den;
    double g  = W2 / den;                          /* lowpass numerator gain */
    q->b0 = g;
    q->b1 = 2.0 * g;
    q->b2 = g;
    q->a1 = a1;
    q->a2 = a2;
    q->z1 = q->z2 = 0.0;
}

/* runtime path stays float, matching the data type used by the other
   modules in this repo */
static void design_section(Biquad* q, double fc_frac, int k) {
    BiquadD d;
    design_section_d(&d, fc_frac, k);
    q->b0 = (float)d.b0; q->b1 = (float)d.b1; q->b2 = (float)d.b2;
    q->a1 = (float)d.a1; q->a2 = (float)d.a2;
    q->z1 = q->z2 = 0.0f;
}

static double biquad_process_d(BiquadD* q, double x) {
    double y = q->b0 * x + q->z1;
    q->z1 = q->b1 * x - q->a1 * y + q->z2;
    q->z2 = q->b2 * x - q->a2 * y;
    return y;
}

static float biquad_process(Biquad* q, float x) {
    /* transposed direct form II */
    float y = q->b0 * x + q->z1;
    q->z1 = q->b1 * x - q->a1 * y + q->z2;
    q->z2 = q->b2 * x - q->a2 * y;
    return y;
}

static void filter_block(Biquad secs[NSEC], const float* in, float* out, int n) {
    for (int i = 0; i < n; i++) {
        float v = in[i];
        for (int s = 0; s < NSEC; s++) v = biquad_process(&secs[s], v);
        out[i] = v;
    }
}

/* ----- radix-2 FFT (same iterative structure as ../fft/src/fft.c) -------- */

typedef struct { double re, im; } C;

static void fft(C* a, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = -2.0 * PI / len;
        double wr = cos(ang), wi = sin(ang);
        for (int i = 0; i < n; i += len) {
            double cwr = 1.0, cwi = 0.0;
            for (int k = 0; k < len / 2; k++) {
                C u = a[i + k];
                double vr = a[i + k + len / 2].re * cwr - a[i + k + len / 2].im * cwi;
                double vi = a[i + k + len / 2].re * cwi + a[i + k + len / 2].im * cwr;
                a[i + k].re = u.re + vr;
                a[i + k].im = u.im + vi;
                a[i + k + len / 2].re = u.re - vr;
                a[i + k + len / 2].im = u.im - vi;
                double nr = cwr * wr - cwi * wi;
                cwi = cwr * wi + cwi * wr;
                cwr = nr;
            }
        }
    }
}

/* ----- response measurement ---------------------------------------------- */

static void measure_response(double fc_frac, double* alias_db,
                             double* phase_rad, double* gd_spread_ms) {
    BiquadD secs[NSEC];
    for (int k = 0; k < NSEC; k++) design_section_d(&secs[k], fc_frac, k);

    /* impulse response, trimmed past the transient */
    double h[NFFT];
    for (int i = 0; i < NFFT; i++) {
        double v = i == 0 ? 1.0 : 0.0;
        for (int s = 0; s < NSEC; s++) v = biquad_process_d(&secs[s], v);
        h[i] = v;
    }
    C H[NFFT];
    for (int i = 0; i < NFFT; i++) { H[i].re = h[i]; H[i].im = 0.0; }
    fft(H, NFFT);

    double mag[NBIN + 1], phase[NBIN + 1];
    for (int b = 0; b <= NBIN; b++) {
        mag[b]   = sqrt(H[b].re * H[b].re + H[b].im * H[b].im);
        phase[b] = atan2(H[b].im, H[b].re);
    }

    /* objective A: energy above the decimation Nyquist folds into band.
       sum |H|^2 over bins 251..511, normalized by DC gain -> dB */
    double e_alias = 0.0;
    for (int b = (int)(F_ALIAS / FS * NFFT) + 1; b < NBIN; b++)
        e_alias += mag[b] * mag[b];
    *alias_db = 10.0 * log10(e_alias / (mag[0] * mag[0]) + 1e-20);

    /* objective P: weighted RMS of unwrapped phase minus best-fit line over
       the fixed signal band 0..125 Hz. weights = |H|^2 so attenuated bins do
       not dominate the residual. */
    int b_hi = (int)(F_SIG / FS * NFFT);          /* bin 128 */
    double unw[NBIN + 1];
    unw[0] = phase[0];
    for (int b = 1; b <= b_hi; b++) {
        double d = phase[b] - phase[b - 1];
        while (d >  PI) d -= 2.0 * PI;
        while (d < -PI) d += 2.0 * PI;
        unw[b] = unw[b - 1] + d;
    }
    double sw = 0.0, sf = 0.0, sp = 0.0, sff = 0.0, sfp = 0.0;
    for (int b = 0; b <= b_hi; b++) {
        double w = mag[b] * mag[b];
        double f = (double)b;                      /* bin index as x */
        sw += w; sf += w * f; sp += w * unw[b];
        sff += w * f * f; sfp += w * f * unw[b];
    }
    double slope = (sw * sfp - sf * sp) / (sw * sff - sf * sf);
    double p0    = (sp - slope * sf) / sw;
    double acc = 0.0;
    for (int b = 0; b <= b_hi; b++) {
        double w = mag[b] * mag[b];
        double r = unw[b] - (p0 + slope * b);
        acc += w * r * r;
    }
    *phase_rad = sqrt(acc / sw);
    /* group delay per bin = -d(phase)/d(omega), domega per bin = 2*pi/NFFT */
    double gmin = 1e9, gmax = -1e9;
    for (int b = 1; b <= b_hi; b++) {
        double dphi = unw[b] - unw[b - 1];
        double gd = -dphi * NFFT / (2.0 * PI);    /* samples */
        if (gd < gmin) gmin = gd;
        if (gd > gmax) gmax = gd;
    }
    *gd_spread_ms = (gmax - gmin) * 1000.0 / FS;
}

/* ----- linear resampler, copied verbatim in logic from
         ../resampler/src/resampler.c (module-local copy) ------------------ */

static int resample_linear(const float* in, int n_in, double ratio,
                           float* out, int cap) {
    if (n_in <= 0 || ratio <= 0.0) return -1;
    int n_out = (int)floor((double)n_in / ratio);
    if (n_out <= 0 || n_out > cap) return -1;
    for (int i = 0; i < n_out; i++) {
        double src = (double)i * ratio;
        int i0 = (int)floor(src);
        double t = src - (double)i0;
        int i1 = i0 + 1 < n_in ? i0 + 1 : n_in - 1;
        out[i] = (float)((1.0 - t) * in[i0] + t * in[i1]);
    }
    return n_out;
}

/* ----- cascade test on a multitone signal --------------------------------
 *
 * Input block NIN = 2048 at fs = 1000 Hz; after 2:1 decimation the output
 * has NOUT = 1024 at fs' = 500 Hz. Both FFTs then share the same bin spacing
 * (1000/2048 = 500/1024 = 0.4883 Hz), so every test tone lands on an exact
 * bin in both domains and no window is needed. Tones below 250 Hz stay on
 * the same bin; tones above 250 Hz fold onto bin NOUT-b, the standard
 * decimation image at fs' - f.
 */

#define NIN  2048
#define NOUT (NIN / 2)
#define DF   (FS / NIN)                 /* 0.4883 Hz per bin, both FFTs */

/* Test bins, chosen so no tone mirrors or collides:
   pass bins 100/200/400  -> 48.8/97.7/195.3 Hz, all below the 250 Hz
   decimation Nyquist, so they stay on the same output bin;
   alias bins 700/750/850/950 -> 342..464 Hz (between 250 and 500 Hz),
   folding onto output bins NOUT-b = 324/274/174/74, disjoint from the
   pass bins. */
#define BIN_PASS1 100
#define BIN_PASS2 200
#define BIN_PASS3 400
#define BIN_AL1   700
#define BIN_AL2   750
#define BIN_AL3   850
#define BIN_AL4   950

static double band_energy_out(const C* Z, int b) {
    return Z[b].re * Z[b].re + Z[b].im * Z[b].im;
}

static void cascade_measure(double fc_frac, double* alias_ratio_db,
                            double* pass_gain_db) {
    static const int bp[3] = { BIN_PASS1, BIN_PASS2, BIN_PASS3 };
    static const int ba[4] = { BIN_AL1, BIN_AL2, BIN_AL3, BIN_AL4 };
    float x[NIN], y[NIN];
    for (int i = 0; i < NIN; i++) {
        double v = 0.0;
        for (int k = 0; k < 3; k++) v += cos(2.0 * PI * bp[k] * DF * i / FS);
        for (int k = 0; k < 4; k++) v += cos(2.0 * PI * ba[k] * DF * i / FS);
        x[i] = (float)v;
    }
    Biquad secs[NSEC];
    for (int k = 0; k < NSEC; k++) design_section(&secs[k], fc_frac, k);
    filter_block(secs, x, y, NIN);

    float z[NOUT];
    int nz = resample_linear(y, NIN, 2.0, z, NIN);   /* 2:1 downsample */

    C Z[NOUT];
    int nfft_out = NOUT;
    for (int i = 0; i < nfft_out; i++) { Z[i].re = 0.0; Z[i].im = 0.0; }
    for (int i = 0; i < nz; i++) Z[i].re = z[i];
    fft(Z, nfft_out);

    double e_in = 0.0;
    for (int k = 0; k < 3; k++) e_in += band_energy_out(Z, bp[k]);
    double e_fold = 0.0;
    for (int k = 0; k < 4; k++) e_fold += band_energy_out(Z, NOUT - ba[k]);
    *alias_ratio_db = 10.0 * log10(e_fold / e_in + 1e-20);

    /* pass-band gain: mean |H| at the in-band test bins, filter response
       measured on a 2048-point impulse so its bin spacing matches */
    BiquadD s2[NSEC];
    for (int k = 0; k < NSEC; k++) design_section_d(&s2[k], fc_frac, k);
    double h[NIN];
    for (int i = 0; i < NIN; i++) {
        double v = i == 0 ? 1.0 : 0.0;
        for (int s = 0; s < NSEC; s++) v = biquad_process_d(&s2[s], v);
        h[i] = v;
    }
    C H[NIN];
    for (int i = 0; i < NIN; i++) { H[i].re = h[i]; H[i].im = 0.0; }
    fft(H, NIN);
    double gsum = 0.0;
    for (int k = 0; k < 3; k++)
        gsum += sqrt(H[bp[k]].re * H[bp[k]].re + H[bp[k]].im * H[bp[k]].im);
    *pass_gain_db = 20.0 * log10(gsum / 3.0);
}

/* ----- Pareto ranking ---------------------------------------------------- */

static int is_dominated(const Candidate* c, const Candidate* all, int n) {
    for (int j = 0; j < n; j++) {
        if (all[j].fc == c->fc) continue;
        if (all[j].alias_db  <= c->alias_db &&
            all[j].phase_rad <= c->phase_rad &&
            (all[j].alias_db  <  c->alias_db ||
             all[j].phase_rad <  c->phase_rad))
            return 1;
    }
    return 0;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* cutoff sweep: 0.02..0.30 of fs, 15 points, every one keeps the
       125 Hz signal band below cutoff */
    Candidate c[NCND];
    for (int i = 0; i < NCND; i++) {
        c[i].fc = 0.02 + 0.28 * i / (NCND - 1);
        measure_response(c[i].fc, &c[i].alias_db, &c[i].phase_rad,
                         &c[i].gd_spread_ms);
        c[i].dominated = 0;
    }
    int n_front = 0;
    for (int i = 0; i < NCND; i++) {
        c[i].dominated = is_dominated(&c[i], c, NCND);
        if (!c[i].dominated) n_front++;
    }

    FILE* f = fopen("sweep.csv", "w");
    fprintf(f, "fc_frac,fc_hz,alias_residual_db,phase_distortion_rad,"
               "groupdelay_spread_ms,pareto\n");
    for (int i = 0; i < NCND; i++)
        fprintf(f, "%.4f,%.1f,%.8f,%.10f,%.6f,%d\n", c[i].fc, c[i].fc * FS,
                c[i].alias_db, c[i].phase_rad, c[i].gd_spread_ms,
                c[i].dominated ? 0 : 1);
    fclose(f);

    /* front ordered by objective A ascending; check no dominated point
       slipped onto the front */
    int front[NCND], m = 0;
    for (int i = 0; i < NCND; i++) if (!c[i].dominated) front[m++] = i;
    for (int i = 0; i < m; i++)
        for (int j = i + 1; j < m; j++)
            if (c[front[j]].alias_db < c[front[i]].alias_db) {
                int t = front[i]; front[i] = front[j]; front[j] = t;
            }
    int front_clean = 1;
    for (int i = 0; i < m; i++)
        if (is_dominated(&c[front[i]], c, NCND)) front_clean = 0;
    /* along a clean 2-objective front the other objective must descend as A
       worsens (fc rises) */
    for (int i = 1; i < m; i++)
        if (!(c[front[i]].phase_rad < c[front[i - 1]].phase_rad))
            front_clean = 0;

    /* three non-dominated points: both ends and the middle of the front */
    int pick[3];
    pick[0] = front[0];
    pick[1] = front[m / 2];
    pick[2] = front[m - 1];
    int picks_distinct = pick[0] != pick[1] && pick[1] != pick[2];

    FILE* g = fopen("cascade.csv", "w");
    fprintf(g, "label,fc_frac,fc_hz,alias_residual_db,phase_distortion_rad,"
               "cascade_fold_vs_pass_db,passband_gain_db\n");
    double baseline_alias;
    /* no-filter baseline: multitone straight into the 2:1 resampler */
    {
        static const int bp[3] = { BIN_PASS1, BIN_PASS2, BIN_PASS3 };
        static const int ba[4] = { BIN_AL1, BIN_AL2, BIN_AL3, BIN_AL4 };
        float x[NIN], z[NOUT];
        for (int i = 0; i < NIN; i++) {
            double v = 0.0;
            for (int k = 0; k < 3; k++) v += cos(2.0 * PI * bp[k] * DF * i / FS);
            for (int k = 0; k < 4; k++) v += cos(2.0 * PI * ba[k] * DF * i / FS);
            x[i] = (float)v;
        }
        int nz = resample_linear(x, NIN, 2.0, z, NIN);
        C Z[NOUT];
        for (int i = 0; i < NOUT; i++) { Z[i].re = 0.0; Z[i].im = 0.0; }
        for (int i = 0; i < nz; i++) Z[i].re = z[i];
        fft(Z, NOUT);
        double e_in = 0.0, e_fold = 0.0;
        for (int k = 0; k < 3; k++) e_in += band_energy_out(Z, bp[k]);
        for (int k = 0; k < 4; k++) e_fold += band_energy_out(Z, NOUT - ba[k]);
        baseline_alias = 10.0 * log10(e_fold / e_in + 1e-20);
    }
    fprintf(g, "no_filter_baseline,0.0000,0.0,0.0000,0.00000,%.4f,0.0000\n",
            baseline_alias);

    static const char* labels[3] = { "P1_min_alias", "P2_balanced", "P3_min_phase" };
    for (int p = 0; p < 3; p++) {
        double fold_db, gain_db;
        cascade_measure(c[pick[p]].fc, &fold_db, &gain_db);
        fprintf(g, "%s,%.4f,%.1f,%.8f,%.10f,%.6f,%.6f\n", labels[p],
                c[pick[p]].fc, c[pick[p]].fc * FS,
                c[pick[p]].alias_db, c[pick[p]].phase_rad, fold_db, gain_db);
    }
    fclose(g);

    /* ----- self tests ----- */
    printf("candidates=%d front_size=%d front_clean=%d\n",
           NCND, n_front, front_clean);
    check("front_has_three_plus", n_front >= 3);
    check("three_picks_non_dominated",
          !c[pick[0]].dominated && !c[pick[1]].dominated &&
          !c[pick[2]].dominated);
    check("three_picks_distinct", picks_distinct);
    {
        /* a dominated candidate must be genuinely dominated by some other
           point (the 100 Hz cutoff here is Pareto-dominated by 80 Hz), and
           the sweep must show the expected direction of the tradeoff */
        int dom_ok = 1;
        for (int i = 0; i < NCND; i++)
            if (c[i].dominated && !is_dominated(&c[i], c, NCND)) dom_ok = 0;
        check("dominated_flags_correct", dom_ok);
        int trend_ok = c[NCND - 1].alias_db > c[0].alias_db &&
                       c[NCND - 1].phase_rad < c[0].phase_rad;
        check("tradeoff_direction", trend_ok);
    }
    check("front_is_pareto", front_clean);
    {
        /* DC must pass at unity, high cutoff should be near 0 dB */
        double a, p, g2;
        measure_response(0.30, &a, &p, &g2);
        Biquad s[NSEC];
        for (int k = 0; k < NSEC; k++) design_section(&s[k], 0.30, k);
        float acc = 0.0f;
        for (int i = 0; i < 4096; i++) {
            float v = 1.0f;
            for (int q = 0; q < NSEC; q++) v = biquad_process(&s[q], v);
            acc = v;
        }
        check("dc_unity_gain", fabs(acc - 1.0f) < 1e-3f);
    }
    {
        /* the filter must improve the cascade over the no-filter baseline */
        double fold_db, gain_db;
        cascade_measure(c[pick[2]].fc, &fold_db, &gain_db);
        check("cascade_beats_baseline", fold_db < baseline_alias);
    }
    {
        /* resampler length derivation, same cases as the resampler module */
        float tiny[64], tout[256];
        for (int i = 0; i < 64; i++) tiny[i] = (float)i;
        check("resampler_down_len",
              resample_linear(tiny, 64, 2.0, tout, 256) == 32);
        check("resampler_bad_ratio",
              resample_linear(tiny, 64, 0.0, tout, 256) == -1);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
