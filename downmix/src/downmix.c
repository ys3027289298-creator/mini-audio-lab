
/* downmix: energy-conserving 5.1 -> stereo downmix.
 *
 * y = M x, with the per-column unit-norm constraint wL_j^2 + wR_j^2 = 1
 * so that each input channel keeps its energy under the uncorrelated-channel
 * assumption. The LFE column is in-phase (same sign in Lt and Rt): an
 * opposite-polarity LFE keeps its stereo energy but vanishes in the mono
 * fold (Lt+Rt)/sqrt(2), which is the "low frequencies cancel" defect.
 *
 * The disabled context is a pure passthrough and produces byte-identical
 * output (memcpy, no arithmetic).
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define CH 6
#define N  256
#define PI 3.1415926535897932

enum { L_CH = 0, R_CH = 1, C_CH = 2, LFE_CH = 3, LS_CH = 4, RS_CH = 5 };

typedef struct { float re, im; } C;

typedef struct {
    int enabled;
    int in_ch;
    int out_ch;
    float m[2][CH];
} Downmix;

static int dm_init(Downmix* d, int enabled) {
    if (d == NULL) return -1;
    memset(d, 0, sizeof(*d));
    d->enabled = enabled ? 1 : 0;
    if (!enabled) { d->in_ch = CH; d->out_ch = CH; return 0; }

    d->in_ch = CH;
    d->out_ch = 2;
    float q = (float)(1.0 / sqrt(2.0));          /* 1/sqrt(2), angle 45 deg */
    float c = (float)cos(PI / 6.0);              /* cos 30 deg */
    float s = 0.5f;                              /* sin 30 deg */

    d->m[0][L_CH]  = 1.0f; d->m[1][L_CH]  = 0.0f;
    d->m[0][R_CH]  = 0.0f; d->m[1][R_CH]  = 1.0f;
    d->m[0][C_CH]  = q;    d->m[1][C_CH]  = q;
    d->m[0][LFE_CH] = q;   d->m[1][LFE_CH] = q;   /* same sign: LF-safe */
    d->m[0][LS_CH] = c;    d->m[1][LS_CH] = s;    /* equal-power pan, 30 deg */
    d->m[0][RS_CH] = s;    d->m[1][RS_CH] = c;    /* equal-power pan, 60 deg */
    return 0;
}

static int dm_process(const Downmix* d, const float* const* in,
                      float* const* out, int frames) {
    if (d == NULL || in == NULL || out == NULL || frames <= 0) return -1;
    for (int c = 0; c < d->in_ch; c++) if (in[c] == NULL) return -1;
    for (int c = 0; c < d->out_ch; c++) if (out[c] == NULL) return -1;

    if (!d->enabled) {
        if (d->in_ch != d->out_ch) return -1;
        for (int c = 0; c < d->in_ch; c++)
            memcpy(out[c], in[c], (size_t)frames * sizeof(float));
        return frames;
    }

    for (int i = 0; i < frames; i++) {
        float sl = 0.0f, sr = 0.0f;
        for (int j = 0; j < CH; j++) {
            sl += d->m[0][j] * in[j][i];
            sr += d->m[1][j] * in[j][i];
        }
        out[0][i] = sl;
        out[1][i] = sr;
    }
    return frames;
}

/* ---- independent reference: explicit per-channel arithmetic in double,
 * no matrix loop, so agreement with dm_process validates the float path. ---- */
static void analytic_downmix(const float* const* in, int n,
                             double* lt, double* rt) {
    const double q = 1.0 / sqrt(2.0);
    const double c = cos(PI / 6.0);
    const double s = 0.5;
    for (int i = 0; i < n; i++) {
        double l  = (double)in[L_CH][i];
        double r  = (double)in[R_CH][i];
        double ce = (double)in[C_CH][i];
        double lf = (double)in[LFE_CH][i];
        double ls = (double)in[LS_CH][i];
        double rs = (double)in[RS_CH][i];
        lt[i] = 1.0 * l + 0.0 * r + q * ce + q * lf + c * ls + s * rs;
        rt[i] = 0.0 * l + 1.0 * r + q * ce + q * lf + s * ls + c * rs;
    }
}

/* buggy variant for demonstration: LFE flipped out of phase */
static int dm_init_buggy_lfe(Downmix* d) {
    if (dm_init(d, 1) != 0) return -1;
    d->m[1][LFE_CH] = -d->m[1][LFE_CH];
    return 0;
}

/* ---- radix-2 FFT, same algorithm as the fft module, for band energy ---- */
static void fft(C* a, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = -2.0 * PI / len;
        C wl = { (float)cos(ang), (float)sin(ang) };
        for (int i = 0; i < n; i += len) {
            C w = { 1.0f, 0.0f };
            for (int k = 0; k < len / 2; k++) {
                C u = a[i + k];
                C v = { a[i + k + len / 2].re * w.re - a[i + k + len / 2].im * w.im,
                        a[i + k + len / 2].re * w.im + a[i + k + len / 2].im * w.re };
                a[i + k].re = u.re + v.re; a[i + k].im = u.im + v.im;
                a[i + k + len / 2].re = u.re - v.re; a[i + k + len / 2].im = u.im - v.im;
                float nr = w.re * wl.re - w.im * wl.im;
                float ni = w.re * wl.im + w.im * wl.re;
                w.re = nr; w.im = ni;
            }
        }
    }
}

/* LF band: bins 0..4 plus their conjugate mirrors. Parseval: E = sum|X|^2 / n */
static double band_energy_lf(const float* x, int n) {
    C a[N];
    for (int i = 0; i < n; i++) { a[i].re = x[i]; a[i].im = 0.0f; }
    fft(a, n);
    static const int bins[] = { 0, 1, 2, 3, 4, N - 4, N - 3, N - 2, N - 1 };
    double e = 0.0;
    for (size_t k = 0; k < sizeof(bins) / sizeof(bins[0]); k++) {
        int b = bins[k];
        e += (double)a[b].re * a[b].re + (double)a[b].im * a[b].im;
    }
    return e / n;
}

static double time_energy(const float* x, int n) {
    double e = 0.0;
    for (int i = 0; i < n; i++) e += (double)x[i] * x[i];
    return e;
}

static const char* names[CH] = { "L", "R", "C", "LFE", "Ls", "Rs" };

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    static float in[CH][N];
    static float lt[N], rt[N], mfold[N];
    static double alt[N], art[N];
    const float* inp[CH];
    float* outp[2];
    for (int c = 0; c < CH; c++) inp[c] = in[c];
    outp[0] = lt; outp[1] = rt;

    /* test tones: distinct DFT bins => orthogonal channels.
     * Only the LFE (bin 2, amplitude 0.8) lives inside the LF band. */
    const int bin[CH] = { 20, 24, 12, 2, 18, 16 };
    const float amp[CH] = { 0.5f, 0.5f, 0.5f, 0.8f, 0.5f, 0.5f };
    for (int c = 0; c < CH; c++)
        for (int i = 0; i < N; i++)
            in[c][i] = amp[c] * cosf(2.0f * PI * bin[c] * i / N);

    Downmix d;
    dm_init(&d, 1);

    /* ---------- material 1: per-channel weight table ---------- */
    printf("[1] downmix weight table (energy-conserving)\n");
    printf("ch    wL            wR            wL^2+wR^2\n");
    int norms_ok = 1, lfe_sign_ok = 0;
    for (int j = 0; j < CH; j++) {
        double wl = d.m[0][j], wr = d.m[1][j];
        double nrm = wl * wl + wr * wr;
        printf("%-4s  % .8f  % .8f  %.12f\n", names[j], wl, wr, nrm);
        if (fabs(nrm - 1.0) > 1e-6) norms_ok = 0;
    }
    lfe_sign_ok = d.m[0][LFE_CH] > 0.0f && d.m[1][LFE_CH] > 0.0f;
    check("column_norm_unit", norms_ok);
    check("lfe_in_phase", lfe_sign_ok);

    /* ---------- run matrix path ---------- */
    check("process_ok", dm_process(&d, inp, outp, N) == N);
    check("bad_args", dm_process(&d, NULL, outp, N) == -1 &&
                      dm_process(&d, inp, outp, 0) == -1);

    /* ---------- material 2: low-frequency energy comparison ---------- */
    double e_lfe_in = band_energy_lf(in[LFE_CH], N);
    double e_lf_stereo = band_energy_lf(lt, N) + band_energy_lf(rt, N);
    for (int i = 0; i < N; i++)
        mfold[i] = (lt[i] + rt[i]) * (float)(1.0 / sqrt(2.0));
    double e_lf_mono = band_energy_lf(mfold, N);

    /* buggy opposite-polarity LFE: same matrix path, one sign flipped */
    static float blt[N], brt[N], bm[N];
    float* boutp[2] = { blt, brt };
    Downmix bad;
    dm_init_buggy_lfe(&bad);
    dm_process(&bad, inp, boutp, N);
    double e_bad_stereo = band_energy_lf(blt, N) + band_energy_lf(brt, N);
    for (int i = 0; i < N; i++)
        bm[i] = (blt[i] + brt[i]) * (float)(1.0 / sqrt(2.0));
    double e_bad_mono = band_energy_lf(bm, N);

    double e_lfe_expect = (double)N / 2.0 * 0.8 * 0.8;   /* cosine, amp 0.8 */
    printf("[2] LF band energy (bins 0..4, Parseval)\n");
    printf("input LFE LF energy              = %.6f  (theory %.6f)\n", e_lfe_in, e_lfe_expect);
    printf("correct matrix: stereo LF energy = %.6f\n", e_lf_stereo);
    printf("correct matrix: mono-fold LF     = %.6f\n", e_lf_mono);
    printf("buggy matrix:   stereo LF energy = %.6f  (sign flip hides here)\n", e_bad_stereo);
    printf("buggy matrix:   mono-fold LF     = %.6f  (canceled)\n", e_bad_mono);
    check("lf_stereo_conserved", fabs(e_lf_stereo - e_lfe_in) / e_lfe_in < 1e-3);
    check("lf_mono_preserved", fabs(e_lf_mono - e_lfe_in) / e_lfe_in < 1e-3);
    check("buggy_mono_cancels", e_bad_mono < 1e-3 * e_lfe_in);
    check("buggy_stereo_hides_loss", fabs(e_bad_stereo - e_lfe_in) / e_lfe_in < 1e-3);

    /* total energy identity: orthogonal inputs => out energy == sum in energy */
    double e_in = 0.0;
    for (int c = 0; c < CH; c++) e_in += time_energy(in[c], N);
    double e_out = time_energy(lt, N) + time_energy(rt, N);
    printf("total energy: in=%.6f out=%.6f relerr=%.3e\n",
           e_in, e_out, fabs(e_out - e_in) / e_in);
    check("total_energy_identity", fabs(e_out - e_in) / e_in < 1e-4);

    /* ---------- material 3: difference vs per-channel analytic compute ---------- */
    analytic_downmix(inp, N, alt, art);
    double maxd = 0.0;
    for (int i = 0; i < N; i++) {
        double dl = fabs((double)lt[i] - alt[i]);
        double dr = fabs((double)rt[i] - art[i]);
        if (dl > maxd) maxd = dl;
        if (dr > maxd) maxd = dr;
    }
    double e_analytic = 0.0;
    for (int i = 0; i < N; i++) e_analytic += alt[i] * alt[i] + art[i] * art[i];
    printf("[3] analytic (double, per-channel explicit) comparison\n");
    printf("max |matrix_float - analytic_double| = %.3e\n", maxd);
    printf("analytic output energy = %.6f, matrix output energy = %.6f\n",
           e_analytic, e_out);
    check("analytic_diff_ulp", maxd < 1e-5);
    check("analytic_energy_match", fabs(e_analytic - e_out) / e_analytic < 1e-5);

    /* dump raw frames so an external (Python) reference can recompute
     * the exact same samples independently. Channel-planar float32. */
    FILE* fi = fopen("downmix_in.f32", "wb");
    FILE* fl = fopen("downmix_lt.f32", "wb");
    FILE* fr = fopen("downmix_rt.f32", "wb");
    if (fi) { for (int c = 0; c < CH; c++) fwrite(in[c], sizeof(float), N, fi); fclose(fi); }
    if (fl) { fwrite(lt, sizeof(float), N, fl); fclose(fl); }
    if (fr) { fwrite(rt, sizeof(float), N, fr); fclose(fr); }

    /* analytic column-norm proof with exact double constants,
     * independent of the float values stored in d.m */
    const double wl_a[CH] = { 1.0, 0.0, 1.0 / sqrt(2.0), 1.0 / sqrt(2.0),
                              cos(PI / 6.0), 0.5 };
    const double wr_a[CH] = { 0.0, 1.0, 1.0 / sqrt(2.0), 1.0 / sqrt(2.0),
                              0.5, cos(PI / 6.0) };
    int analytic_norms = 1;
    for (int j = 0; j < CH; j++)
        if (fabs(wl_a[j] * wl_a[j] + wr_a[j] * wr_a[j] - 1.0) > 1e-12)
            analytic_norms = 0;
    check("analytic_column_norms", analytic_norms);

    /* ---------- material 4: bypass is byte-identical ---------- */
    static float bin_buf[CH][N], bout_buf[CH][N];
    srand(1234);
    for (int c = 0; c < CH; c++)
        for (int i = 0; i < N; i++)
            bin_buf[c][i] = rand() / (float)RAND_MAX * 2.0f - 1.0f;
    const float* binp[CH];
    float* boutp2[CH];
    for (int c = 0; c < CH; c++) { binp[c] = bin_buf[c]; boutp2[c] = bout_buf[c]; }

    Downmix off;
    dm_init(&off, 0);
    check("bypass_counts", off.in_ch == CH && off.out_ch == CH);
    check("bypass_ok", dm_process(&off, binp, boutp2, N) == N);
    int byte_same = memcmp(bin_buf, bout_buf, sizeof(bin_buf)) == 0;
    check("bypass_byte_identical", byte_same);

    FILE* f1 = fopen("bypass_in.f32", "wb");
    FILE* f2 = fopen("bypass_out.f32", "wb");
    if (f1 && f2) {
        for (int c = 0; c < CH; c++) {
            fwrite(bin_buf[c], sizeof(float), N, f1);
            fwrite(bout_buf[c], sizeof(float), N, f2);
        }
    }
    if (f1) fclose(f1);
    if (f2) fclose(f2);
    printf("[4] bypass raw dumps: bypass_in.f32 / bypass_out.f32 "
           "(Makefile verifies with cmp)\n");

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
