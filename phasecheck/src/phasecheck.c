/* phasecheck: detect phase cancellation when the same material is stacked twice.
 *
 * ------------------------------------------------------------------ theory
 * Stack two takes A and B.  Per band, over a window of W frames and the
 * band's B frequency bins, form the band-limited cross-spectrum and powers:
 *
 *   C   = sum_f sum_k Xa[k,f] conj(Xb[k,f])      (Re/Im retained)
 *   P_a = sum_f sum_k |Xa[k,f]|^2
 *   P_b = sum_f sum_k |Xb[k,f]|^2
 *   rho = Re(C) / sqrt(P_a P_b)                   signed coherence
 *   gam = |C|   / sqrt(P_a P_b)                   magnitude coherence
 *
 * The summed (coherent) signal power is
 *   P_sum = P_a + P_b + 2 Re(C)
 * so, normalising by an incoherent power sum (two uncorrelated sources),
 *   G = P_sum / (P_a + P_b) = 1 + rho
 *   dL = 10 log10(1 + rho)               dB relative to incoherent summing.
 *
 *   rho = +1  identical takes          -> +3.01 dB
 *   rho =  0  uncorrelated takes       ->  0 dB
 *   rho = -1  inverted copy            -> -inf dB  (printed with -60 dB floor)
 *
 * So the signed coherence *is* the level shortfall: "stacked twice but
 * lighter" is exactly rho < 0.  The relation assumes comparable take power
 * inside the band; the report prints both take levels per band so the
 * assumption can be checked by eye.
 *
 * --------------------------------------------------------------- threshold
 * Under H0 (the takes are independent in this band) the band-window statistic
 * has nu = B*W independent complex samples.  The magnitude-squared coherence
 * of two independent complex Gaussian vectors of length nu is Beta(1, nu-1),
 * hence
 *   gam^2 ~ Beta(1, nu-1),   and since |rho| <= gam,
 *   P(rho <= -t) <= P(gam >= t) = (1 - t^2)^(nu-1).
 * Choose alpha_band = alpha / n_active (Bonferroni across bands) and solve:
 *   t_band = sqrt(1 - alpha_band^(1/(nu-1))).                 (nu >= 2)
 * A band-window is flagged when rho_band < -t_band.  `alpha` is the
 * per-window, across-bands family-wise false-alarm probability; the report
 * prints the expected number of false-alarm windows, nwin * alpha.
 *
 * The bound is exact given independent complex-Gaussian samples.  That holds
 * well for noise-like material.  It does *not* hold for stationary tonal
 * content, where a pure tone keeps its phase across frames, nu_eff collapses
 * to the number of independent bins, and the true false-alarm rate exceeds
 * nominal.  --calib measures both regimes and prints the ratio; see README.
 *
 * --------------------------------------------------------------- behaviour
 * Detection is never silent.  Cancellation is always printed and echoed to
 * stderr.  --out refuses to write the mix unless --accept-cancellation is
 * given.  Exit: 0 clean, 3 cancellation detected, 4 mix refused.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_BANDS 8
#define PC_PI    3.14159265358979323846
#define PC_SQRT2 1.41421356237309504880
#define PC_EPS   1e-30
#define GDB_FLOOR (-60.0)

typedef struct { double re, im; } D;

static const double BAND_CENTER[MAX_BANDS] =
    { 125, 250, 500, 1000, 2000, 4000, 8000, 16000 };

typedef struct { int fs, seg_len, win, hop, max_lag; double alpha; } Cfg;

typedef struct {
    int active, k0, k1, bins;
    int nu;                     /* B*W, independent samples in a band-window */
    double thr;                 /* detection threshold on rho, see header */
    long windows, flagged;
    double rho_total, g_total_db;
    double p_a, p_b;            /* band powers over all frames */
} BandMeta;

typedef struct {                /* one analysis frame (one FFT segment) */
    double pa[MAX_BANDS], pb[MAX_BANDS], cr[MAX_BANDS], ci[MAX_BANDS];
    double pa_t, pb_t, cr_t, ci_t;
} Frame;

typedef struct {                /* one detection window (W frames) */
    int f0, f1;
    double t0, t1;
    double rho[MAX_BANDS], g_db[MAX_BANDS], gamma2[MAX_BANDS];
    int cancel[MAX_BANDS], nband_cancel;
    double rho_bb, g_bb_db;
} Win;

typedef struct {
    int fs, L, nfft, nseg, W, hop, nwin, n_active;
    double alpha, alpha_band;
    BandMeta band[MAX_BANDS];
    Frame* frame;
    Win* win;
    long flagged_windows_total;
    double expected_false_alarms;
    double rms_a, rms_b, rms_sum, peak_a, peak_b, peak_sum;
    int best_lag, max_lag;
    double best_lag_sec, peak_corr, corr_zero;
    double* xcorr;              /* normalized corr for m = -max_lag..+max_lag */
} Result;

/* ------------------------------------------------------------------ utils */

static int prev_pow2(int n) {
    int p = 1;
    while ((p << 1) <= n && p < (1 << 24)) p <<= 1;
    return p;
}
static int next_pow2(int n) {
    int p = 1;
    while (p < n && p < (1 << 24)) p <<= 1;
    return p;
}
static int clampi(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }
static double db10(double x) {
    if (!(x > 0.0)) return -200.0;
    return 10.0 * log10(x);
}
static double gain_db(double rho) {
    double g = 10.0 * log10(1.0 + rho);
    return g < GDB_FLOOR ? GDB_FLOOR : g;
}
static double now_ms(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec * 1000.0 + (double)t.tv_nsec / 1e6;
}

typedef struct { uint64_t s; } Rng;
static void rng_seed(Rng* r, uint64_t seed) { r->s = seed ? seed : 1; }
static double rng_uniform(Rng* r) {
    uint64_t x = r->s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    r->s = x;
    return (double)((x * 2685821657736338717ULL) >> 11) / 9007199254740992.0;
}
static double rng_gauss(Rng* r) {
    double u1 = rng_uniform(r), u2 = rng_uniform(r);
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PC_PI * u2);
}

/* -------------------------------------------------------------------- fft */

static void fft_radix2(D* a, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { D t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * PC_PI / (double)len * (inverse ? 1.0 : -1.0);
        D wl = { cos(ang), sin(ang) };
        for (int i = 0; i < n; i += len) {
            D w = { 1.0, 0.0 };
            for (int k = 0; k < len / 2; k++) {
                D u = a[i + k];
                D v = { a[i + k + len / 2].re * w.re - a[i + k + len / 2].im * w.im,
                        a[i + k + len / 2].re * w.im + a[i + k + len / 2].im * w.re };
                a[i + k].re = u.re + v.re; a[i + k].im = u.im + v.im;
                a[i + k + len / 2].re = u.re - v.re;
                a[i + k + len / 2].im = u.im - v.im;
                double nr = w.re * wl.re - w.im * wl.im;
                double ni = w.re * wl.im + w.im * wl.re;
                w.re = nr; w.im = ni;
            }
        }
    }
    if (inverse) for (int i = 0; i < n; i++) { a[i].re /= n; a[i].im /= n; }
}

/* -------------------------------------------------------------------- wav */

typedef struct { double* x; int n; int fs; int channels; char err[160]; } Wav;

static uint32_t rd32(const unsigned char* p) {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}
static uint16_t rd16(const unsigned char* p) {
    return (uint16_t)((uint16_t)p[0] | ((uint16_t)p[1] << 8));
}

/* Loads one channel of a 16-bit PCM RIFF/WAVE file as doubles in [-1,1).
 * Channels are selected, never averaged: averaging channels could itself
 * cancel and would corrupt the measurement. */
static int load_wav(const char* path, int ch, Wav* w) {
    memset(w, 0, sizeof(*w));
    FILE* f = fopen(path, "rb");
    if (!f) { snprintf(w->err, sizeof(w->err), "cannot open %s", path); return -1; }
    unsigned char hdr[12];
    if (fread(hdr, 1, 12, f) != 12 || memcmp(hdr, "RIFF", 4) || memcmp(hdr + 8, "WAVE", 4)) {
        snprintf(w->err, sizeof(w->err), "%s: not a RIFF/WAVE file", path);
        fclose(f); return -1;
    }
    int fmt_seen = 0, have_data = 0;
    long data_off = 0, data_len = 0;
    unsigned char ck[8];
    while (fread(ck, 1, 8, f) == 8) {
        uint32_t sz = rd32(ck + 4);
        if (!memcmp(ck, "fmt ", 4)) {
            unsigned char fb[40];
            uint32_t want = sz < sizeof(fb) ? sz : (uint32_t)sizeof(fb);
            if (want < 16 || fread(fb, 1, want, f) != want) break;
            int tag = rd16(fb), channels = rd16(fb + 2);
            uint32_t fs = rd32(fb + 4), br = rd32(fb + 8);
            int align = rd16(fb + 12), bits = rd16(fb + 14);
            if (tag == 0xFFFE) {
                snprintf(w->err, sizeof(w->err), "%s: WAVE_FORMAT_EXTENSIBLE unsupported", path);
                fclose(f); return -1;
            }
            if (tag != 1 || bits != 16 || channels < 1 ||
                align != channels * 2 || br != fs * (uint32_t)align) {
                snprintf(w->err, sizeof(w->err),
                         "%s: need 16-bit PCM, got tag=%d bits=%d ch=%d", path, tag, bits, channels);
                fclose(f); return -1;
            }
            if (sz > want) fseek(f, (long)(sz - want), SEEK_CUR);
            w->fs = (int)fs; w->channels = channels; fmt_seen = 1;
        } else if (!memcmp(ck, "data", 4)) {
            data_off = ftell(f); data_len = (long)sz;
            if (sz == 0 || sz == 0xFFFFFFFFu) {
                fseek(f, 0, SEEK_END); data_len = ftell(f) - data_off;
            }
            have_data = 1;
            break;
        } else {
            fseek(f, (long)(sz + (sz & 1)), SEEK_CUR);
        }
    }
    if (!fmt_seen || !have_data || data_len <= 0) {
        snprintf(w->err, sizeof(w->err), "%s: missing fmt/data chunk", path);
        fclose(f); return -1;
    }
    if (ch < 0) ch = 0;
    if (ch >= w->channels) {
        snprintf(w->err, sizeof(w->err), "%s: channel %d out of range (have %d)",
                 path, ch, w->channels);
        fclose(f); return -1;
    }
    long frames = data_len / (2L * w->channels);
    double* x = (double*)malloc(sizeof(double) * (size_t)(frames > 0 ? frames : 1));
    unsigned char* buf = (unsigned char*)malloc(65536);
    if (!x || !buf) { free(x); free(buf); fclose(f);
        snprintf(w->err, sizeof(w->err), "out of memory"); return -1; }
    fseek(f, data_off, SEEK_SET);
    long fr = 0;
    size_t frame_bytes = 2u * (size_t)w->channels;
    while (fr < frames) {
        size_t chunk = (size_t)(frames - fr);
        size_t maxf = 65536 / frame_bytes;
        if (chunk > maxf) chunk = maxf;
        size_t got = fread(buf, frame_bytes, chunk, f);
        if (got == 0) break;
        for (size_t i = 0; i < got; i++) {
            int16_t v = (int16_t)rd16(buf + i * frame_bytes + 2u * (size_t)ch);
            x[fr + (long)i] = (double)v / 32768.0;
        }
        fr += (long)got;
    }
    free(buf);
    fclose(f);
    if (fr <= 0) { free(x); snprintf(w->err, sizeof(w->err), "%s: no samples", path); return -1; }
    w->x = x; w->n = (int)fr;
    return 0;
}

static int write_wav16(const char* path, const float* x, int n, int fs, int* clipped) {
    FILE* f = fopen(path, "wb");
    if (!f) return -1;
    uint32_t data_bytes = (uint32_t)n * 2u;
    uint32_t riff = 36u + data_bytes;
    unsigned char h[44];
    memset(h, 0, sizeof(h));
    memcpy(h, "RIFF", 4);
    h[4] = (unsigned char)riff; h[5] = (unsigned char)(riff >> 8);
    h[6] = (unsigned char)(riff >> 16); h[7] = (unsigned char)(riff >> 24);
    memcpy(h + 8, "WAVEfmt ", 8);
    h[16] = 16; h[20] = 1; h[22] = 1;
    h[24] = (unsigned char)fs; h[25] = (unsigned char)(fs >> 8);
    h[26] = (unsigned char)(fs >> 16); h[27] = (unsigned char)(fs >> 24);
    uint32_t br = (uint32_t)fs * 2u;
    h[28] = (unsigned char)br; h[29] = (unsigned char)(br >> 8);
    h[30] = (unsigned char)(br >> 16); h[31] = (unsigned char)(br >> 24);
    h[32] = 2; h[34] = 16;
    memcpy(h + 36, "data", 4);
    h[40] = (unsigned char)data_bytes; h[41] = (unsigned char)(data_bytes >> 8);
    h[42] = (unsigned char)(data_bytes >> 16); h[43] = (unsigned char)(data_bytes >> 24);
    if (fwrite(h, 1, 44, f) != 44) { fclose(f); return -1; }
    int clip = 0;
    for (int i = 0; i < n; i++) {
        double v = x[i];
        if (v > 0.999969) { v = 0.999969; clip++; }
        if (v < -1.0) { v = -1.0; clip++; }
        int16_t s = (int16_t)lrint(v * 32767.0);
        unsigned char b[2] = { (unsigned char)(s & 0xFF), (unsigned char)((s >> 8) & 0xFF) };
        if (fwrite(b, 1, 2, f) != 2) { fclose(f); return -1; }
    }
    fclose(f);
    if (clipped) *clipped = clip;
    return 0;
}

/* ------------------------------------------------------------------ bands */

static int map_bands(int nfft, int fs, BandMeta* band) {
    int n = 0;
    for (int bi = 0; bi < MAX_BANDS; bi++) {
        memset(&band[bi], 0, sizeof(band[bi]));
        band[bi].g_total_db = GDB_FLOOR;
        if (BAND_CENTER[bi] >= (double)fs * 0.5) continue;   /* above Nyquist */
        int k0 = (int)lround(BAND_CENTER[bi] / PC_SQRT2 * (double)nfft / (double)fs);
        int k1 = (int)lround(BAND_CENTER[bi] * PC_SQRT2 * (double)nfft / (double)fs);
        k0 = clampi(k0, 1, nfft / 2 - 1);
        k1 = clampi(k1, 1, nfft / 2 - 1);
        if (k1 < k0) { int t = k0; k0 = k1; k1 = t; }
        if (k1 < 1 || k0 > nfft / 2 - 1) continue;
        band[bi].k0 = k0; band[bi].k1 = k1;
        band[bi].bins = k1 - k0 + 1;
        band[bi].active = 1;
        n++;
    }
    return n;
}

/* --------------------------------------------------------------- analysis */

static double corr_zero_(const double* a, const double* b, int n) {
    double sa = 0, sb = 0, sab = 0;
    for (int i = 0; i < n; i++) { sa += a[i]*a[i]; sb += b[i]*b[i]; sab += a[i]*b[i]; }
    double den = sqrt(sa * sb);
    return den > PC_EPS ? sab / den : 0.0;
}

static void free_result(Result* R) {
    free(R->frame); free(R->win); free(R->xcorr);
    R->frame = NULL; R->win = NULL; R->xcorr = NULL;
}

static int analyze(const double* xa, int na, const double* xb, int nb,
                   const Cfg* cfg, Result* R) {
    memset(R, 0, sizeof(*R));
    int L = na < nb ? na : nb;
    if (L < 2048) return -1;

    int nfft = cfg->seg_len > 0 ? prev_pow2(cfg->seg_len) : prev_pow2(L / 64);
    if (nfft < 512) nfft = 512;
    if (nfft > 16384) nfft = 16384;
    int lim = prev_pow2(L);
    if (nfft > lim) nfft = lim;
    int nseg = L / nfft;
    if (nseg < 4) return -1;

    R->fs = cfg->fs; R->L = L; R->nfft = nfft; R->nseg = nseg;
    R->alpha = (cfg->alpha > 0.0 && cfg->alpha < 1.0) ? cfg->alpha : 0.01;
    R->n_active = map_bands(nfft, cfg->fs, R->band);
    if (R->n_active < 1) return -1;

    int W = cfg->win > 0 ? cfg->win : 32;
    if (W > nseg) W = nseg;
    if (W < 2) W = 2;
    int hop = cfg->hop > 0 ? cfg->hop : (W / 2 > 0 ? W / 2 : 1);
    if (hop > W) hop = W;
    if (hop < 1) hop = 1;
    int nwin = (nseg - W) / hop + 1;
    if (nwin < 1) return -1;
    R->W = W; R->hop = hop; R->nwin = nwin;
    R->alpha_band = R->alpha / (double)R->n_active;
    R->expected_false_alarms = (double)nwin * R->alpha;

    /* per-band threshold from nu = bins * W independent complex samples */
    for (int bi = 0; bi < MAX_BANDS; bi++) {
        BandMeta* b = &R->band[bi];
        if (!b->active) continue;
        int nu = b->bins * W;
        b->nu = nu;
        if (nu < 2) { b->thr = 1.0; continue; }
        double q = pow(R->alpha_band, 1.0 / (double)(nu - 1));
        double t2 = 1.0 - q;
        if (t2 < 0.0) t2 = 0.0;
        if (t2 > 1.0) t2 = 1.0;
        b->thr = sqrt(t2);
    }

    /* ---------------- per-frame spectra ---------------- */
    R->frame = (Frame*)calloc((size_t)nseg, sizeof(Frame));
    D* fa = (D*)calloc((size_t)nfft, sizeof(D));
    D* fb = (D*)calloc((size_t)nfft, sizeof(D));
    if (!R->frame || !fa || !fb) { free_result(R); free(fa); free(fb); return -1; }
    for (int si = 0; si < nseg; si++) {
        Frame* F = &R->frame[si];
        long off = (long)si * nfft;
        for (int i = 0; i < nfft; i++) {
            fa[i].re = xa[off + i]; fa[i].im = 0;
            fb[i].re = xb[off + i]; fb[i].im = 0;
        }
        fft_radix2(fa, nfft, 0);
        fft_radix2(fb, nfft, 0);
        for (int k = 1; k < nfft / 2; k++) {
            F->pa_t += fa[k].re * fa[k].re + fa[k].im * fa[k].im;
            F->pb_t += fb[k].re * fb[k].re + fb[k].im * fb[k].im;
            F->cr_t += fa[k].re * fb[k].re + fa[k].im * fb[k].im;
            F->ci_t += fa[k].re * fb[k].im - fa[k].im * fb[k].re;
        }
        for (int bi = 0; bi < MAX_BANDS; bi++) {
            if (!R->band[bi].active) continue;
            double pa = 0, pb = 0, cr = 0, ci = 0;
            for (int k = R->band[bi].k0; k <= R->band[bi].k1; k++) {
                pa += fa[k].re * fa[k].re + fa[k].im * fa[k].im;
                pb += fb[k].re * fb[k].re + fb[k].im * fb[k].im;
                cr += fa[k].re * fb[k].re + fa[k].im * fb[k].im;
                ci += fa[k].re * fb[k].im - fa[k].im * fb[k].re;
            }
            F->pa[bi] = pa; F->pb[bi] = pb; F->cr[bi] = cr; F->ci[bi] = ci;
        }
    }
    free(fa); free(fb);

    /* per-band totals over all frames (the coherence table) */
    for (int bi = 0; bi < MAX_BANDS; bi++) {
        BandMeta* b = &R->band[bi];
        if (!b->active) continue;
        double pa = 0, pb = 0, cr = 0;
        for (int si = 0; si < nseg; si++) {
            pa += R->frame[si].pa[bi]; pb += R->frame[si].pb[bi]; cr += R->frame[si].cr[bi];
        }
        b->p_a = pa; b->p_b = pb;
        double den = sqrt(pa * pb);
        b->rho_total = den > PC_EPS ? cr / den : 0.0;
        b->g_total_db = den > PC_EPS ? gain_db(b->rho_total) : GDB_FLOOR;
    }

    /* ---------------- detection windows ---------------- */
    R->win = (Win*)calloc((size_t)nwin, sizeof(Win));
    if (!R->win) { free_result(R); return -1; }
    for (int wi = 0; wi < nwin; wi++) {
        Win* V = &R->win[wi];
        V->f0 = wi * hop;
        V->f1 = V->f0 + W;
        V->t0 = (double)((long)V->f0 * nfft) / (double)cfg->fs;
        V->t1 = (double)((long)V->f1 * nfft) / (double)cfg->fs;
        double pa_t = 0, pb_t = 0, cr_t = 0, ci_t = 0;
        for (int si = V->f0; si < V->f1; si++) {
            const Frame* F = &R->frame[si];
            pa_t += F->pa_t; pb_t += F->pb_t; cr_t += F->cr_t; ci_t += F->ci_t;
        }
        double den_t = sqrt(pa_t * pb_t);
        V->rho_bb = den_t > PC_EPS ? cr_t / den_t : 0.0;
        V->g_bb_db = den_t > PC_EPS ? gain_db(V->rho_bb) : GDB_FLOOR;
        for (int bi = 0; bi < MAX_BANDS; bi++) {
            BandMeta* b = &R->band[bi];
            if (!b->active) continue;
            double pa = 0, pb = 0, cr = 0, ci = 0;
            for (int si = V->f0; si < V->f1; si++) {
                pa += R->frame[si].pa[bi]; pb += R->frame[si].pb[bi];
                cr += R->frame[si].cr[bi]; ci += R->frame[si].ci[bi];
            }
            double den = sqrt(pa * pb);
            V->cancel[bi] = 0;
            if (den <= PC_EPS) { V->rho[bi] = 0.0; V->g_db[bi] = GDB_FLOOR; V->gamma2[bi] = 0.0; continue; }
            V->gamma2[bi] = (cr * cr + ci * ci) / (pa * pb);
            V->rho[bi] = cr / den;
            V->g_db[bi] = gain_db(V->rho[bi]);
            b->windows++;
            if (V->rho[bi] < -b->thr) {
                V->cancel[bi] = 1;
                V->nband_cancel++;
                b->flagged++;
                R->flagged_windows_total++;
            }
        }
    }

    /* ---------------- levels over the analysed region ---------------- */
    int used = nseg * nfft;
    double sa = 0, sb = 0, ss = 0, pka = 0, pkb = 0, pks = 0;
    for (int i = 0; i < used; i++) {
        double va = xa[i], vb = xb[i], vs = va + vb;
        sa += va * va; sb += vb * vb; ss += vs * vs;
        if (fabs(va) > pka) pka = fabs(va);
        if (fabs(vb) > pkb) pkb = fabs(vb);
        if (fabs(vs) > pks) pks = fabs(vs);
    }
    R->rms_a = sqrt(sa / used);
    R->rms_b = sqrt(sb / used);
    R->rms_sum = sqrt(ss / used);
    R->peak_a = pka; R->peak_b = pkb; R->peak_sum = pks;
    R->corr_zero = corr_zero_(xa, xb, used);
    R->max_lag = cfg->max_lag;

    /* ---------------- normalized cross-correlation vs lag ---------------- */
    R->best_lag = 0; R->best_lag_sec = 0.0; R->peak_corr = 0.0;
    if (cfg->max_lag > 0) {
        int nwin_lag = used < 32768 ? used : 32768;
        int N = next_pow2(2 * nwin_lag);
        D* fa2 = (D*)calloc((size_t)N, sizeof(D));
        D* fb2 = (D*)calloc((size_t)N, sizeof(D));
        double* ea = (double*)malloc(sizeof(double) * (size_t)(nwin_lag + 1));
        double* eb = (double*)malloc(sizeof(double) * (size_t)(nwin_lag + 1));
        if (fa2 && fb2 && ea && eb) {
            for (int i = 0; i < nwin_lag; i++) {
                fa2[i].re = xa[i]; fb2[i].re = xb[i];
                ea[i + 1] = ea[i] + xa[i] * xa[i];
                eb[i + 1] = eb[i] + xb[i] * xb[i];
            }
            fft_radix2(fa2, N, 0);
            fft_radix2(fb2, N, 0);
            D* xc = (D*)calloc((size_t)N, sizeof(D));
            if (xc) {
                for (int k = 0; k < N; k++)
                    xc[k] = (D){ fa2[k].re * fb2[k].re + fa2[k].im * fb2[k].im,
                                 fa2[k].re * fb2[k].im - fa2[k].im * fb2[k].re };
                fft_radix2(xc, N, 1);
                int ml = cfg->max_lag;
                if (ml > nwin_lag / 4) ml = nwin_lag / 4;
                R->max_lag = ml;
                R->xcorr = (double*)calloc((size_t)(2 * ml + 1), sizeof(double));
                double best_abs = -2.0;
                for (int m = -ml; m <= ml; m++) {
                    int mi = m & (N - 1);
                    int lo = m > 0 ? 0 : -m;
                    int hi = m > 0 ? nwin_lag - m : nwin_lag;
                    double r = 0.0;
                    if (hi > lo) {
                        double den = sqrt((ea[hi] - ea[lo]) * (eb[hi + m] - eb[lo + m]));
                        r = den > PC_EPS ? xc[mi].re / den : 0.0;
                    }
                    if (R->xcorr) R->xcorr[m + ml] = r;
                    if (fabs(r) > best_abs) { best_abs = fabs(r); R->best_lag = m; }
                }
                R->best_lag_sec = (double)R->best_lag / (double)cfg->fs;
                R->peak_corr = best_abs;
                free(xc);
            }
        }
        free(fa2); free(fb2); free(ea); free(eb);
    }
    return 0;
}

/* ---------------------------------------------------------------- report */

typedef struct {
    int any;                    /* any band flagged in this window */
    double t0, t1;
} AnySpan;

static void print_report(const Result* R, const char* name_a, const char* name_b,
                         double ms_analysis) {
    printf("phasecheck: A=%s B=%s\n", name_a, name_b);
    printf("fs=%d\n", R->fs);
    printf("samples_used=%d\n", R->nseg * R->nfft);
    printf("tail_samples_dropped=%d\n", R->L - R->nseg * R->nfft);
    printf("seg_len=%d\n", R->nfft);
    printf("frames=%d\n", R->nseg);
    printf("window_frames=%d\n", R->W);
    printf("window_hop=%d\n", R->hop);
    printf("windows=%d\n", R->nwin);
    printf("window_sec=%.3f\n", (double)R->W * (double)R->nfft / (double)R->fs);
    printf("bands_active=%d\n", R->n_active);
    printf("alpha=%.4f\n", R->alpha);
    printf("alpha_band=%.6f\n", R->alpha_band);
    printf("expected_false_alarm_windows=%.2f\n", R->expected_false_alarms);

    printf("\n[levels-before-and-after]\n");
    printf("rms_a=%.6f\n", R->rms_a);
    printf("rms_b=%.6f\n", R->rms_b);
    printf("rms_sum=%.6f\n", R->rms_sum);
    printf("peak_a=%.6f\n", R->peak_a);
    printf("peak_b=%.6f\n", R->peak_b);
    printf("peak_sum=%.6f\n", R->peak_sum);
    /* Actual level change vs the single take, and the change that would have
     * been expected had the takes been uncorrelated (power sum).  The
     * shortfall is the gap between them: negative means the stack is lighter
     * than an incoherent sum, i.e. cancellation.  Perfect inversion makes the
     * actual delta -200 dB and the shortfall correspondingly floored. */
    double pa = R->rms_a * R->rms_a, pb = R->rms_b * R->rms_b;
    double sum_p = R->rms_sum * R->rms_sum;
    double d_single = db10(sum_p / (pa > PC_EPS ? pa : 1.0));
    double d_pred = db10((pa + pb) / (pa > PC_EPS ? pa : 1.0));
    double shortfall = d_single - d_pred;
    if (shortfall < GDB_FLOOR) shortfall = GDB_FLOOR;
    printf("level_delta_vs_single_take_db=%.3f\n", d_single);
    printf("level_delta_if_incoherent_db=%.3f\n", d_pred);
    printf("shortfall_vs_incoherent_db=%.3f\n", shortfall);
    printf("sum_is_silent=%d\n", sum_p <= PC_EPS ? 1 : 0);

    printf("\n[cross-correlation]\n");
    printf("corr_zero=%.6f\n", R->corr_zero);
    printf("peak_corr=%.6f\n", R->peak_corr);
    printf("best_lag_samples=%d\n", R->best_lag);
    printf("best_lag_ms=%.3f\n", R->best_lag_sec * 1000.0);
    printf("max_lag_searched=%d\n", R->max_lag);
    if (R->xcorr && R->max_lag > 0) {
        printf("\n[cross-correlation-lags]\n  lag   ms      corr\n");
        int step = R->max_lag / 10;
        if (step < 1) step = 1;
        for (int m = -R->max_lag; m <= R->max_lag; m += step)
            printf("%5d  %6.2f  %7.4f\n", m, (double)m * 1000.0 / R->fs,
                   R->xcorr[m + R->max_lag]);
        printf("%5d  %6.2f  %7.4f  (peak)\n", R->best_lag,
               R->best_lag_sec * 1000.0, R->xcorr[R->best_lag + R->max_lag]);
    }

    printf("\n[coherence-by-band]\n");
    printf("band_hz  bins  nu     thr     windows  flagged  flag_pct  rho_total  gain_total_db\n");
    for (int bi = 0; bi < MAX_BANDS; bi++) {
        const BandMeta* b = &R->band[bi];
        if (!b->active) continue;
        double pct = b->windows ? 100.0 * (double)b->flagged / (double)b->windows : 0.0;
        printf("%7.0f  %4d  %5d  %6.4f  %7ld  %7ld  %7.1f  %9.4f  %13.2f\n",
               BAND_CENTER[bi], b->bins, b->nu, b->thr,
               b->windows, b->flagged, pct, b->rho_total, b->g_total_db);
    }

    /* per-band cancellation intervals */
    printf("\n[cancel-intervals]\n");
    int nint = 0;
    for (int bi = 0; bi < MAX_BANDS; bi++) {
        const BandMeta* b = &R->band[bi];
        if (!b->active) continue;
        int wi = 0;
        while (wi < R->nwin) {
            if (!R->win[wi].cancel[bi]) { wi++; continue; }
            int start = wi;
            double worst = 0.0, sum_rho = 0.0;
            int cnt = 0;
            while (wi < R->nwin && R->win[wi].cancel[bi]) {
                if (R->win[wi].g_db[bi] < worst) worst = R->win[wi].g_db[bi];
                sum_rho += R->win[wi].rho[bi];
                cnt++;
                wi++;
            }
            printf("interval_%d  band=%.0fHz  t=%.3f..%.3f s  windows=%d..%d  "
                   "mean_rho=%.3f  worst_gain_db=%.2f\n",
                   ++nint, BAND_CENTER[bi], R->win[start].t0, R->win[wi - 1].t1,
                   start, wi - 1, sum_rho / cnt, worst);
        }
    }
    if (nint == 0) printf("(none)\n");
    printf("cancel_intervals=%d\n", nint);

    /* broadband intervals: any band flagged, merged */
    printf("\n[cancel-intervals-any-band]\n");
    int nany = 0;
    int wi = 0;
    while (wi < R->nwin) {
        if (R->win[wi].nband_cancel == 0) { wi++; continue; }
        int start = wi;
        int union_bands = 0;
        while (wi < R->nwin && R->win[wi].nband_cancel > 0) {
            for (int bi = 0; bi < MAX_BANDS; bi++)
                if (R->win[wi].cancel[bi]) union_bands |= 1 << bi;
            wi++;
        }
        printf("span_%d  t=%.3f..%.3f s  windows=%d..%d  bands=",
               ++nany, R->win[start].t0, R->win[wi - 1].t1, start, wi - 1);
        int first = 1;
        for (int bi = 0; bi < MAX_BANDS; bi++) {
            if (!(union_bands & (1 << bi))) continue;
            printf("%s%.0f", first ? "" : ",", BAND_CENTER[bi]);
            first = 0;
        }
        printf("\n");
    }
    if (nany == 0) printf("(none)\n");
    printf("cancel_spans_any_band=%d\n", nany);
    printf("flagged_windows=%ld\n", R->flagged_windows_total);
    printf("analysis_ms=%.2f\n", ms_analysis);
    printf("CANCEL_DETECTED=%d\n", nany > 0 ? 1 : 0);
}

/* --------------------------------------------------------------- selftest */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
static void fill_noise(double* x, int n, Rng* r) {
    for (int i = 0; i < n; i++) x[i] = rng_gauss(r);
}

static void selftest(void) {
    Rng rng; rng_seed(&rng, 0x9E3779B97F4A7C15ULL);
    const int FS = 48000;

    BandMeta bm[MAX_BANDS];
    check("bands_mapped", map_bands(1024, FS, bm) == 8);

    /* threshold closed form for a specific band */
    {
        int W = 32;
        double alpha = 0.01;
        double alpha_band = alpha / 8.0;
        int nu = bm[0].bins * W;
        double want = sqrt(1.0 - pow(alpha_band, 1.0 / (nu - 1)));
        Cfg c = { FS, 1024, W, 0, 0, alpha };
        int n = FS * 2;
        double* a = (double*)malloc(sizeof(double) * (size_t)n);
        fill_noise(a, n, &rng);
        Result R;
        check("threshold_analyze", analyze(a, n, a, n, &c, &R) == 0);
        check("threshold_formula", fabs(R.band[0].thr - want) < 1e-12);
        check("threshold_in_range", R.band[0].thr > 0.01 && R.band[0].thr < 0.6);
        check("nu_reported", R.band[0].nu == nu);
        check("expected_fa_reported", R.expected_false_alarms > 0.0);
        free_result(&R); free(a);
    }

    /* exact inversion: rho = -1 in every band, every window flagged */
    {
        int n = FS * 2;
        double* a = (double*)malloc(sizeof(double) * (size_t)n);
        double* b = (double*)malloc(sizeof(double) * (size_t)n);
        fill_noise(a, n, &rng);
        for (int i = 0; i < n; i++) b[i] = -a[i];
        Cfg c = { FS, 1024, 32, 0, 0, 0.01 };
        Result R;
        check("invert_analyze", analyze(a, n, b, n, &c, &R) == 0);
        check("invert_corr_zero", fabs(R.corr_zero + 1.0) < 1e-9);
        check("invert_rho_neg", R.band[0].rho_total < -0.999);
        check("invert_all_windows_flagged", R.flagged_windows_total == R.nwin * R.n_active);
        check("invert_gain_floor", R.band[0].g_total_db <= -59.0);
        free_result(&R); free(a); free(b);
    }

    /* identical takes: rho = +1, +3.01 dB, nothing flagged */
    {
        int n = FS * 2;
        double* a = (double*)malloc(sizeof(double) * (size_t)n);
        fill_noise(a, n, &rng);
        Cfg c = { FS, 1024, 32, 0, 0, 0.01 };
        Result R;
        check("same_analyze", analyze(a, n, a, n, &c, &R) == 0);
        check("same_corr_zero", fabs(R.corr_zero - 1.0) < 1e-12);
        check("same_rho_pos", fabs(R.band[0].rho_total - 1.0) < 1e-9);
        check("same_gain_3db", fabs(R.band[0].g_total_db - 3.0103) < 0.01);
        check("same_not_flagged", R.flagged_windows_total == 0);
        free_result(&R); free(a);
    }

    /* the core relation: P_sum/(P_a+P_b) == 1 + rho, broadband */
    {
        int n = FS * 2;
        double *a = (double*)malloc(sizeof(double) * (size_t)n);
        double *b = (double*)malloc(sizeof(double) * (size_t)n);
        fill_noise(a, n, &rng);
        for (int i = 0; i < n; i++)
            b[i] = (i >= 3 ? 0.6 * a[i - 3] : 0.0) + 0.8 * rng_gauss(&rng);
        Cfg c = { FS, 1024, 32, 0, 0, 0.01 };
        Result R;
        check("relation_analyze", analyze(a, n, b, n, &c, &R) == 0);
        int used = R.nseg * R.nfft;
        double pa_ = 0, pb_ = 0, ps = 0;
        for (int i = 0; i < used; i++) {
            pa_ += a[i]*a[i]; pb_ += b[i]*b[i]; ps += (a[i]+b[i])*(a[i]+b[i]);
        }
        /* broadband rho computed from the same frames */
        double cr = 0, pa2 = 0, pb2 = 0;
        for (int si = 0; si < R.nseg; si++) {
            cr += R.frame[si].cr_t; pa2 += R.frame[si].pa_t; pb2 += R.frame[si].pb_t;
        }
        double rho_bb = cr / sqrt(pa2 * pb2);
        check("gain_equals_one_plus_rho", fabs(ps / (pa_ + pb_) - (1.0 + rho_bb)) < 5e-3);
        check("band_gain_matches_rho",
              fabs(R.band[0].g_total_db - gain_db(R.band[0].rho_total)) < 1e-9);
        free_result(&R); free(a); free(b);
    }

    /* delay sign against a brute-force reference */
    {
        int n = 8192, d = 7;
        double *a = (double*)calloc((size_t)n, sizeof(double));
        double *b = (double*)calloc((size_t)n, sizeof(double));
        for (int i = 0; i < n; i++) a[i] = rng_gauss(&rng);
        for (int i = d; i < n; i++) b[i] = a[i - d];
        int bf_best = 0; double bf_val = -2.0;
        for (int m = -64; m <= 64; m++) {
            double sab = 0, sa = 0, sb = 0;
            for (int i = 0; i < n; i++) {
                int j = i + m;
                if (j < 0 || j >= n) continue;
                sab += a[i] * b[j]; sa += a[i]*a[i]; sb += b[j]*b[j];
            }
            double r = sab / sqrt(sa * sb);
            if (r > bf_val) { bf_val = r; bf_best = m; }
        }
        Cfg c = { FS, 512, 0, 0, 64, 0.01 };
        Result R;
        check("delay_analyze", analyze(a, n, b, n, &c, &R) == 0);
        check("delay_sign_and_value", R.best_lag == d);
        check("delay_matches_bruteforce", R.best_lag == bf_best);
        check("xcorr_table_present", R.xcorr != NULL && R.max_lag == 64);
        free_result(&R); free(a); free(b);
    }

    /* silence must not crash or flag */
    {
        int n = FS;
        double* a = (double*)calloc((size_t)n, sizeof(double));
        Cfg c = { FS, 1024, 32, 0, 0, 0.01 };
        Result R;
        check("silence_analyze", analyze(a, n, a, n, &c, &R) == 0);
        check("silence_no_flag", R.flagged_windows_total == 0);
        free_result(&R); free(a);
    }

    /* detection power: strong partial cancellation caught, mild not */
    {
        int n = FS * 3;
        double *a = (double*)malloc(sizeof(double) * (size_t)n);
        double *b = (double*)malloc(sizeof(double) * (size_t)n);
        Cfg c = { FS, 1024, 32, 0, 0, 0.01 };
        Rng r2; rng_seed(&r2, 777);
        fill_noise(a, n, &r2);
        /* rho = -0.9: b = -0.9 a + sqrt(0.19) n  =>  C=-0.9P, Pb=0.81P+0.19P */
        for (int i = 0; i < n; i++) b[i] = -0.9 * a[i] + sqrt(0.19) * rng_gauss(&r2);
        Result R;
        check("partial_analyze", analyze(a, n, b, n, &c, &R) == 0);
        check("partial_rho_near_m90", R.band[0].rho_total < -0.85);
        check("partial_detected", R.flagged_windows_total > 0);
        free_result(&R);

        fill_noise(a, n, &r2);
        for (int i = 0; i < n; i++) b[i] = -0.02 * a[i] + sqrt(0.9996) * rng_gauss(&r2);
        check("mild_analyze", analyze(a, n, b, n, &c, &R) == 0);
        /* a single clip is not a fair test: with many windows even a tiny
         * negative rho triggers occasionally.  Measure the rate instead. */
        free_result(&R);
        int mild_windows = 0, mild_hits = 0;
        for (int t = 0; t < 20; t++) {
            fill_noise(a, n, &r2);
            for (int i = 0; i < n; i++)
                b[i] = -0.02 * a[i] + sqrt(0.9996) * rng_gauss(&r2);
            if (analyze(a, n, b, n, &c, &R) == 0) {
                mild_windows += R.nwin;
                mild_hits += (R.flagged_windows_total > 0) ? 1 : 0;
                free_result(&R);
            }
        }
        printf("mild_flag_rate=%.3f\n", mild_windows ? (double)mild_hits / mild_windows : 1.0);
        check("mild_flag_rate_low", mild_windows > 0 &&
              (double)mild_hits / mild_windows <= 0.25);
        free(a); free(b);
    }

    /* empirical false-alarm rate, noise-like material */
    {
        const int TRIALS = 120, n = FS * 2;
        double *a = (double*)malloc(sizeof(double) * (size_t)n);
        double *b = (double*)malloc(sizeof(double) * (size_t)n);
        Cfg c = { FS, 1024, 32, 0, 0, 0.01 };
        long tested = 0, flagged = 0;
        double nominal = 0.0;
        for (int tr = 0; tr < TRIALS; tr++) {
            fill_noise(a, n, &rng);
            fill_noise(b, n, &rng);
            Result R;
            if (analyze(a, n, b, n, &c, &R) == 0) {
                tested += (long)R.nwin * R.n_active;
                flagged += R.flagged_windows_total;
                nominal = R.alpha_band;
                free_result(&R);
            }
        }
        double empirical = tested ? (double)flagged / (double)tested : 1.0;
        printf("false_alarm_nominal=%.6f\n", nominal);
        printf("false_alarm_empirical=%.6f\n", empirical);
        printf("false_alarm_tests=%ld\n", tested);
        check("false_alarm_within_bound", empirical <= nominal + 0.002);
        free(a); free(b);
    }
}

/* ------------------------------------------------------------------- calib */

/* Measure false-alarm rate over independent take pairs, for noise-like and
 * tonal material, sweeping W and alpha.  Tonal material violates the
 * independence assumption: a stationary tone keeps its phase across frames,
 * so nu_eff is far below B*W and the rate rises.  Reporting the gap is the
 * point. */
typedef struct { double band_nominal, band_empirical, win_empirical, nominal_any; } CalibRow;

static void synth_noise_pair(double* a, double* b, int n, Rng* r) {
    fill_noise(a, n, r);
    fill_noise(b, n, r);
}
static void synth_tonal_pair(double* a, double* b, int n, int fs, Rng* r) {
    memset(a, 0, sizeof(double) * (size_t)n);
    memset(b, 0, sizeof(double) * (size_t)n);
    /* several stationary tones, independent random phases per take */
    for (int s = 0; s < 5; s++) {
        double f = 200.0 * pow(2.0, s);
        double wa = 2.0 * PC_PI * f / fs, wb = 2.0 * PC_PI * f / fs;
        double pa_ = rng_uniform(r) * 2.0 * PC_PI;
        double pb_ = rng_uniform(r) * 2.0 * PC_PI;
        for (int i = 0; i < n; i++) {
            a[i] += 0.2 * sin(pa_ + wa * i) + 0.01 * rng_gauss(r);
            b[i] += 0.2 * sin(pb_ + wb * i) + 0.01 * rng_gauss(r);
        }
    }
}

static CalibRow calib_run(int tonal, int W, double alpha, int trials, int fs) {
    CalibRow row = { 0, 0, 0, alpha };
    int n = fs * 2;
    double* a = (double*)malloc(sizeof(double) * (size_t)n);
    double* b = (double*)malloc(sizeof(double) * (size_t)n);
    if (!a || !b) { free(a); free(b); return row; }
    Rng r; rng_seed(&r, 424242 + (uint64_t)(tonal ? 1 : 0) + (uint64_t)W);
    long tested = 0, flagged = 0, wins = 0, win_hits = 0;
    Cfg c = { fs, 1024, W, 0, 0, alpha };
    for (int t = 0; t < trials; t++) {
        if (tonal) synth_tonal_pair(a, b, n, fs, &r);
        else synth_noise_pair(a, b, n, &r);
        Result R;
        if (analyze(a, n, b, n, &c, &R) == 0) {
            tested += (long)R.nwin * R.n_active;
            flagged += R.flagged_windows_total;
            wins += R.nwin;
            win_hits += (R.flagged_windows_total > 0) ? 1 : 0;
            row.band_nominal = R.alpha_band;
            free_result(&R);
        }
    }
    row.band_empirical = tested ? (double)flagged / (double)tested : 0.0;
    row.win_empirical = wins ? (double)win_hits / (double)wins : 0.0;
    free(a); free(b);
    return row;
}

/* Detection power.  Both takes are white noise; take B reuses take A's
 * spectrum *only inside one octave band*, with signed coherence rho0, and is
 * independent everywhere else.  This injects exactly the per-band effect the
 * detector tests, so the measured rate is the detector's power for that rho. */
static void run_power(int trials, int W, double alpha) {
    int fs = 48000, n = fs * 2;
    int N = 2048;
    double* a = (double*)malloc(sizeof(double) * (size_t)n);
    double* b = (double*)malloc(sizeof(double) * (size_t)n);
    if (!a || !b) { free(a); free(b); return; }
    double rhos[] = { -0.1, -0.2, -0.3, -0.4, -0.5, -0.7, -0.9 };
    Cfg c = { fs, 1024, W, 0, 0, alpha };
    printf("# detection power, W=%d alpha=%.2f, %d trials per cell\n", W, alpha, trials);
    printf("# target_hz,rho_true,band_thr,mean_rho,detect_rate\n");
    int targets[] = { 500, 1000, 4000 };
    for (int ti = 0; ti < 3; ti++) {
        int fc = targets[ti];
        int lo = (int)(fc / PC_SQRT2 * N / fs), hi = (int)(fc * PC_SQRT2 * N / fs);
        for (int ri = 0; ri < 7; ri++) {
            double rho0 = rhos[ri];
            Rng r; rng_seed(&r, 9000 + (uint64_t)ri * 31 + (uint64_t)fc);
            long hits = 0, ok = 0;
            double thr = 0, sum_rho = 0;
            for (int t = 0; t < trials; t++) {
                /* one frame's worth of independent spectra, tiled */
                D* fa = (D*)calloc((size_t)N, sizeof(D));
                D* fb = (D*)calloc((size_t)N, sizeof(D));
                for (int i = 0; i < N; i++) {
                    fa[i].re = rng_gauss(&r); fa[i].im = rng_gauss(&r);
                }
                for (int k = 0; k < N / 2; k++) {
                    if (k >= lo && k <= hi) {
                        /* b's bin = rho0*a's bin + sqrt(1-rho0^2)*Z */
                        double s = sqrt(fabs(1.0 - rho0 * rho0));
                        fb[k].re = rho0 * fa[k].re + s * rng_gauss(&r);
                        fb[k].im = rho0 * fa[k].im + s * rng_gauss(&r);
                    } else {
                        fb[k].re = rng_gauss(&r); fb[k].im = rng_gauss(&r);
                    }
                }
                /* real time-domain signals: take Re of inverse FFT of a
                 * Hermitian-symmetric spectrum */
                for (int k = 1; k < N / 2; k++) {
                    fa[N - k].re = fa[k].re; fa[N - k].im = -fa[k].im;
                    fb[N - k].re = fb[k].re; fb[N - k].im = -fb[k].im;
                }
                fa[0].im = fa[N/2].im = fb[0].im = fb[N/2].im = 0;
                fft_radix2(fa, N, 1);
                fft_radix2(fb, N, 1);
                int reps = (n + N - 1) / N;
                for (int i = 0, p = 0; p < reps; p++)
                    for (int j = 0; j < N && i < n; j++, i++) {
                        a[i] = fa[j].re; b[i] = fb[j].re;
                    }
                free(fa); free(fb);
                Result R;
                if (analyze(a, n, b, n, &c, &R) == 0) {
                    hits += (R.flagged_windows_total > 0);
                    ok++;
                    /* find the BandMeta for this centre */
                    for (int bi = 0; bi < MAX_BANDS; bi++) {
                        if (R.band[bi].active && fabs(BAND_CENTER[bi] - fc) < 1.0) {
                            thr = R.band[bi].thr;
                            sum_rho += R.band[bi].rho_total;
                        }
                    }
                    free_result(&R);
                }
            }
            printf("%d,%.2f,%.4f,%.4f,%.4f\n", fc, rho0, thr,
                   ok ? sum_rho / ok : 0.0, ok ? (double)hits / ok : 0.0);
        }
    }
    free(a); free(b);
}

static void run_calib(int trials) {
    int Wlist[] = { 16, 32, 64 };
    double alist[] = { 0.01, 0.05, 0.10 };
    printf("# false-alarm calibration: %d independent take pairs per cell\n", trials);
    printf("# material,W,alpha,band_nominal,band_measured,band_ratio,"
           "window_nominal(window_hits>0),window_measured\n");
    for (int mat = 0; mat < 2; mat++) {
        for (int wi = 0; wi < 3; wi++) {
            for (int ai = 0; ai < 3; ai++) {
                CalibRow r = calib_run(mat, Wlist[wi], alist[ai], trials, 48000);
                printf("%s,%d,%.2f,%.6f,%.6f,%.3f,%.6f,%.6f\n",
                       mat ? "tonal" : "noise", Wlist[wi], alist[ai],
                       r.band_nominal, r.band_empirical,
                       r.band_nominal > 0 ? r.band_empirical / r.band_nominal : 0.0,
                       r.nominal_any, r.win_empirical);
            }
        }
    }
    printf("\n");
    run_power(trials, 32, 0.01);
}

/* ------------------------------------------------------------------- bench */

/* Build material from `nsources` independent sources.  Take B is the exact
 * inverse of take A, so the cancellation structure (rho = -1 in every band)
 * is identical for both benchmark rows; only the source count differs.  The
 * analysis cost is governed by duration and FFT size, not by source count,
 * so this isolates the thing the benchmark is meant to measure. */
static double synth_material(double* a, double* b, int n, int nsources, int fs, Rng* r) {
    double t0 = now_ms();
    memset(a, 0, sizeof(double) * (size_t)n);
    double scale = 1.0 / sqrt((double)nsources);
    for (int s = 0; s < nsources; s++) {
        double f = 110.0 * pow(2.0, (double)(s % 12) * 0.4);
        double ph = 2.0 * PC_PI * rng_uniform(r);
        double w = 2.0 * PC_PI * f / fs;
        for (int i = 0; i < n; i++)
            a[i] += scale * (0.4 * sin(ph + w * i) + 0.6 * rng_gauss(r));
    }
    for (int i = 0; i < n; i++) b[i] = -a[i];
    return now_ms() - t0;
}

static int run_bench(int n1, int n2, double secs, int fs, const Cfg* base) {
    int counts[2] = { n1, n2 };
    for (int ci = 0; ci < 2; ci++) {
        int ns = counts[ci] < 1 ? 1 : counts[ci];
        int n = (int)(secs * fs);
        double* a = (double*)malloc(sizeof(double) * (size_t)n);
        double* b = (double*)malloc(sizeof(double) * (size_t)n);
        if (!a || !b) { free(a); free(b); return -1; }
        Rng r; rng_seed(&r, 12345 + (uint64_t)ns);
        double synth_ms = synth_material(a, b, n, ns, fs, &r);
        Cfg c = *base;
        c.fs = fs;
        double t0 = now_ms();
        Result R;
        int rc = analyze(a, n, b, n, &c, &R);
        double ms = now_ms() - t0;
        if (rc == 0) {
            printf("bench_sources=%d  seconds=%.1f  samples=%d  synth_ms=%.2f  "
                   "analysis_ms=%.2f  realtime_x=%.1f  flagged_windows=%ld\n",
                   ns, secs, n, synth_ms, ms, ms > 0 ? (secs * 1000.0 / ms) : 0.0,
                   R.flagged_windows_total);
            free_result(&R);
        } else {
            printf("bench_sources=%d  seconds=%.1f  analysis_ms=%.2f  status=failed\n",
                   ns, secs, ms);
        }
        free(a); free(b);
    }
    return 0;
}

/* -------------------------------------------------------------------- main */

static void usage(const char* argv0) {
    fprintf(stderr,
        "usage: %s [options] A.wav B.wav [ch_a ch_b]\n"
        "  --fs N        sample rate when input is raw (WAV header wins)\n"
        "  --seg N       segment length in samples, rounded down to pow2 (default 1024)\n"
        "  --win N       frames aggregated per detection test (default 32)\n"
        "  --hop N       window hop in frames (default win/2)\n"
        "  --max-lag N   max |delay| searched, samples (default 4096)\n"
        "  --alpha X     per-window across-bands false-alarm rate (default 0.01)\n"
        "  --accept-cancellation  allow --out to write when cancellation is found\n"
        "  --out FILE    write A+B to FILE; refused if cancellation is detected\n"
        "  --calib N     measure false-alarm rates over N take pairs and exit\n"
        "  --bench N1,N2 time analysis for two source counts (default 2,8)\n"
        "  --bench-sec X seconds of synthetic material per benchmark (default 10)\n"
        "  --selftest    run built-in checks and exit\n", argv0);
}

int main(int argc, char** argv) {
    Cfg cfg = { 0, 1024, 32, 0, 4096, 0.01 };
    int accept = 0, selftest_flag = 0, bench_flag = 0, calib_trials = 0;
    int ch_a = 0, ch_b = 0;
    int bench_n1 = 2, bench_n2 = 8;
    double bench_sec = 10.0;
    const char* out = NULL;
    const char* files[2] = { NULL, NULL };
    int nfiles = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--selftest")) selftest_flag = 1;
        else if (!strcmp(argv[i], "--accept-cancellation")) accept = 1;
        else if (!strcmp(argv[i], "--calib")) {
            calib_trials = 200;
            if (i + 1 < argc && argv[i + 1][0] != '-') calib_trials = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "--bench")) {
            bench_flag = 1;
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                int a = 0, b = 0;
                if (sscanf(argv[++i], "%d,%d", &a, &b) == 2) { bench_n1 = a; bench_n2 = b; }
            }
        }
        else if (!strcmp(argv[i], "--bench-sec") && i + 1 < argc) bench_sec = atof(argv[++i]);
        else if (!strcmp(argv[i], "--fs") && i + 1 < argc) cfg.fs = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seg") && i + 1 < argc) cfg.seg_len = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--win") && i + 1 < argc) cfg.win = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--hop") && i + 1 < argc) cfg.hop = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--max-lag") && i + 1 < argc) cfg.max_lag = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--alpha") && i + 1 < argc) cfg.alpha = atof(argv[++i]);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) { usage(argv[0]); return 0; }
        else if (argv[i][0] != '-' && nfiles < 2) files[nfiles++] = argv[i];
        else if (argv[i][0] != '-' && nfiles == 2) {
            if (ch_a == 0) ch_a = atoi(argv[i]); else ch_b = atoi(argv[i]);
        } else { usage(argv[0]); return 2; }
    }

    if (selftest_flag) { selftest(); printf("all_ok=%d\n", failures == 0); return failures ? 1 : 0; }
    if (calib_trials > 0) { run_calib(calib_trials); return 0; }
    if (bench_flag) {
        if (!cfg.fs) cfg.fs = 48000;
        return run_bench(bench_n1, bench_n2, bench_sec, cfg.fs, &cfg) == 0 ? 0 : 1;
    }
    if (nfiles < 2) { usage(argv[0]); return 2; }

    Wav wa, wb;
    if (load_wav(files[0], ch_a, &wa)) { fprintf(stderr, "error: %s\n", wa.err); return 1; }
    if (load_wav(files[1], ch_b, &wb)) {
        fprintf(stderr, "error: %s\n", wb.err); free(wa.x); return 1;
    }
    if (!cfg.fs) cfg.fs = wa.fs;
    if (wa.fs != wb.fs) {
        fprintf(stderr, "error: sample-rate mismatch %d vs %d\n", wa.fs, wb.fs);
        free(wa.x); free(wb.x); return 1;
    }

    double t0 = now_ms();
    Result R;
    int rc = analyze(wa.x, wa.n, wb.x, wb.n, &cfg, &R);
    double ms = now_ms() - t0;
    if (rc) {
        fprintf(stderr, "error: inputs unusable (need >= 2048 samples and >= 4 frames)\n");
        free(wa.x); free(wb.x); return 1;
    }
    print_report(&R, files[0], files[1], ms);

    /* ---- alignment re-check ----
     * A pure delay makes low bands look negatively coherent (comb filtering)
     * even when the takes are not polarity-inverted.  If a nonzero lag was
     * found and anything was flagged, re-run after shifting B back by that
     * lag.  If the flags disappear, the cause was delay/drift, not inversion;
     * if they survive, it is genuine cancellation.  Either way it is printed. */
    if (R.best_lag != 0 && R.flagged_windows_total > 0) {
        int L2 = R.L;
        int lag = R.best_lag;
        double* bshift = (double*)calloc((size_t)L2, sizeof(double));
        if (bshift) {
            /* align B to A: B'[i] = B[i + lag] (lag > 0 means B lags A) */
            for (int i = 0; i < L2; i++) {
                int j = i + lag;
                bshift[i] = (j >= 0 && j < wb.n) ? wb.x[j] : 0.0;
            }
            Result Ra;
            Cfg c2 = cfg;
            if (analyze(wa.x, wa.n, bshift, L2, &c2, &Ra) == 0) {
                printf("\n[alignment-recheck]\n");
                printf("shifted_b_by_samples=%d\n", -lag);
                printf("flagged_windows_before=%ld\n", R.flagged_windows_total);
                printf("flagged_windows_after=%ld\n", Ra.flagged_windows_total);
                int delay_cause = R.flagged_windows_total > 0 && Ra.flagged_windows_total == 0;
                printf("delay_explains_cancellation=%d\n", delay_cause ? 1 : 0);
                printf("genuine_inversion_after_align=%d\n",
                       Ra.flagged_windows_total > 0 ? 1 : 0);
                free_result(&Ra);
            }
            free(bshift);
        }
    }


    int cancel = R.flagged_windows_total > 0;
    int refused = 0;
    if (out) {
        if (cancel && !accept) {
            printf("mix_written=0\nmix_refused_reason=cancellation_detected\n");
            refused = 1;
        } else {
            int n = R.nseg * R.nfft;
            float* mix = (float*)malloc(sizeof(float) * (size_t)n);
            if (!mix) fprintf(stderr, "error: out of memory\n");
            else {
                for (int i = 0; i < n; i++) mix[i] = (float)(wa.x[i] + wb.x[i]);
                int clipped = 0;
                if (write_wav16(out, mix, n, cfg.fs, &clipped))
                    fprintf(stderr, "error: cannot write %s\n", out);
                else
                    printf("mix_written=1\nmix_clipped_samples=%d\n", clipped);
                free(mix);
            }
        }
    }

    free_result(&R); free(wa.x); free(wb.x);
    if (cancel) fprintf(stderr, "phasecheck: CANCELLATION DETECTED - see report\n");
    if (refused) return 4;
    return cancel ? 3 : 0;
}
