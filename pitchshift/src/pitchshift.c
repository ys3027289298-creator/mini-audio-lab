
/* pitchshift: two-path pitch shift for footsteps.
 *
 * A whole-sample varispeed shift scales every frequency by the ratio,
 * including the low thud that gives the step its weight, and shortens the
 * amplitude envelope by 1/ratio.
 *
 * Here the shift is split in two at a crossover frequency:
 *   low path  - the band below the crossover, kept as-is, so the heel
 *               strike keeps its pitch and its envelope duration;
 *   high path - the band above the crossover, time-stretched by the ratio
 *               then resampled by the ratio: pitch scales, duration does not.
 * The crossover is the source's spectral energy centroid, held fixed across
 * speeds.  Re-deriving it from the processed material is the failure mode:
 * the centroid rides up with the pitch and the low band collapses.
 *
 * Usage: ./pitchshift [file.wav]
 *   With no argument a deterministic synthetic footstep is used.
 *   WAV input: mono or multi-channel (averaged), PCM16 or float32.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define PI 3.14159265358979323846
#define NSPEED 3
#define WIN 1024
#define HOP 256
#define FMIN 20.0
#define FMAX 20000.0

static int g_sr = 48000;

static const double RATIO[NSPEED] = { 1.00, 1.25, 1.50 };
static const char* NAME[NSPEED] = { "slow", "normal", "fast" };

static void* xcalloc(size_t n, size_t sz) {
    void* p = calloc(n, sz);
    if (p == NULL) { fprintf(stderr, "out of memory\n"); exit(2); }
    return p;
}

/* ---------- FFT (radix-2, double) ---------- */

static int fft_ok(int n) { return n >= 2 && (n & (n - 1)) == 0; }

static void fft(double* re, double* im, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * PI / (double)len * (inverse ? 1.0 : -1.0);
        double wr = cos(ang), wi = sin(ang);
        for (int i = 0; i < n; i += len) {
            double cwr = 1.0, cwi = 0.0;
            for (int k = 0; k < len / 2; k++) {
                int a = i + k, b = a + len / 2;
                double vr = re[b] * cwr - im[b] * cwi;
                double vi = re[b] * cwi + im[b] * cwr;
                re[b] = re[a] - vr; im[b] = im[a] - vi;
                re[a] += vr;         im[a] += vi;
                double nwr = cwr * wr - cwi * wi;
                double nwi = cwr * wi + cwi * wr;
                cwr = nwr; cwi = nwi;
            }
        }
    }
    if (inverse) {
        for (int i = 0; i < n; i++) { re[i] /= (double)n; im[i] /= (double)n; }
    }
}

/* ---------- spectral measures ---------- */

/* Energy centroid and the share of energy below `edge`, both restricted to
 * 20 Hz .. 20 kHz, Hann-windowed and zero-padded to a power of two. */
static void spectral_stats(const double* x, int n, double edge,
                           double* centroid, double* share) {
    if (n < 2) { *centroid = 0.0; *share = 0.0; return; }
    int nfft = 4096;
    while (nfft < n) nfft <<= 1;
    double* re = xcalloc((size_t)nfft, sizeof(double));
    double* im = xcalloc((size_t)nfft, sizeof(double));
    for (int i = 0; i < n; i++) {
        double w = 0.5 - 0.5 * cos(2.0 * PI * (double)i / (double)(n - 1));
        re[i] = x[i] * w;
    }
    fft(re, im, nfft, 0);
    double sumP = 0.0, sumfP = 0.0, sumLow = 0.0;
    for (int k = 1; k < nfft / 2; k++) {
        double f = (double)k * (double)g_sr / (double)nfft;
        if (f < FMIN || f > FMAX) continue;
        double P = re[k] * re[k] + im[k] * im[k];
        sumP += P;
        sumfP += f * P;
        if (f < edge) sumLow += P;
    }
    *centroid = sumP > 0.0 ? sumfP / sumP : 0.0;
    *share = sumP > 0.0 ? sumLow / sumP : 0.0;
    free(re); free(im);
}

/* ---------- complementary crossover split ---------- */

/* Low-path gain: 1 below fc*2^-1/3, 0 above fc*2^1/3, cos^2 ramp between.
 * High-path gain is 1 - low, so low + high reconstructs the input exactly and
 * the -6 dB point sits on fc.  Zero-phase, applied in the frequency domain. */
static double low_gain(double f, double fc) {
    double a = fc * pow(2.0, -1.0 / 3.0);
    double b = fc * pow(2.0, 1.0 / 3.0);
    if (f >= b) return 0.0;
    if (f <= a) return 1.0;
    double u = (f - a) / (b - a);
    return 0.5 + 0.5 * cos(PI * u);
}

static void split_bands(const double* x, int n, double fc, double* low, double* high) {
    int nfft = 4096;
    while (nfft < 2 * n) nfft <<= 1;
    double* lr = xcalloc((size_t)nfft, sizeof(double));
    double* li = xcalloc((size_t)nfft, sizeof(double));
    double* hr = xcalloc((size_t)nfft, sizeof(double));
    double* hi = xcalloc((size_t)nfft, sizeof(double));
    memcpy(lr, x, (size_t)n * sizeof(double));
    fft(lr, li, nfft, 0);
    memcpy(hr, lr, (size_t)nfft * sizeof(double));
    memcpy(hi, li, (size_t)nfft * sizeof(double));
    for (int k = 0; k < nfft; k++) {
        int bin = k <= nfft / 2 ? k : nfft - k;   /* Hermitian mirror */
        double f = (double)bin * (double)g_sr / (double)nfft;
        double g = low_gain(f, fc);
        lr[k] *= g;          li[k] *= g;
        hr[k] *= (1.0 - g);  hi[k] *= (1.0 - g);
    }
    fft(lr, li, nfft, 1);
    fft(hr, hi, nfft, 1);
    memcpy(low, lr, (size_t)n * sizeof(double));
    memcpy(high, hr, (size_t)n * sizeof(double));
    free(lr); free(li); free(hr); free(hi);
}

/* ---------- time stretch / resample ---------- */

static double ola_window(int j) {
    /* non-zero at both ends so the windowed overlap-add sum never hits zero */
    return 0.5 - 0.5 * cos(2.0 * PI * ((double)j + 0.5) / (double)WIN);
}

/* Overlap-add time stretch by `ratio` (ratio > 1 makes the signal longer),
 * synthesis hop HOP, analysis hop HOP/ratio.  Returns output length.
 * At ratio == 1 this is an exact reconstruction. */
static int ola_stretch(const double* x, int n, double ratio, double* y, int cap) {
    double ha = (double)HOP / ratio;
    int nout = (int)((double)n / ha) * HOP + WIN;
    if (nout > cap) nout = cap;
    memset(y, 0, (size_t)nout * sizeof(double));
    double* wsum = xcalloc((size_t)nout, sizeof(double));
    for (int k = 0; ; k++) {
        double pos = (double)k * ha;
        int i0 = (int)floor(pos);
        double t = pos - (double)i0;
        if (i0 + WIN + 1 >= n) break;
        int o = k * HOP;
        if (o + WIN > nout) break;
        for (int j = 0; j < WIN; j++) {
            double s = (1.0 - t) * x[i0 + j] + t * x[i0 + j + 1];
            double w = ola_window(j);
            y[o + j] += s * w;
            wsum[o + j] += w;
        }
    }
    for (int i = 0; i < nout; i++)
        if (wsum[i] > 1e-9) y[i] /= wsum[i];
        else y[i] = 0.0;
    free(wsum);
    return nout;
}

static double interp_at(const double* x, int n, double pos) {
    if (n <= 0) return 0.0;
    if (pos < 0.0) pos = 0.0;
    int i0 = (int)floor(pos);
    if (i0 >= n - 1) return x[n - 1];
    double t = pos - (double)i0;
    return (1.0 - t) * x[i0] + t * x[i0 + 1];
}

/* Pitch shift that keeps duration: stretch by ratio, then resample by ratio. */
static int pitch_shift_hold_duration(const double* x, int n, double ratio, double* out) {
    int np = n + WIN + 2;
    double* xp = xcalloc((size_t)np, sizeof(double));
    memcpy(xp, x, (size_t)n * sizeof(double));
    int cap = (int)((double)np * ratio) + 2 * WIN + 2 * HOP;
    double* ys = xcalloc((size_t)cap, sizeof(double));
    int m = ola_stretch(xp, np, ratio, ys, cap);
    for (int i = 0; i < n; i++)
        out[i] = interp_at(ys, m, (double)i * ratio);
    free(xp); free(ys);
    return n;
}

/* Whole-sample varispeed: pitch and tempo both scale by ratio. */
static int varispeed(const double* x, int n, double ratio, double* out, int cap) {
    int m = (int)floor((double)n / ratio);
    if (m > cap) m = cap;
    for (int i = 0; i < m; i++)
        out[i] = interp_at(x, n, (double)i * ratio);
    return m;
}

/* ---------- synthetic footstep ---------- */

typedef struct { double b0, b1, b2, a1, a2; } Biquad;

static Biquad bq_make(int kind, double fc, double q) {
    double w0 = 2.0 * PI * fc / (double)g_sr;
    double cw = cos(w0), sw = sin(w0);
    double alpha = sw / (2.0 * q);
    double a0 = 1.0 + alpha;
    Biquad b;
    if (kind == 0) {          /* lowpass  */
        b.b0 = (1.0 - cw) * 0.5; b.b1 = 1.0 - cw; b.b2 = (1.0 - cw) * 0.5;
    } else if (kind == 1) {   /* highpass */
        b.b0 = (1.0 + cw) * 0.5; b.b1 = -(1.0 + cw); b.b2 = (1.0 + cw) * 0.5;
    } else {                  /* bandpass */
        b.b0 = alpha; b.b1 = 0.0; b.b2 = -alpha;
    }
    b.b0 /= a0; b.b1 /= a0; b.b2 /= a0;
    b.a1 = -2.0 * cw / a0; b.a2 = (1.0 - alpha) / a0;
    return b;
}

static double bq_run(const Biquad* b, double x, double* z1, double* z2) {
    double y = b->b0 * x + *z1;
    *z1 = b->b1 * x - b->a1 * y + *z2;
    *z2 = b->b2 * x - b->a2 * y;
    return y;
}

static void filter_chain(double* x, int n, Biquad b, int stages) {
    for (int s = 0; s < stages; s++) {
        double z1 = 0.0, z2 = 0.0;
        for (int i = 0; i < n; i++) x[i] = bq_run(&b, x[i], &z1, &z2);
    }
}

/* Heel strike: a low thud with a short downward glide (the weight), a low-mid
 * body, a bright click and a short scuff.  Deterministic, so every run is
 * comparable. */
static int synth_footstep(double** out) {
    int n = (int)(0.40 * (double)g_sr);
    double* w = xcalloc((size_t)n, sizeof(double));
    unsigned s = 12345u;
    for (int i = 0; i < n; i++) {
        s = s * 1664525u + 1013904223u;
        w[i] = (double)s / 2147483648.0 - 1.0;
    }
    double* body = xcalloc((size_t)n, sizeof(double));
    double* click = xcalloc((size_t)n, sizeof(double));
    double* scuff = xcalloc((size_t)n, sizeof(double));
    memcpy(body, w, (size_t)n * sizeof(double));
    memcpy(click, w, (size_t)n * sizeof(double));
    memcpy(scuff, w, (size_t)n * sizeof(double));
    filter_chain(body, n, bq_make(2, 280.0, 0.9), 2);
    filter_chain(click, n, bq_make(1, 3000.0, 0.7071), 2);
    filter_chain(scuff, n, bq_make(2, 1200.0, 1.0), 2);

    double* x = xcalloc((size_t)n, sizeof(double));
    double ph = 0.0;
    int delay = (int)(0.035 * (double)g_sr);
    for (int i = 0; i < n; i++) {
        double t = (double)i / (double)g_sr;
        double f = 52.0 + 10.0 * exp(-t / 0.030);
        ph += 2.0 * PI * f / (double)g_sr;
        double thud = 0.85 * exp(-t / 0.075) * sin(ph);
        double bd = 0.55 * exp(-t / 0.045) * body[i];
        double ck = 0.22 * exp(-t / 0.006) * click[i];
        double sc = i >= delay ? 0.12 * exp(-t / 0.060) * scuff[i - delay] : 0.0;
        x[i] = thud + bd + ck + sc;
    }
    double peak = 0.0;
    for (int i = 0; i < n; i++) { double a = fabs(x[i]); if (a > peak) peak = a; }
    if (peak > 0.0) for (int i = 0; i < n; i++) x[i] *= 0.90 / peak;
    free(w); free(body); free(click); free(scuff);
    *out = x;
    return n;
}

/* ---------- WAV input ---------- */

static uint32_t rd32(const unsigned char* p) {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}
static uint16_t rd16(const unsigned char* p) {
    return (uint16_t)((uint16_t)p[0] | ((uint16_t)p[1] << 8));
}

/* Average all channels to mono.  Returns 0 on success. */
static int load_wav(const char* path, double** out, int* out_n) {
    FILE* f = fopen(path, "rb");
    if (f == NULL) { fprintf(stderr, "cannot open %s\n", path); return -1; }
    unsigned char hdr[12];
    if (fread(hdr, 1, 12, f) != 12 || memcmp(hdr, "RIFF", 4) != 0 ||
        memcmp(hdr + 8, "WAVE", 4) != 0) {
        fprintf(stderr, "%s: not a RIFF/WAVE file\n", path);
        fclose(f);
        return -1;
    }
    int fmt = 0;
    unsigned channels = 0, bits = 0;
    unsigned long rate = 0;
    long data_off = -1, data_len = 0;
    for (;;) {
        unsigned char ch[8];
        if (fread(ch, 1, 8, f) != 8) break;
        uint32_t sz = rd32(ch + 4);
        if (memcmp(ch, "fmt ", 4) == 0) {
            unsigned char b[16];
            if (sz < 16 || fread(b, 1, 16, f) != 16) break;
            fmt = rd16(b);
            channels = rd16(b + 2);
            rate = rd32(b + 4);
            bits = rd16(b + 14);
            if (sz > 16) fseek(f, (long)(sz - 16), SEEK_CUR);
        } else if (memcmp(ch, "data", 4) == 0) {
            data_off = ftell(f);
            data_len = (long)sz;
            fseek(f, (long)sz + (sz & 1u), SEEK_CUR);
            break;
        } else {
            fseek(f, (long)sz + (sz & 1u), SEEK_CUR);
        }
    }
    if (data_off < 0 || channels == 0 || rate == 0 ||
        !((fmt == 1 && bits == 16) || (fmt == 3 && bits == 32))) {
        fprintf(stderr, "%s: need mono/stereo PCM16 or float32 WAV\n", path);
        fclose(f);
        return -1;
    }
    int bytes = (int)(bits / 8);
    long frames = data_len / (bytes * (int)channels);
    if (frames <= 1) { fprintf(stderr, "%s: empty data chunk\n", path); fclose(f); return -1; }
    unsigned char* raw = xcalloc((size_t)frames * (size_t)channels * (size_t)bytes, 1);
    fseek(f, data_off, SEEK_SET);
    if (fread(raw, (size_t)bytes, (size_t)frames * channels, f) !=
        (size_t)frames * channels) {
        fprintf(stderr, "%s: short read\n", path);
        free(raw);
        fclose(f);
        return -1;
    }
    fclose(f);
    double* x = xcalloc((size_t)frames, sizeof(double));
    for (long i = 0; i < frames; i++) {
        double acc = 0.0;
        for (unsigned c = 0; c < channels; c++) {
            const unsigned char* p = raw + ((size_t)i * channels + c) * (size_t)bytes;
            if (fmt == 1) acc += (double)(int16_t)rd16(p) / 32768.0;
            else { float v; memcpy(&v, p, sizeof(v)); acc += (double)v; }
        }
        x[i] = acc / (double)channels;
    }
    free(raw);
    g_sr = (int)rate;
    *out = x;
    *out_n = (int)frames;
    return 0;
}

/* ---------- measures ---------- */

static double envelope_duration_ms(const double* x, int n) {
    double c = exp(-1.0 / (0.005 * (double)g_sr));   /* 5 ms one-pole */
    double* e = xcalloc((size_t)n, sizeof(double));
    double acc = 0.0, peak = 0.0;
    for (int i = 0; i < n; i++) {
        acc = c * acc + (1.0 - c) * fabs(x[i]);
        e[i] = acc;
        if (acc > peak) peak = acc;
    }
    double thr = 0.05 * peak;
    int first = -1, last = -1;
    for (int i = 0; i < n; i++) {
        if (e[i] > thr) { if (first < 0) first = i; last = i; }
    }
    free(e);
    if (first < 0) return 0.0;
    return (double)(last - first + 1) * 1000.0 / (double)g_sr;
}

typedef struct {
    double ratio;
    double env_ms;
    double low_share;
    double crossover_hz;    /* crossover used, or derived from the output */
    double out_centroid_hz;
    double peak;
    int n_out;
} Group;

static void measure(const double* x, int n, double fc0, Group* g) {
    g->n_out = n;
    g->env_ms = envelope_duration_ms(x, n);
    double cen = 0.0, share = 0.0;
    spectral_stats(x, n, fc0, &cen, &share);
    g->out_centroid_hz = cen;
    g->low_share = share;
    double pk = 0.0;
    for (int i = 0; i < n; i++) { double a = fabs(x[i]); if (a > pk) pk = a; }
    g->peak = pk;
}

/* Output band-limited to the source crossover, so it can be measured and
 * compared with the source on the same terms. */
static void run_anchored(const double* x, int n, double fc0, double ratio, Group* g) {
    double* low = xcalloc((size_t)n, sizeof(double));
    double* high = xcalloc((size_t)n, sizeof(double));
    double* hs = xcalloc((size_t)n, sizeof(double));
    double* out = xcalloc((size_t)n, sizeof(double));
    split_bands(x, n, fc0, low, high);
    pitch_shift_hold_duration(high, n, ratio, hs);
    for (int i = 0; i < n; i++) out[i] = low[i] + hs[i];
    g->ratio = ratio;
    g->crossover_hz = fc0;
    measure(out, n, fc0, g);
    free(low); free(high); free(hs); free(out);
}

/* Whole-sample varispeed, with the crossover a "derive it from the material"
 * pass would pick up on the sped-up result.  The low share is still measured
 * against the source crossover, which is where the weight of the step lives. */
static void run_naive(const double* x, int n, double fc0, double ratio, Group* g) {
    double* out = xcalloc((size_t)n, sizeof(double));
    int m = varispeed(x, n, ratio, out, n);
    g->ratio = ratio;
    measure(out, m, fc0, g);
    g->crossover_hz = g->out_centroid_hz;
    free(out);
}

/* ---------- test harness ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* Structural checks hold for any input.  The material checks below encode
 * what the synthetic reference footstep is expected to do; they are skipped
 * for external WAVs because their thresholds depend on the material. */
static void checks_structural(const double* x, int ns, double fc0, const Group* anch) {
    check("fft_pow2", fft_ok(4096) && !fft_ok(4095));
    int nfft = 64;
    double* re = xcalloc((size_t)nfft, sizeof(double));
    double* im = xcalloc((size_t)nfft, sizeof(double));
    for (int i = 0; i < nfft; i++) re[i] = 1.0;
    fft(re, im, nfft, 0);
    check("fft_dc", fabs(re[0] - (double)nfft) < 1e-9);
    free(re); free(im);

    double* low = xcalloc((size_t)ns, sizeof(double));
    double* high = xcalloc((size_t)ns, sizeof(double));
    split_bands(x, ns, fc0, low, high);
    double maxrec = 0.0;
    for (int i = 0; i < ns; i++) {
        double d = fabs(low[i] + high[i] - x[i]);
        if (d > maxrec) maxrec = d;
    }
    check("split_reconstructs", maxrec < 1e-9);
    free(low); free(high);

    double* low2 = xcalloc((size_t)ns, sizeof(double));
    double* high2 = xcalloc((size_t)ns, sizeof(double));
    double* hs = xcalloc((size_t)ns, sizeof(double));
    split_bands(x, ns, fc0, low2, high2);
    pitch_shift_hold_duration(high2, ns, 1.0, hs);
    double id_err = 0.0;
    for (int i = 0; i < ns; i++) {
        double d = fabs(low2[i] + hs[i] - x[i]);
        if (d > id_err) id_err = d;
    }
    check("anchored_identity_ratio1", id_err < 1e-6);
    free(low2); free(high2); free(hs);

    check("anchored_crossover_fixed", fabs(anch[2].crossover_hz - fc0) < 1e-9);
    check("anchored_crossover_stable", fabs(anch[2].out_centroid_hz / fc0 - 1.0) < 0.10);
}

static void checks_reference_material(const Group* src, const Group* anch, const Group* naive) {
    double fc0 = src->crossover_hz;
    check("anchored_env_held_fast", fabs(anch[2].env_ms / src->env_ms - 1.0) < 0.10);
    check("anchored_low_share_held", anch[2].low_share > 0.60);

    check("naive_low_share_collapses", naive[2].low_share < 0.01);
    check("naive_crossover_drifts", naive[2].crossover_hz > 1.35 * fc0 &&
                                    naive[2].crossover_hz < 1.65 * fc0);
    check("naive_env_compresses", fabs(naive[2].env_ms / (src->env_ms / 1.5) - 1.0) < 0.20);
    check("naive_low_share_by_speed", naive[1].low_share < naive[0].low_share &&
                                      naive[2].low_share < naive[0].low_share);
}

int main(int argc, char** argv) {
    double* x = NULL;
    int ns = 0;
    int external = 0;
    if (argc > 1) {
        if (load_wav(argv[1], &x, &ns) != 0) return 2;
        printf("input: %s  sr=%d  n=%d\n", argv[1], g_sr, ns);
        external = 1;
    } else {
        ns = synth_footstep(&x);
        printf("input: synthetic footstep  sr=%d  n=%d\n", g_sr, ns);
    }
    if (ns < WIN + 2) { fprintf(stderr, "input too short\n"); free(x); return 2; }

    double fc0 = 0.0, dummy = 0.0;
    spectral_stats(x, ns, 0.0, &fc0, &dummy);
    Group src;
    src.ratio = 1.0;
    src.crossover_hz = fc0;
    measure(x, ns, fc0, &src);

    Group anch[NSPEED], naive[NSPEED];
    for (int i = 0; i < NSPEED; i++) {
        run_anchored(x, ns, fc0, RATIO[i], &anch[i]);
        run_naive(x, ns, fc0, RATIO[i], &naive[i]);
    }

    printf("source: env_ms=%.1f low_share=%.4f centroid_hz=%.1f peak=%.3f\n\n",
           src.env_ms, src.low_share, src.out_centroid_hz, src.peak);

    printf("== anchored two-path (low band held, high band pitch x ratio) ==\n");
    printf("%-7s %6s %9s %10s %13s %16s %7s\n",
           "speed", "ratio", "env_ms", "low_share", "crossover_hz", "out_centroid_hz", "peak");
    for (int i = 0; i < NSPEED; i++)
        printf("%-7s %6.2f %9.1f %10.4f %13.1f %16.1f %7.3f\n",
               NAME[i], anch[i].ratio, anch[i].env_ms, anch[i].low_share,
               anch[i].crossover_hz, anch[i].out_centroid_hz, anch[i].peak);

    printf("\n== naive varispeed (whole-sample proportional scaling) ==\n");
    printf("%-7s %6s %9s %10s %13s %16s %7s\n",
           "speed", "ratio", "env_ms", "low_share", "crossover_hz", "out_centroid_hz", "peak");
    for (int i = 0; i < NSPEED; i++)
        printf("%-7s %6.2f %9.1f %10.4f %13.1f %16.1f %7.3f\n",
               NAME[i], naive[i].ratio, naive[i].env_ms, naive[i].low_share,
               naive[i].crossover_hz, naive[i].out_centroid_hz, naive[i].peak);

    printf("\n== three groups, three numbers (anchored) ==\n");
    for (int i = 0; i < NSPEED; i++)
        printf("%-7s env_ms=%.1f low_share=%.4f crossover_hz=%.1f\n",
               NAME[i], anch[i].env_ms, anch[i].low_share, anch[i].crossover_hz);

    printf("\n== fastest group: crossover and low share before/after ==\n");
    printf("%-20s %14s %11s\n", "", "crossover_hz", "low_share");
    printf("%-20s %14.1f %11.4f\n", "before (source)", fc0, src.low_share);
    printf("%-20s %14.1f %11.4f\n", "after anchored", anch[2].crossover_hz, anch[2].low_share);
    printf("%-20s %14.1f %11.4f\n", "after naive", naive[2].crossover_hz, naive[2].low_share);

    int collapsed = naive[2].low_share < 0.01;
    printf("\ndiagnosis: fastest low_share=%.4f -> %s\n", naive[2].low_share,
           collapsed ? "collapsed; crossover followed the pitch"
                     : "low band retained; crossover held");

    /* ---- self checks ---- */
    checks_structural(x, ns, fc0, anch);
    if (!external) checks_reference_material(&src, anch, naive);
    else printf("note: material checks skipped for external input\n");

    printf("all_ok=%d\n", failures == 0);
    free(x);
    return failures == 0 ? 0 : 1;
}
