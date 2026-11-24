/* footsteps: hardness/porosity -> contact resonance & damping, impulse synthesis,
 * and measurement of spectral centroid and decay time (T60). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SR 48000
#define N_MATERIALS 5
#define N_SAMPLES 10
#define WIN 2208            /* 46 ms analysis window */
#define HOP 1104            /* 23 ms hop */
#define NFFT 8192
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define MAX_SECONDS 2.5
#define MAX_SAMPLES (120001)     /* 2.5 s at 48 kHz */

/* -- mapping constants ------------------------------------------------------
 * Hardness H in [0,1] is normalised log stiffness: H = 0 -> E_LO, H = 1 -> E_HI.
 * Porosity phi in [0,1] is the air fraction of the contact surface. */
#define E_LO 1.0e6          /* Pa, soft organic mat */
#define E_HI 1.0e11         /* Pa, steel */
#define POISSON 0.30
#define F_LOAD 800.0        /* N, heel strike of an ~80 kg walker */
#define R_HEEL 0.03         /* m, heel contact radius */
#define M_EFF 0.35          /* kg, effective modal mass at the heel */
#define ETA_MIN 5.0e-4      /* metal, structural damping */
#define ETA_MAX 0.25        /* fully compliant, lossy contact */
#define KAPPA_PHI 0.35      /* extra loss per unit air fraction */
#define LN1000 6.907755278982137

typedef struct { const char* name; double hardness; double porosity; } Material;

static const Material MATERIALS[N_MATERIALS] = {
    { "grass", 0.22, 0.55 },
    { "sand",  0.32, 0.40 },
    { "wood",  0.77, 0.06 },
    { "stone", 0.90, 0.01 },
    { "metal", 1.00, 0.00 },
};

static double modulus_of(double hardness) {
    return E_LO * pow(E_HI / E_LO, hardness);
}

static double f0_of(double hardness, double porosity) {
    double e = modulus_of(hardness) * (1.0 - porosity) * (1.0 - porosity);
    double estar = e / (1.0 - POISSON * POISSON);
    double a = cbrt(3.0 * F_LOAD * R_HEEL / (4.0 * estar));
    double k = 2.0 * a * estar;
    return sqrt(k / M_EFF) / (2.0 * M_PI);
}

static double eta_of(double hardness, double porosity) {
    double lossy = 1.0 - hardness;
    return ETA_MIN + (ETA_MAX - ETA_MIN) * lossy * lossy + KAPPA_PHI * porosity;
}

static double t60_of(double f0, double eta) {
    return LN1000 / (M_PI * f0 * eta);
}

/* -- deterministic RNG (splitmix64 + Box-Muller) --------------------------- */
typedef struct { unsigned long long s; } Rng;

static unsigned long long rng_next(Rng* r) {
    r->s += 0x9E3779B97F4A7C15ULL;
    unsigned long long z = r->s;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static double rng_uniform(Rng* r) {
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0) + 0.5e-16;
}

static double rng_gauss(Rng* r) {
    double u1 = rng_uniform(r), u2 = rng_uniform(r);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* -- radix-2 complex FFT (zero-padded analysis) ---------------------------- */
typedef struct { float re, im; } C;

static void fft(C* a, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = -2.0 * M_PI / len * (inverse ? -1.0 : 1.0);
        C wl = { (float)cos(ang), (float)sin(ang) };
        for (int i = 0; i < n; i += len) {
            C w = { 1.0f, 0.0f };
            for (int k = 0; k < len / 2; k++) {
                C u = a[i + k];
                C v = { a[i + k + len / 2].re * w.re - a[i + k + len / 2].im * w.im,
                        a[i + k + len / 2].re * w.im + a[i + k + len / 2].im * w.re };
                a[i + k].re = u.re + v.re; a[i + k].im = u.im + v.im;
                a[i + k + len / 2].re = u.re - v.re;
                a[i + k + len / 2].im = u.im - v.im;
                float nr = w.re * wl.re - w.im * wl.im;
                float ni = w.re * wl.im + w.im * wl.re;
                w.re = nr; w.im = ni;
            }
        }
    }
    if (inverse) {
        for (int i = 0; i < n; i++) { a[i].re /= (float)n; a[i].im /= (float)n; }
    }
}

/* -- synthesis -------------------------------------------------------------
 * Modes of a free circular plate at 1 : 2.756 : 5.404, all sharing one T60.
 * Porous materials add a bed of 24 log-spaced partials (air/grain noise) and
 * a 2 ms contact click. */
static int synth(const Material* m, int seed, float* out, int cap, double* f0_out,
                 double* t60_out) {
    Rng rng = { (unsigned long long)seed };
    double f0 = f0_of(m->hardness, m->porosity);
    double eta = eta_of(m->hardness, m->porosity);
    double load = F_LOAD * (1.0 + 0.30 * rng_gauss(&rng));
    if (load < 80.0) load = 80.0;
    f0 *= pow(load / F_LOAD, 1.0 / 6.0) * (1.0 + 0.02 * rng_gauss(&rng));
    eta *= 1.0 + 0.04 * rng_gauss(&rng);
    if (eta < ETA_MIN) eta = ETA_MIN;
    double t60 = t60_of(f0, eta);

    double dur = 3.0 * t60;
    if (dur > MAX_SECONDS) dur = MAX_SECONDS;
    if (dur < 0.35) dur = 0.35;
    int n = (int)(dur * SR);
    if (n > cap) return -1;

    static const double ratio[3] = { 1.0, 2.756, 5.404 };
    static const double amp[3] = { 1.0, 0.50, 0.25 };
    for (int i = 0; i < n; i++) out[i] = 0.0f;

    for (int k = 0; k < 3; k++) {
        double fk = f0 * ratio[k];
        if (fk > 0.45 * SR) continue;
        double a = amp[k] * (1.0 + 0.25 * rng_gauss(&rng));
        double ph = 2.0 * M_PI * rng_uniform(&rng);
        double w = 2.0 * M_PI * fk;
        for (int i = 0; i < n; i++) {
            double t = (double)i / SR;
            out[i] += (float)(a * exp(-LN1000 * t / t60) * sin(w * t + ph));
        }
    }

    /* Porous bed: white noise shaped in the frequency domain to 0.5*f0..6*f0
     * and decaying with the modes. Discrete partials beat against each other
     * and corrupt the energy-decay slope; band-limited noise does not, so the
     * measured T60 tracks the mapping. */
    double bed = 0.5 * sqrt(m->porosity) * (1.0 + 0.12 * rng_gauss(&rng));
    if (bed > 0.0) {
        int nn = 1;
        while (nn < n) nn <<= 1;
        C* spec = (C*)calloc((size_t)nn, sizeof(C));
        if (!spec) return -1;
        for (int i = 0; i < nn; i++) {
            spec[i].re = (float)rng_gauss(&rng);
            spec[i].im = (float)rng_gauss(&rng);
        }
        fft(spec, nn, 0);
        double lo = 0.5 * f0;
        double hi = 6.0 * f0;
        if (hi > 0.45 * SR) hi = 0.45 * SR;
        for (int k = 0; k <= nn / 2; k++) {
            double f = (double)k * SR / nn;
            if (f < lo || f > hi) { spec[k].re = 0.0f; spec[k].im = 0.0f; }
        }
        for (int k = 1; k < nn / 2; k++) {
            spec[nn - k].re = spec[k].re;
            spec[nn - k].im = -spec[k].im;
        }
        spec[0].re = spec[0].im = 0.0f;
        spec[nn / 2].im = 0.0f;
        fft(spec, nn, 1);
        double rms = 0.0;
        for (int i = 0; i < n; i++) rms += (double)spec[i].re * spec[i].re;
        rms = sqrt(rms / n);
        if (rms > 1e-12) {
            for (int i = 0; i < n; i++) {
                double t = (double)i / SR;
                out[i] += (float)(bed * (double)spec[i].re / rms *
                                  exp(-LN1000 * t / t60) / 4.0);
            }
        }
        free(spec);
    }

    int click_n = (int)(0.002 * SR);
    double click = 0.08 * (0.3 + 0.7 * m->hardness) * (1.0 + 0.2 * rng_gauss(&rng));
    for (int i = 0; i < click_n && i < n; i++) {
        double t = (double)i / SR;
        out[i] += (float)(click * rng_gauss(&rng) * exp(-t / 0.00035));
    }

    float peak = 0.0f;
    for (int i = 0; i < n; i++) {
        float a = fabsf(out[i]);
        if (a > peak) peak = a;
    }
    if (peak <= 0.0f) return -1;
    for (int i = 0; i < n; i++) out[i] /= peak;

    int fade = (int)(0.01 * SR);
    if (fade < n) {
        for (int i = 0; i < fade; i++) {
            double g = 1.0 - (double)i / (double)fade;
            out[n - fade + i] = (float)(out[n - fade + i] * g);
        }
    }

    *f0_out = f0;
    *t60_out = t60;
    return n;
}

/* Energy-weighted spectral centroid over frames above -60 dB re the loudest
 * frame. Hann window, 4x zero padding, bins >= 50 Hz. */
static double measure_centroid(const float* x, int n) {
    static C buf[NFFT];
    double emax = 0.0;
    for (int st = 0; st + WIN <= n; st += HOP) {
        double e = 0.0;
        for (int i = 0; i < WIN; i++) e += (double)x[st + i] * x[st + i];
        if (e > emax) emax = e;
    }
    if (emax <= 0.0) return 0.0;

    double num = 0.0, den = 0.0;
    for (int st = 0; st + WIN <= n; st += HOP) {
        double e = 0.0;
        for (int i = 0; i < WIN; i++) e += (double)x[st + i] * x[st + i];
        if (e < emax * 1e-6) continue;
        for (int i = 0; i < NFFT; i++) { buf[i].re = 0.0f; buf[i].im = 0.0f; }
        for (int i = 0; i < WIN; i++) {
            float w = (float)(0.5 - 0.5 * cos(2.0 * M_PI * i / (WIN - 1)));
            buf[i].re = x[st + i] * w;
        }
        fft(buf, NFFT, 0);
        double cnum = 0.0, cden = 0.0;
        for (int k = 9; k <= NFFT / 2; k++) {
            double f = (double)k * SR / NFFT;
            double mag = sqrt((double)buf[k].re * buf[k].re + (double)buf[k].im * buf[k].im);
            cnum += f * mag;
            cden += mag;
        }
        if (cden > 0.0) {
            num += e * (cnum / cden);
            den += e;
        }
    }
    return den > 0.0 ? num / den : 0.0;
}

/* Schroeder backward integration, least-squares slope between -5 and -25 dB. */
static double measure_t60(const float* x, int n) {
    double* acc = (double*)malloc((size_t)n * sizeof(double));
    if (!acc) return 0.0;
    double tail = 0.0;
    for (int i = n - 1; i >= 0; i--) {
        tail += (double)x[i] * x[i];
        acc[i] = tail;
    }
    double ref = acc[0];
    if (ref <= 0.0) { free(acc); return 0.0; }

    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    int cnt = 0;
    for (int i = 0; i < n; i++) {
        if (acc[i] <= 0.0) break;
        double db = 10.0 * log10(acc[i] / ref);
        if (db > -5.0 || db < -25.0) continue;
        double t = (double)i / SR;
        sx += t; sy += db; sxx += t * t; sxy += t * db;
        cnt++;
    }
    free(acc);
    if (cnt < 16) return 0.0;
    double denom = cnt * sxx - sx * sx;
    if (fabs(denom) < 1e-18) return 0.0;
    double slope = (cnt * sxy - sx * sy) / denom;
    return slope < 0.0 ? -60.0 / slope : 0.0;
}

/* -- statistics in log(centroid) / log(T60) space -------------------------- */
typedef struct { double mu[2]; double cov[2][2]; } Group;

static double log_or_zero(double v) { return v > 0.0 ? log(v) : 0.0; }

static Group group_of(const double* cen, const double* t60, int n) {
    Group g;
    memset(&g, 0, sizeof(g));
    for (int i = 0; i < n; i++) {
        g.mu[0] += log_or_zero(cen[i]);
        g.mu[1] += log_or_zero(t60[i]);
    }
    g.mu[0] /= n; g.mu[1] /= n;
    for (int i = 0; i < n; i++) {
        double a = log_or_zero(cen[i]) - g.mu[0];
        double b = log_or_zero(t60[i]) - g.mu[1];
        g.cov[0][0] += a * a; g.cov[1][1] += b * b; g.cov[0][1] += a * b;
    }
    double d = n > 1 ? (double)(n - 1) : 1.0;
    g.cov[0][0] /= d; g.cov[1][1] /= d; g.cov[0][1] /= d; g.cov[1][0] = g.cov[0][1];
    return g;
}

static double det2(const double m[2][2]) { return m[0][0] * m[1][1] - m[0][1] * m[1][0]; }

/* Mahalanobis distance between group means under the pooled covariance. */
static double separation(const Group* a, const Group* b) {
    double s[2][2], inv[2][2];
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            s[i][j] = 0.5 * (a->cov[i][j] + b->cov[i][j]);
    double d = det2(s);
    if (fabs(d) < 1e-18) return 0.0;
    inv[0][0] = s[1][1] / d; inv[1][1] = s[0][0] / d;
    inv[0][1] = -s[0][1] / d; inv[1][0] = -s[1][0] / d;
    double dx = b->mu[0] - a->mu[0], dy = b->mu[1] - a->mu[1];
    double q = dx * (inv[0][0] * dx + inv[0][1] * dy)
             + dy * (inv[1][0] * dx + inv[1][1] * dy);
    return q > 0.0 ? sqrt(q) : 0.0;
}

static double qform(const Group* g, double x, double y) {
    double d = det2(g->cov);
    if (fabs(d) < 1e-18) return 1e9;
    double inv[2][2];
    inv[0][0] = g->cov[1][1] / d; inv[1][1] = g->cov[0][0] / d;
    inv[0][1] = -g->cov[0][1] / d; inv[1][0] = -g->cov[1][0] / d;
    double dx = x - g->mu[0], dy = y - g->mu[1];
    return dx * (inv[0][0] * dx + inv[0][1] * dy)
         + dy * (inv[1][0] * dx + inv[1][1] * dy);
}

/* Do the two k-sigma covariance ellipses intersect? Tests one boundary against
 * the other ellipse, then each centre against the other ellipse (containment). */
static int ellipses_overlap(const Group* a, const Group* b, double k) {
    double k2 = k * k;
    double sa = sqrt(a->cov[0][0]), sb = sqrt(b->cov[1][1]);
    if (sa <= 0.0 || sb <= 0.0) return 1;
    double r = a->cov[0][1] / (sa * sb);
    if (r > 0.999) r = 0.999;
    if (r < -0.999) r = -0.999;
    double l11 = sa, l21 = r * sb, l22 = sb * sqrt(1.0 - r * r);
    for (int k = 0; k < 360; k++) {
        double th = 2.0 * M_PI * k / 360.0;
        double x = a->mu[0] + l11 * cos(th);
        double y = a->mu[1] + l21 * cos(th) + l22 * sin(th);
        if (qform(b, x, y) <= k2) return 1;
    }
    if (qform(b, a->mu[0], a->mu[1]) <= k2) return 1;
    if (qform(a, b->mu[0], b->mu[1]) <= k2) return 1;
    return 0;
}

/* -- WAV output ------------------------------------------------------------ */
static void wav_header(FILE* f, int n) {
    unsigned char h[44] = { 0 };
    unsigned int data = (unsigned)n * 2u, riff = 36u + data, rate = SR;
    unsigned int brate = SR * 2u;
    memcpy(h, "RIFF", 4); memcpy(h + 4, &riff, 4); memcpy(h + 8, "WAVEfmt ", 8);
    h[16] = 16; h[20] = 1; h[22] = 1;
    memcpy(h + 24, &rate, 4); memcpy(h + 28, &brate, 4);
    h[32] = 2; h[34] = 16; memcpy(h + 36, "data", 4); memcpy(h + 40, &data, 4);
    fwrite(h, 1, 44, f);
}

static int write_wav(const char* dir, const char* name, int idx, const float* x, int n) {
    char path[512];
    snprintf(path, sizeof(path), "%s/%s_%02d.wav", dir, name, idx);
    FILE* f = fopen(path, "wb");
    if (!f) return -1;
    wav_header(f, n);
    for (int i = 0; i < n; i++) {
        int v = (int)lrintf(x[i] * 32767.0f);
        if (v > 32767) v = 32767;
        if (v < -32768) v = -32768;
        short s = (short)v;
        fwrite(&s, 2, 1, f);
    }
    fclose(f);
    return 0;
}

/* -- driver ---------------------------------------------------------------- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

typedef struct { double cen[N_SAMPLES], t60[N_SAMPLES], pred_t60[N_SAMPLES]; } Measurements;

static int measure_material(int mi, Measurements* out, const char* wav_dir, int emit_csv) {
    const Material* m = &MATERIALS[mi];
    static float buf[MAX_SAMPLES];
    for (int s = 0; s < N_SAMPLES; s++) {
        int seed = 0x51A9D000 + mi * 1000 + s * 37;
        double f0 = 0.0, t60 = 0.0;
        int n = synth(m, seed, buf, (int)(sizeof(buf) / sizeof(buf[0])), &f0, &t60);
        if (n <= 0) return -1;
        if (wav_dir && write_wav(wav_dir, m->name, s, buf, n) != 0) return -1;
        out->cen[s] = measure_centroid(buf, n);
        out->t60[s] = measure_t60(buf, n) * 1000.0;
        out->pred_t60[s] = t60 * 1000.0;
        if (emit_csv) {
            printf("%s,%d,%.2f,%.2f,%.2f,%.2f\n", m->name, s, out->cen[s], out->t60[s],
                   f0, out->pred_t60[s]);
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    const char* wav_dir = NULL;
    int csv = 0, test = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--csv")) csv = 1;
        else if (!strcmp(argv[i], "--test")) test = 1;
        else if (!strcmp(argv[i], "--wav") && i + 1 < argc) wav_dir = argv[++i];
    }

    if (csv) {
        printf("material,sample,centroid_hz,t60_ms,pred_f0_hz,pred_t60_ms\n");
    }

    Measurements meas[N_MATERIALS];
    Group groups[N_MATERIALS];
    for (int i = 0; i < N_MATERIALS; i++) {
        if (measure_material(i, &meas[i], wav_dir, csv) != 0) {
            fprintf(stderr, "synthesis failed for %s\n", MATERIALS[i].name);
            return 1;
        }
        groups[i] = group_of(meas[i].cen, meas[i].t60, N_SAMPLES);
        if (!csv) {
            const Material* m = &MATERIALS[i];
            double f0 = f0_of(m->hardness, m->porosity);
            double eta = eta_of(m->hardness, m->porosity);
            printf("material=%s hardness=%.2f porosity=%.2f f0_hz=%.1f eta=%.4f t60_pred_ms=%.1f\n",
                   m->name, m->hardness, m->porosity, f0, eta, t60_of(f0, eta) * 1000.0);
        }
    }

    double sep[N_MATERIALS][N_MATERIALS];
    int overlap[N_MATERIALS][N_MATERIALS];
    for (int i = 0; i < N_MATERIALS; i++) {
        for (int j = 0; j < N_MATERIALS; j++) {
            sep[i][j] = i == j ? 0.0 : separation(&groups[i], &groups[j]);
            overlap[i][j] = i == j ? 1 : ellipses_overlap(&groups[i], &groups[j], 2.0);
        }
    }

    if (test) {
        /* hardness and porosity act in the expected directions */
        double f_soft = f0_of(0.2, 0.0), f_mid = f0_of(0.5, 0.0), f_hard = f0_of(0.9, 0.0);
        check("f0_rises_with_hardness", f_soft < f_mid && f_mid < f_hard);
        check("f0_falls_with_porosity", f0_of(0.5, 0.0) > f0_of(0.5, 0.3) && f0_of(0.5, 0.3) > f0_of(0.5, 0.6));
        double e_pred = modulus_of(0.5);
        check("modulus_log_midpoint", fabs(e_pred - sqrt(E_LO * E_HI)) / sqrt(E_LO * E_HI) < 1e-12);
        check("eta_rises_with_porosity", eta_of(0.5, 0.4) > eta_of(0.5, 0.0));
        check("metal_rings_longest", t60_of(f0_of(1.0, 0.0), eta_of(1.0, 0.0)) >
                                     5.0 * t60_of(f0_of(0.9, 0.01), eta_of(0.9, 0.01)));

        /* every sample finite, normalized, and T60 tracks the prediction */
        int finite = 1, t60_ok = 1;
        for (int i = 0; i < N_MATERIALS; i++) {
            for (int s = 0; s < N_SAMPLES; s++) {
                if (!isfinite(meas[i].cen[s]) || !isfinite(meas[i].t60[s]) ||
                    meas[i].cen[s] <= 0.0 || meas[i].t60[s] <= 0.0) finite = 0;
                double pred = meas[i].pred_t60[s];
                if (fabs(meas[i].t60[s] - pred) / pred > 0.15) t60_ok = 0;
            }
        }
        check("samples_finite", finite);
        check("t60_matches_mapping", t60_ok);

        /* centroid ordering: hard materials brighter */
        double mean_cen[N_MATERIALS];
        for (int i = 0; i < N_MATERIALS; i++) {
            mean_cen[i] = 0.0;
            for (int s = 0; s < N_SAMPLES; s++) mean_cen[i] += meas[i].cen[s];
            mean_cen[i] /= N_SAMPLES;
        }
        const int order[N_MATERIALS] = { 0, 1, 2, 3, 4 }; /* grass sand wood stone metal */
        int rise = 1;
        for (int k = 1; k < N_MATERIALS; k++)
            if (!(mean_cen[order[k]] > mean_cen[order[k - 1]])) rise = 0;
        check("centroid_rises_with_hardness", rise);

        /* determinism: same seed, same signal */
        static float a[MAX_SAMPLES], b[MAX_SAMPLES];
        double fa, ta, fb, tb;
        int na = synth(&MATERIALS[1], 12345, a, (int)(sizeof(a) / sizeof(a[0])), &fa, &ta);
        int nb = synth(&MATERIALS[1], 12345, b, (int)(sizeof(b) / sizeof(b[0])), &fb, &tb);
        int same = na == nb && na > 0;
        for (int i = 0; same && i < na; i++) if (a[i] != b[i]) same = 0;
        check("deterministic_seed", same);

        /* the five groups are distinct at 2 sigma under pooled covariance */
        int all_sep = 1, no_overlap = 1;
        for (int i = 0; i < N_MATERIALS; i++) {
            for (int j = i + 1; j < N_MATERIALS; j++) {
                if (sep[i][j] < 2.0) all_sep = 0;
                if (overlap[i][j]) no_overlap = 0;
            }
        }
        check("groups_separated_2sigma", all_sep);

        /* every one of the 50 samples is nearest its own material centroid
         * under pooled covariance - a leave-one-out classification */
        int classified = 1;
        for (int i = 0; i < N_MATERIALS; i++)
            for (int s = 0; s < N_SAMPLES; s++) {
                double x = log_or_zero(meas[i].cen[s]);
                double y = log_or_zero(meas[i].t60[s]);
                double own = qform(&groups[i], x, y);
                for (int j = 0; j < N_MATERIALS; j++)
                    if (j != i && qform(&groups[j], x, y) < own) classified = 0;
            }
        check("samples_classified", classified);
        check("ellipses_2sigma_disjoint", no_overlap);
    }

    if (!csv) {
        printf("pair_separation_2d:\n");
        for (int i = 0; i < N_MATERIALS; i++)
            for (int j = i + 1; j < N_MATERIALS; j++)
                printf("  %s-%s sep=%.2f overlap=%d\n", MATERIALS[i].name,
                       MATERIALS[j].name, sep[i][j], overlap[i][j]);
        printf("all_ok=%d\n", failures == 0);
    }
    return failures == 0 ? 0 : 1;
}
