
/* tfres: time-frequency resolution.
 * 1) Gabor/Heisenberg lower bound dt*dw >= 1/2  (dt*df >= 1/(4*pi)).
 * 2) window-type deviation from the bound (RMS widths, measured).
 * 3) three candidate window lengths on a two-tone signal; select one.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

typedef struct { double re, im; } C;

static void fft(C* a, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * PI / len;
        C wl = { cos(ang), -sin(ang) };
        for (int i = 0; i < n; i += len) {
            C w = { 1.0, 0.0 };
            for (int k = 0; k < len / 2; k++) {
                C u = a[i + k];
                C v = { a[i+k+len/2].re*w.re - a[i+k+len/2].im*w.im,
                        a[i+k+len/2].re*w.im + a[i+k+len/2].im*w.re };
                a[i+k]         = (C){ u.re + v.re, u.im + v.im };
                a[i+k+len/2]   = (C){ u.re - v.re, u.im - v.im };
                w = (C){ w.re*wl.re - w.im*wl.im, w.re*wl.im + w.im*wl.re };
            }
        }
    }
}

/* windows on x in [-0.5, 0.5] */
enum { W_GAUSS, W_HANN, W_HAMMING, W_RECT };
static double win(int type, double x) {
    switch (type) {
    case W_GAUSS:   { double s = x / 0.12; return exp(-0.5 * s * s); }
    case W_HANN:    return 0.5 + 0.5 * cos(2.0 * PI * x);
    case W_HAMMING: return 0.54 + 0.46 * cos(2.0 * PI * x);
    case W_RECT:    return 1.0;
    }
    return 0.0;
}

/* RMS widths of a window: L samples across support, K-point DTFT grid.
 * dt in samples, dw in rad/sample; product dt*dw is grid-independent
 * for windows whose RMS bandwidth converges. */
static double tf_product(int type, int L, int K, double* dt_out, double* dw_out) {
    C* W = calloc((size_t)K, sizeof(C));
    for (int n = 0; n < L; n++)
        W[n].re = win(type, (n + 0.5) / L - 0.5);
    fft(W, K);

    double e = 0.0, mt = 0.0;
    for (int n = 0; n < L; n++) {
        double x = (n + 0.5) / L - 0.5, w = win(type, x);
        e += w * w; mt += x * w * w;
    }
    double tbar = mt / e, vt = 0.0;
    for (int n = 0; n < L; n++) {
        double x = (n + 0.5) / L - 0.5 - tbar, w = win(type, x);
        vt += x * x * w * w;
    }
    double dt = sqrt(vt / e) * L;          /* samples */

    double ef = 0.0, vf = 0.0;
    for (int k = 0; k < K; k++) {
        double om = 2.0 * PI * k / K;
        if (om > PI) om -= 2.0 * PI;
        double m2 = W[k].re * W[k].re + W[k].im * W[k].im;
        ef += m2; vf += om * om * m2;
    }
    free(W);
    *dt_out = dt;
    *dw_out = sqrt(vf / ef);
    return dt * (*dw_out);
}

/* hann-windowed two-tone segment of length N, zero-padded to K: magnitude dB */
static double* tone_spectrum(int N, int K, double fs, double f1, double f2) {
    C* a = calloc((size_t)K, sizeof(C));
    for (int n = 0; n < N; n++) {
        double s = cos(2.0 * PI * f1 * n / fs) + cos(2.0 * PI * f2 * n / fs);
        a[n].re = s * win(W_HANN, (n + 0.5) / N - 0.5);
    }
    fft(a, K);
    double* db = malloc((size_t)(K / 2) * sizeof(double));
    for (int k = 0; k < K / 2; k++)
        db[k] = 20.0 * log10(hypot(a[k].re, a[k].im) + 1e-12);
    free(a);
    return db;
}

/* depth of the valley between the two tone peaks, dB */
static double tone_dip(const double* db, int K, double fs, double f1, double f2) {
    int span = (int)(40.0 * K / fs);
    int k1 = (int)(f1 * K / fs), k2 = (int)(f2 * K / fs);
    double p1 = -1e9, p2 = -1e9;
    int kp1 = k1, kp2 = k2;
    for (int k = k1 - span; k <= k1 + span; k++) if (db[k] > p1) { p1 = db[k]; kp1 = k; }
    for (int k = k2 - span; k <= k2 + span; k++) if (db[k] > p2) { p2 = db[k]; kp2 = k; }
    double val = 1e9;
    for (int k = kp1; k <= kp2; k++) if (db[k] < val) val = db[k];
    double pk = p1 < p2 ? p1 : p2;
    return pk - val;
}

static void plot(const double* db, int K, double fs, double flo, double fhi,
                 double step, const char* title) {
    double mx = -1e9;
    for (double f = flo; f <= fhi + 1e-9; f += step) {
        int k = (int)(f * K / fs + 0.5);
        if (db[k] > mx) mx = db[k];
    }
    printf("%s\n", title);
    for (double f = flo; f <= fhi + 1e-9; f += step) {
        int k = (int)(f * K / fs + 0.5);
        double v = db[k] - mx;
        if (v < -48.0) v = -48.0;
        int bars = (int)((v + 48.0) / 48.0 * 40.0 + 0.5);
        printf("%5.0f |", f);
        for (int i = 0; i < bars; i++) putchar('#');
        putchar('\n');
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    printf("== bound ==\n");
    printf("dt*dw >= 1/2 (rad)   <=>   dt*df >= 1/(4*pi) = %.5f s*Hz\n\n", 1.0 / (4.0 * PI));

    /* ---- window-type deviation ---- */
    const int L = 4096, K = 64 * L;
    const char* names[4] = { "gauss", "hann", "hamming", "rect" };
    double prods[4], dts[4], dws[4];
    printf("== window deviation (L=%d, DTFT grid=64L) ==\n", L);
    printf("%-9s %8s %10s %9s %8s\n", "window", "dt/T", "dw*T", "dt*dw", "x bound");
    for (int t = 0; t < 4; t++) {
        prods[t] = tf_product(t, L, K, &dts[t], &dws[t]);
        printf("%-9s %8.4f %10.4f %9.4f %8.4f\n",
               names[t], dts[t] / L, dws[t] * L, prods[t], prods[t] / 0.5);
    }
    /* RMS bandwidth of a window with an edge jump diverges as sqrt(L):
     * dt ~ L, dw^2 ~ |jump|^2/L, so dt*dw ~ |jump|*sqrt(L). */
    double r1, r2, r3, d_, w_;
    r1 = tf_product(W_RECT, 1024,  16 * 1024,  &d_, &w_);
    r2 = tf_product(W_RECT, 4096,  16 * 4096,  &d_, &w_);
    r3 = tf_product(W_RECT, 16384, 16 * 16384, &d_, &w_);
    printf("rect dt*dw at L=1024/4096/16384: %.3f / %.3f / %.3f"
           "  (~sqrt(L): RMS bandwidth diverges, edge jump=1)\n", r1, r2, r3);
    printf("hamming diverges the same way, scaled by its edge jump 0.08\n\n");

    /* ---- three candidate window lengths ---- */
    const double fs = 48000.0, f1 = 2000.0, f2 = 2120.0;
    const double spacing = f2 - f1;
    const double enbw = 1.5;                 /* hann, bins */
    const double dtT = dts[W_HANN] / L;      /* hann RMS duration, frac of T */
    const double dwT = dws[W_HANN] * L;      /* hann RMS bandwidth, rad per T */
    const int Ns[3] = { 256, 1024, 4096 };
    double dips[3];
    int K2[3];

    printf("== window length candidates (fs=48000, hann, tones %.0f/%.0f Hz) ==\n", f1, f2);
    printf("%-6s %11s %11s %12s %12s %9s\n",
           "N", "dt_rms ms", "df_rms Hz", "dt*df s.Hz", "ENBW df Hz", "dip dB");
    for (int i = 0; i < 3; i++) {
        int N = Ns[i], k = 1;
        while (k < 16 * N) k <<= 1;
        K2[i] = k;
        double* db = tone_spectrum(N, k, fs, f1, f2);
        dips[i] = tone_dip(db, k, fs, f1, f2);
        double dtms = dtT * N / fs * 1000.0;
        double df = dwT / N * fs / (2.0 * PI);
        printf("%-6d %11.2f %11.2f %12.5f %12.1f %9.2f\n",
               N, dtms, df, dtms / 1000.0 * df, enbw * fs / N, dips[i]);
        free(db);
    }
    putchar('\n');

    /* ---- three-point spectrum comparison ---- */
    for (int i = 0; i < 3; i++) {
        double* db = tone_spectrum(Ns[i], K2[i], fs, f1, f2);
        char title[80];
        snprintf(title, sizeof(title), "-- N=%d (dip %.1f dB) --", Ns[i], dips[i]);
        plot(db, K2[i], fs, 1800.0, 2300.0, 25.0, title);
        free(db);
    }

    /* ---- selection: smallest N resolving the tone spacing ---- */
    int sel = -1;
    for (int i = 0; i < 3; i++)
        if (enbw * fs / Ns[i] <= 0.6 * spacing) { sel = i; break; }
    printf("== selection ==\n");
    if (sel >= 0)
        printf("selected N=%d: ENBW df %.1f Hz <= 0.6*spacing (%.0f Hz); "
               "smaller N cannot separate the tones, larger N only adds "
               "temporal smear (dt_rms %.2f ms -> %.2f ms at N=%d)\n",
               Ns[sel], enbw * fs / Ns[sel], 0.6 * spacing,
               dtT * Ns[sel] / fs * 1000.0, dtT * Ns[2] / fs * 1000.0, Ns[2]);
    else
        printf("no candidate resolves the spacing\n");

    check("bound_gauss", fabs(prods[W_GAUSS] - 0.5) < 0.01);
    check("hann_near_bound", prods[W_HANN] / 0.5 > 1.0 && prods[W_HANN] / 0.5 < 1.06);
    check("rect_diverges", r3 > 3.0 * r1);
    check("n256_merged", dips[0] < 1.0);
    check("n1024_resolved", dips[1] > 4.0);
    check("n4096_deeper", dips[2] > dips[1]);
    check("selected_1024", sel == 1);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
