
/* blockfft: block-streaming spectrum analyzer (Welch, Hann window, 50% overlap).
 *
 * BLOCKING=1 (default): spectrum is estimated block by block as samples
 *   arrive; the caller pays at most one block of algorithmic latency.
 * BLOCKING=0: blocking is disabled; the whole signal is buffered and one
 *   full-length FFT is taken (pre-blocking, per-sample buffering behavior).
 *
 * The same FFT kernel, windowing and normalization are used in both paths
 * so the only differences are block length B and the number of averages. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef BLOCKING
#define BLOCKING 1
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { float re, im; } C;

/* iterative radix-2 FFT, same kernel as fft/src/fft.c */
static void fft(C* a, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * M_PI / len * (inverse ? 1.0 : -1.0);
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
    if (inverse) { for (int i = 0; i < n; i++) { a[i].re /= n; a[i].im /= n; } }
}

/* ---- block-streaming analyzer: compiled out when blocking is disabled ---- */
static int is_pow2(int n) { return n > 0 && (n & (n - 1)) == 0; }

#if BLOCKING
typedef struct {
    int B;              /* block length (FFT size), power of two */
    int hop;            /* samples between blocks, B/2 */
    int have;           /* valid samples in window */
    int navg;           /* number of averaged periodograms */
    float fs;           /* sample rate (Hz), PSD density normalization */
    float* win;         /* Hann window */
    float* buf;         /* sample buffer (B) */
    double* psd;        /* accumulated one-sided PSD, B/2+1 bins */
    float wsum;         /* sum of squared window weights (power normalization) */
} BlkFFT;

static int blk_init(BlkFFT* st, int B, float fs) {
    if (!is_pow2(B) || B < 4 || fs <= 0.0f) return -1;
    memset(st, 0, sizeof(*st));
    st->B = B;
    st->hop = B / 2;
    st->fs = fs;
    st->win = (float*)malloc((size_t)B * sizeof(float));
    st->buf = (float*)malloc((size_t)B * sizeof(float));
    st->psd = (double*)calloc((size_t)(B / 2 + 1), sizeof(double));
    if (!st->win || !st->buf || !st->psd) { free(st->win); free(st->buf); free(st->psd); return -1; }
    st->wsum = 0.0f;
    for (int i = 0; i < B; i++) {
        st->win[i] = 0.5f - 0.5f * (float)cos(2.0 * M_PI * i / (B - 1));
        st->wsum += st->win[i] * st->win[i];
    }
    return 0;
}

static void blk_free(BlkFFT* st) {
    free(st->win); free(st->buf); free(st->psd);
    memset(st, 0, sizeof(*st));
}

/* transform one full window and fold its periodogram into the average.
 * one-sided density convention P = |Xw|^2 / (fs * sum w^2), units 1/Hz;
 * interior bins carry the one-sided x2 (DC and Nyquist do not), matching
 * scipy.signal.welch(return_onesided=True). */
static void blk_emit(BlkFFT* st) {
    C* z = (C*)malloc((size_t)st->B * sizeof(C));
    for (int i = 0; i < st->B; i++) { z[i].re = st->buf[i] * st->win[i]; z[i].im = 0.0f; }
    fft(z, st->B, 0);
    double scale = 1.0 / ((double)st->fs * (double)st->wsum);
    int half = st->B / 2;
    for (int k = 0; k <= half; k++) {
        double p = (double)z[k].re * z[k].re + (double)z[k].im * z[k].im;
        st->psd[k] += p * scale * (k == 0 || k == half ? 1.0 : 2.0);
    }
    st->navg++;
    free(z);
}

/* push n samples through the streaming analyzer */
static void blk_push(BlkFFT* st, const float* x, int n) {
    for (int i = 0; i < n; i++) {
        st->buf[st->have++] = x[i];
        if (st->have == st->B) {
            blk_emit(st);
            memmove(st->buf, st->buf + st->hop, (size_t)(st->B - st->hop) * sizeof(float));
            st->have = st->B - st->hop;
        }
    }
}

/* close the estimate: the residual partial window is discarded (the last
 * full window already covered those samples; zero-padding it would only
 * dilute the average). only the per-block normalization remains. */
static void blk_flush(BlkFFT* st) {
    st->have = 0;
    if (st->navg > 0)
        for (int k = 0; k <= st->B / 2; k++) st->psd[k] /= st->navg;
}
#endif /* BLOCKING */

/* reference path: blocking disabled — buffer everything, one full-length FFT.
 * identical Hann window and density normalization as the block path. */
static int full_spectrum(const float* x, int N, float fs, double* psd_out) {
    if (!is_pow2(N) || fs <= 0.0f) return -1;
    float* win = (float*)malloc((size_t)N * sizeof(float));
    C* z = (C*)malloc((size_t)N * sizeof(C));
    if (!win || !z) { free(win); free(z); return -1; }
    float wsum = 0.0f;
    for (int i = 0; i < N; i++) {
        win[i] = 0.5f - 0.5f * (float)cos(2.0 * M_PI * i / (N - 1));
        wsum += win[i] * win[i];
    }
    for (int i = 0; i < N; i++) { z[i].re = x[i] * win[i]; z[i].im = 0.0f; }
    fft(z, N, 0);
    double scale = 1.0 / ((double)fs * (double)wsum);
    for (int k = 0; k <= N / 2; k++) {
        double p = (double)z[k].re * z[k].re + (double)z[k].im * z[k].im;
        psd_out[k] = p * scale * (k == 0 || k == N / 2 ? 1.0 : 2.0);
    }
    free(win); free(z);
    return 0;
}
