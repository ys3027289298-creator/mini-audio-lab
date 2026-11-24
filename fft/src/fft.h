
/* fft: iterative radix-2 complex FFT. */
#ifndef MINI_FFT_H
#define MINI_FFT_H

#include <math.h>

typedef struct { float re, im; } FftC;

static void fft(FftC* a, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { FftC t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * 3.14159265358979 / len * (inverse ? 1.0 : -1.0);
        FftC wl = { (float)cos(ang), (float)sin(ang) };
        for (int i = 0; i < n; i += len) {
            FftC w = { 1.0f, 0.0f };
            for (int k = 0; k < len / 2; k++) {
                FftC u = a[i + k];
                FftC v = { a[i + k + len / 2].re * w.re - a[i + k + len / 2].im * w.im,
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

#endif
