/* fft_core.c — iterative radix-2 complex FFT, shared by fft.c and spectral.c. */
#include <math.h>

#include "fft_internal.h"

void fft_forward(C *a, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = -2.0 * 3.14159265358979 / len;
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

void fft_inverse(C *a, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * 3.14159265358979 / len;
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
    for (int i = 0; i < n; i++) { a[i].re /= n; a[i].im /= n; }
}
