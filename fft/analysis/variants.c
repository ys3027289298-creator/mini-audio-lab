/* variants: isolate the rounding-error source in the repo's float FFT.
 *   mode 0: baseline        - float data, float recurrent twiddle (repo code)
 *   mode 1: direct twiddle  - float data, w recomputed per k as float cos/sin
 *   mode 2: double twiddle  - float data, recurrence kept in double, applied as float
 *   mode 3: double data     - everything double (reference sanity floor)
 * Prints: mode,n,rms_rel   (rms_abs / rms_ref vs double reference)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct { float re, im; } Cf;
typedef struct { double re, im; } Cd;

static void fft_var(Cf* a, int n, int mode) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { Cf t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * 3.14159265358979 / len * -1.0;
        Cf wl = { (float)cos(ang), (float)sin(ang) };
        double wlr = cos(ang), wli = sin(ang);
        for (int i = 0; i < n; i += len) {
            Cf w = { 1.0f, 0.0f };
            double wr = 1.0, wi = 0.0;
            for (int k = 0; k < len / 2; k++) {
                float wre, wim;
                if (mode == 1) {          /* direct per-k twiddle */
                    wre = (float)cos(ang * k); wim = (float)sin(ang * k);
                } else if (mode == 2) {   /* double recurrence, float apply */
                    wre = (float)wr; wim = (float)wi;
                } else {                  /* baseline float recurrence */
                    wre = w.re; wim = w.im;
                }
                Cf u = a[i + k];
                Cf v = { a[i + k + len / 2].re * wre - a[i + k + len / 2].im * wim,
                         a[i + k + len / 2].re * wim + a[i + k + len / 2].im * wre };
                a[i + k].re = u.re + v.re; a[i + k].im = u.im + v.im;
                a[i + k + len / 2].re = u.re - v.re; a[i + k + len / 2].im = u.im - v.im;
                float nr = w.re * wl.re - w.im * wl.im;
                float ni = w.re * wl.im + w.im * wl.re;
                w.re = nr; w.im = ni;
                double dnr = wr * wlr - wi * wli;
                double dni = wr * wli + wi * wlr;
                wr = dnr; wi = dni;
            }
        }
    }
}

static void fft_d(Cd* a, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { Cd t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * 3.14159265358979 / len * -1.0;
        Cd wl = { cos(ang), sin(ang) };
        for (int i = 0; i < n; i += len) {
            Cd w = { 1.0, 0.0 };
            for (int k = 0; k < len / 2; k++) {
                Cd u = a[i + k];
                Cd v = { a[i + k + len / 2].re * w.re - a[i + k + len / 2].im * w.im,
                         a[i + k + len / 2].re * w.im + a[i + k + len / 2].im * w.re };
                a[i + k].re = u.re + v.re; a[i + k].im = u.im + v.im;
                a[i + k + len / 2].re = u.re - v.re; a[i + k + len / 2].im = u.im - v.im;
                double nr = w.re * wl.re - w.im * wl.im;
                double ni = w.re * wl.im + w.im * wl.re;
                w.re = nr; w.im = ni;
            }
        }
    }
}

static unsigned long long rng_state;
static double frand(void) {
    rng_state = rng_state * 1664525ULL + 1013904223ULL;
    return (double)(rng_state >> 11) / 9007199254740992.0 * 2.0 - 1.0;
}

int main(void) {
    printf("mode,n,rms_rel\n");
    for (int mode = 0; mode <= 2; mode++) {
        for (int p = 6; p <= 20; p += 2) {
            int n = 1 << p;
            Cf* af = malloc(sizeof(Cf) * n);
            Cd* ad = malloc(sizeof(Cd) * n);
            rng_state = 0x9e3779b97f4a7c15ULL;
            for (int i = 0; i < n; i++) {
                double re = frand(), im = frand();
                af[i].re = (float)re; af[i].im = (float)im;
                ad[i].re = (double)af[i].re; ad[i].im = (double)af[i].im;
            }
            fft_var(af, n, mode);
            fft_d(ad, n);
            double sa = 0.0, sr = 0.0;
            for (int i = 0; i < n; i++) {
                double dre = (double)af[i].re - ad[i].re;
                double dim = (double)af[i].im - ad[i].im;
                sa += dre * dre + dim * dim;
                sr += ad[i].re * ad[i].re + ad[i].im * ad[i].im;
            }
            printf("%d,%d,%.9g\n", mode, n, sqrt(sa / sr));
            free(af); free(ad);
        }
    }
    return 0;
}
