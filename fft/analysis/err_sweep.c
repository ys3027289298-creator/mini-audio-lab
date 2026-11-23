/* err_sweep: characterize rounding-error accumulation of the repo's float
 * FFT (fft/src/fft.c) vs FFT length N, against a double-precision reference
 * of identical algorithmic structure.
 *
 * The reference reuses the SAME recurrent-twiddle algorithm in double, so
 * its own error (~N*eps_double ~ 1e-10 at N=2^20) is negligible against the
 * float path's (~N*eps_float). Input is generated in double, quantized to
 * float, and the reference transforms the SAME quantized values, so only
 * in-transform rounding is measured (input quantization excluded).
 *
 * Usage:
 *   err_sweep sweep            -> CSV: n,max_abs,rms_abs,rms_ref,rel_max
 *   err_sweep bins N file.csv  -> per-bin diffs: bin,re_d,im_d,abs_d,mag_ref
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct { float re, im; } Cf;
typedef struct { double re, im; } Cd;

/* --- algorithm under test: verbatim structure from fft/src/fft.c --- */
static void fft_f(Cf* a, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { Cf t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * 3.14159265358979 / len * (inverse ? 1.0 : -1.0);
        Cf wl = { (float)cos(ang), (float)sin(ang) };
        for (int i = 0; i < n; i += len) {
            Cf w = { 1.0f, 0.0f };
            for (int k = 0; k < len / 2; k++) {
                Cf u = a[i + k];
                Cf v = { a[i + k + len / 2].re * w.re - a[i + k + len / 2].im * w.im,
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

/* --- reference: identical structure in double --- */
static void fft_d(Cd* a, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { Cd t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * 3.14159265358979 / len * (inverse ? 1.0 : -1.0);
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
    if (inverse) { for (int i = 0; i < n; i++) { a[i].re /= n; a[i].im /= n; } }
}

static unsigned long long rng_state = 0x9e3779b97f4a7c15ULL;
static double frand(void) { /* uniform in [-1,1) */
    rng_state = rng_state * 1664525ULL + 1013904223ULL;
    return (double)(rng_state >> 11) / 9007199254740992.0 * 2.0 - 1.0;
}

/* run one length; fills max_abs, rms_abs, rms_ref; optionally dumps per-bin */
static int run_one(int n, double* max_abs, double* rms_abs, double* rms_ref,
                   const char* binfile) {
    Cf* af = malloc(sizeof(Cf) * n);
    Cd* ad = malloc(sizeof(Cd) * n);
    if (!af || !ad) return 1;
    rng_state = 0x9e3779b97f4a7c15ULL; /* same input for every N */
    for (int i = 0; i < n; i++) {
        double re = frand(), im = frand();
        af[i].re = (float)re; af[i].im = (float)im;
        ad[i].re = (double)af[i].re; ad[i].im = (double)af[i].im; /* quantized input */
    }
    fft_f(af, n, 0);
    fft_d(ad, n, 0);

    FILE* bf = binfile ? fopen(binfile, "w") : NULL;
    if (binfile && !bf) { free(af); free(ad); return 1; }
    if (bf) fprintf(bf, "bin,re_diff,im_diff,abs_diff,mag_ref\n");

    double ma = 0.0, sa = 0.0, sr = 0.0;
    for (int i = 0; i < n; i++) {
        double dre = (double)af[i].re - ad[i].re;
        double dim = (double)af[i].im - ad[i].im;
        double d = sqrt(dre * dre + dim * dim);
        double m = sqrt(ad[i].re * ad[i].re + ad[i].im * ad[i].im);
        if (d > ma) ma = d;
        sa += d * d;
        sr += m * m;
        if (bf) fprintf(bf, "%d,%.9g,%.9g,%.9g,%.9g\n", i, dre, dim, d, m);
    }
    if (bf) fclose(bf);
    *max_abs = ma;
    *rms_abs = sqrt(sa / n);
    *rms_ref = sqrt(sr / n);
    free(af); free(ad);
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 2 && argv[1][0] == 'b') { /* bins N file */
        if (argc < 4) { fprintf(stderr, "usage: err_sweep bins N file.csv\n"); return 2; }
        int n = atoi(argv[2]);
        double ma, ra, rr;
        if (run_one(n, &ma, &ra, &rr, argv[3])) return 1;
        fprintf(stderr, "n=%d max_abs=%.6g rms_abs=%.6g rms_ref=%.6g rel_max=%.6g\n",
                n, ma, ra, rr, ma / rr);
        return 0;
    }
    /* sweep mode */
    printf("n,max_abs,rms_abs,rms_ref,rel_max\n");
    for (int p = 4; p <= 22; p++) {
        int n = 1 << p;
        double ma, ra, rr;
        if (run_one(n, &ma, &ra, &rr, NULL)) return 1;
        printf("%d,%.9g,%.9g,%.9g,%.9g\n", n, ma, ra, rr, ma / rr);
        fflush(stdout);
    }
    return 0;
}
