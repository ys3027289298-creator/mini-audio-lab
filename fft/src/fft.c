
/* fft: iterative radix-2 complex FFT with magnitude and Parseval check. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define NFFT 64

typedef struct { float re, im; } C;

static void fft(C* a, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * 3.14159265358979 / len * (inverse ? 1.0 : -1.0);
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

static float energy_time(const C* a, int n) {
    float s = 0.0f;
    for (int i = 0; i < n; i++) s += a[i].re * a[i].re + a[i].im * a[i].im;
    return s;
}

static float energy_freq(const C* a, int n) {
    float s = 0.0f;
    for (int i = 0; i < n; i++) s += a[i].re * a[i].re + a[i].im * a[i].im;
    return s / (float)n;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    C a[NFFT];
    memset(a, 0, sizeof(a));

    /* DC signal */
    for (int i = 0; i < NFFT; i++) a[i].re = 1.0f;
    float e0 = energy_time(a, NFFT);
    fft(a, NFFT, 0);
    check("dc_bin0", fabsf(a[0].re - (float)NFFT) < 1e-3f);
    check("dc_others_zero", fabsf(a[5].re) < 1e-3f);
    check("parseval", fabsf(energy_freq(a, NFFT) - e0) / e0 < 1e-3f);

    /* single tone at bin 4 */
    for (int i = 0; i < NFFT; i++) {
        double ph = 2.0 * 3.14159265358979 * 4.0 * i / NFFT;
        a[i].re = (float)cos(ph); a[i].im = 0.0f;
    }
    fft(a, NFFT, 0);
    float mag4 = sqrtf(a[4].re * a[4].re + a[4].im * a[4].im);
    float mag3 = sqrtf(a[3].re * a[3].re + a[3].im * a[3].im);
    check("tone_peak_bin", mag4 > 10.0f * mag3);
    check("tone_amplitude", fabsf(mag4 / (NFFT / 2) - 1.0f) < 1e-2f);

    /* round trip */
    C b[NFFT];
    for (int i = 0; i < NFFT; i++) { b[i].re = (float)(i % 5); b[i].im = 0.0f; }
    C orig[NFFT];
    memcpy(orig, b, sizeof(orig));
    fft(b, NFFT, 0);
    fft(b, NFFT, 1);
    float maxd = 0.0f;
    for (int i = 0; i < NFFT; i++) {
        float d = fabsf(b[i].re - orig[i].re);
        if (d > maxd) maxd = d;
    }
    check("roundtrip", maxd < 1e-3f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
