
/* fft: self-test for fft.h — magnitude and Parseval check. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fft.h"

#define NFFT 64

static float energy_time(const FftC* a, int n) {
    float s = 0.0f;
    for (int i = 0; i < n; i++) s += a[i].re * a[i].re + a[i].im * a[i].im;
    return s;
}

static float energy_freq(const FftC* a, int n) {
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
    FftC a[NFFT];
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
    FftC b[NFFT];
    for (int i = 0; i < NFFT; i++) { b[i].re = (float)(i % 5); b[i].im = 0.0f; }
    FftC orig[NFFT];
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
