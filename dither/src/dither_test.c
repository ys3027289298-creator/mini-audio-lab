/* dither: self-test for the TPDF dither / 16-bit quantizer. */
#include "dither.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define N 48000
#define QSCALE 32768.0f

int main(void) {
    static float in[N], out_d[N], out_u[N];
    DitherReport rd, ru;
    DitherRng rng;

    /* -60 dBFS 1 kHz tone: the textbook dither case. */
    for (int i = 0; i < N; i++)
        in[i] = 0.001f * sinf(2.0f * 3.14159265f * 1000.0f * (float)i / 48000.0f);

    dither_rng_seed(&rng, 1234, 54);
    check("tpdf_ok", dither_quantize16(in, N, out_d, 1, &rng, &rd) == N);
    check("plain_ok", dither_quantize16(in, N, out_u, 0, NULL, &ru) == N);
    check("no_clip", rd.clipped == 0 && ru.clipped == 0);

    /* TPDF error variance should land near 1/4 LSB^2. */
    check("tpdf_err_rms", rd.err_rms > 0.45f && rd.err_rms < 0.55f);
    check("tpdf_err_dc", fabs(rd.err_dc) < 0.01);

    /* A -60 dBFS sine is only 32 LSB tall; dither must move many codes. */
    check("dither_moves_codes", rd.changed > N / 4);

    /* Determinism: same seed => identical bytes; different seed => different. */
    static float out_again[N];
    dither_rng_seed(&rng, 1234, 54);
    dither_quantize16(in, N, out_again, 1, &rng, NULL);
    check("same_seed_identical", memcmp(out_d, out_again, sizeof(float) * N) == 0);

    dither_rng_seed(&rng, 99, 54);
    dither_quantize16(in, N, out_again, 1, &rng, NULL);
    check("diff_seed_differs", memcmp(out_d, out_again, sizeof(float) * N) != 0);

    /* Dither off => plain rounding, bit for bit, with no RNG involved. */
    static float out_plain[N];
    check("dither_off_ok", dither_quantize16(in, N, out_plain, 0, NULL, NULL) == N);
    check("dither_off_plain", memcmp(out_plain, out_u, sizeof(float) * N) == 0);

    /* Silence: only dither makes an idle LSB toggle. */
    memset(in, 0, sizeof(float) * N);
    dither_rng_seed(&rng, 7, 54);
    dither_quantize16(in, N, out_d, 1, &rng, &rd);
    dither_quantize16(in, N, out_u, 0, NULL, &ru);
    int nz = 0;
    for (int i = 0; i < N; i++) nz += out_d[i] != 0.0f;
    check("silence_dither_idles", nz > 0);
    check("silence_plain_zero", ru.err_rms == 0.0f && rd.err_rms > 0.45f);

    /* 0.3 LSB DC is below the rounding threshold: plain rounding collapses it
       to silence, TPDF dither keeps the mean as a duty cycle of codes. */
    for (int i = 0; i < N; i++) in[i] = 0.3f / QSCALE;
    dither_rng_seed(&rng, 3, 54);
    dither_quantize16(in, N, out_d, 1, &rng, &rd);
    dither_quantize16(in, N, out_u, 0, NULL, &ru);
    float mean_d = 0.0f;
    for (int i = 0; i < N; i++) mean_d += out_d[i];
    mean_d = mean_d / (float)N * QSCALE;
    check("sub_lsb_plain_stuck", out_u[0] == 0.0f && out_u[N - 1] == 0.0f);
    check("sub_lsb_dither_mean", fabsf(mean_d - 0.3f) < 0.02f);

    /* Clipping is reported, not hidden. TPDF adds up to +1 LSB, so +1.0 clips
       on the positive rail; 32766/32768 is the dither-safe peak. */
    for (int i = 0; i < N; i++) in[i] = (i % 2 ? 1.0f : -1.0f);
    dither_rng_seed(&rng, 5, 54);
    dither_quantize16(in, N, out_d, 1, &rng, &rd);
    check("fullscale_clip_pos_rail", rd.clipped > 0);

    for (int i = 0; i < N; i++) in[i] = (i % 2 ? 1.0f : -1.0f) * (32766.0f / QSCALE);
    dither_rng_seed(&rng, 5, 54);
    dither_quantize16(in, N, out_d, 1, &rng, &rd);
    check("near_full_no_clip", rd.clipped == 0);

    for (int i = 0; i < N; i++) in[i] = 1.5f;
    dither_quantize16(in, N, out_d, 1, &rng, &rd);
    check("over_clip_reported", rd.clipped == N);

    check("bad_n", dither_quantize16(in, 0, out_d, 1, &rng, &rd) == -1);
    check("bad_mode", dither_quantize16(in, N, out_d, 2, &rng, &rd) == -1);
    check("missing_rng", dither_quantize16(in, N, out_d, 1, NULL, &rd) == -1);
    in[0] = NAN;
    check("bad_nan", dither_quantize16(in, N, out_d, 1, &rng, &rd) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
