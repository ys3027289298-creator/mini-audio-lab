/* dither_core.c: TPDF dither + 16-bit quantization, shared by the self-test and the CLI. */
#include "dither.h"
#include <math.h>
#include <stdlib.h>

/*
 * Quantizer convention: x maps to v = x * 32768; the integer is
 * floor(v + dither + 0.5) clamped to [-32768, 32767]; the decoded float is
 * q / 32768. Using 32768 (not 32767) makes one quantization step exactly
 * 1.0 LSB in the decoded domain and keeps the error mean-zero at full scale.
 */
#define QSCALE 32768.0f

void dither_rng_seed(DitherRng* r, uint64_t seed, uint64_t stream) {
    r->state = 0u;
    r->inc = (stream << 1u) | 1u;
    r->state = r->state + seed;
    (void)dither_rng_next(r);
}

uint32_t dither_rng_next(DitherRng* r) {
    uint64_t old = r->state;
    r->state = old * 6364136223846793005ull + r->inc;
    uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = (uint32_t)(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((32u - rot) & 31u));
}

/* 24-bit uniform: values are k / 2^24 - 0.5, so exactly 2^24 evenly spaced. */
float dither_rng_uni(DitherRng* r) {
    return (float)(dither_rng_next(r) >> 8) * (1.0f / 16777216.0f) - 0.5f;
}

float dither_rng_tpdf(DitherRng* r) {
    return dither_rng_uni(r) + dither_rng_uni(r);
}

int dither_quantize16_code(const float* in, int n, int16_t* codes, int dither,
                           DitherRng* rng, DitherReport* rep) {
    if (in == NULL || codes == NULL || n <= 0) return -1;
    if (dither != 0 && dither != 1) return -1;
    if (dither == 1 && rng == NULL) return -1;

    float peak = 0.0f;
    double err_sq = 0.0, err_sum = 0.0;
    int64_t clipped = 0, changed = 0;

    for (int i = 0; i < n; i++) {
        float x = in[i];
        if (!isfinite(x)) return -1;
        if (fabsf(x) > peak) peak = fabsf(x);

        float v = x * QSCALE;
        float d = dither ? dither_rng_tpdf(rng) : 0.0f;
        float r = floorf(v + d + 0.5f);
        int plain = (int)floorf(v + 0.5f);   /* dither-free code, for the changed counter */
        int q;
        if (r > 32767.0f) { q = 32767; clipped++; }
        else if (r < -32768.0f) { q = -32768; clipped++; }
        else q = (int)r;
        if (q != plain) changed++;

        codes[i] = (int16_t)q;
        double e = (double)((float)q / QSCALE) - (double)x;
        err_sq += e * e;
        err_sum += e;
    }

    if (rep) {
        rep->peak = peak;
        rep->err_rms = (float)sqrt(err_sq / (double)n) * QSCALE;
        rep->err_dc = err_sum / (double)n * QSCALE;
        rep->clipped = clipped;
        rep->changed = changed;
    }
    return n;
}

int dither_quantize16(const float* in, int n, float* out, int dither,
                      DitherRng* rng, DitherReport* rep) {
    /* float-domain convenience wrapper, used by the self-test; not in the
       export path, which writes the codes directly. */
    if (in == NULL || out == NULL || n <= 0) return -1;
    int16_t* codes = (int16_t*)malloc(sizeof(int16_t) * (size_t)n);
    if (!codes) return -1;
    int rc = dither_quantize16_code(in, n, codes, dither, rng, rep);
    if (rc == n)
        for (int i = 0; i < n; i++) out[i] = (float)codes[i] / QSCALE;
    free(codes);
    return rc;
}
