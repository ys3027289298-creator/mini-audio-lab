/* dither.h: TPDF dither and 16-bit quantization. */
#ifndef DITHER_H
#define DITHER_H

#include <stdint.h>
#include <stddef.h>

/* PCG32 (XSH RR 64/32). Fixed seed => fixed dither sequence => reproducible bytes. */
typedef struct { uint64_t state, inc; } DitherRng;

void dither_rng_seed(DitherRng* r, uint64_t seed, uint64_t stream);
uint32_t dither_rng_next(DitherRng* r);
float dither_rng_uni(DitherRng* r);    /* uniform on [-0.5, 0.5) */
float dither_rng_tpdf(DitherRng* r);   /* TPDF on (-1, 1) LSB */

typedef struct {
    float peak;       /* max |input|, full scale = 1.0 */
    float err_rms;    /* rms of (decoded - input), in LSB */
    double err_dc;    /* mean of (decoded - input), in LSB */
    int64_t clipped;  /* samples pinned at -32768 or 32767 */
    int64_t changed;  /* samples whose code differs from plain rounding */
} DitherReport;

/*
 * Quantize n interleaved samples to 16-bit codes (int16_t, the exact values
 * written to a PCM file). dither: 1 = TPDF, 0 = plain rounding (rng may be
 * NULL then). Returns n, or -1 on bad input (n <= 0, bad mode, non-finite
 * sample). The report's error figures are in LSB.
 */
int dither_quantize16_code(const float* in, int n, int16_t* codes, int dither,
                           DitherRng* rng, DitherReport* rep);

/*
 * Same, but writes the decoded float (q / 32768) instead of the codes, so
 * callers can compare output against input directly in the float domain.
 */
int dither_quantize16(const float* in, int n, float* out, int dither,
                      DitherRng* rng, DitherReport* rep);

#endif
