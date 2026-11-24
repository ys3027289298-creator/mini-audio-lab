/*
 * stream.h - unified output-length derivation, single source of truth.
 *
 * Every stage in the streaming pipeline derives its output length from the
 * same whole-stream formulas instead of rounding each block independently.
 *
 * ratio = sr_in / sr_out (identical convention to resample_linear():
 * output index i reads source position i * ratio).
 *
 * For a stream of N input samples split into blocks of B samples:
 *   total resampler output      M = floor(N / ratio)
 *   output of block b           floor((b+1)B / ratio) - floor(bB / ratio)
 *   stereo pairs over M frames  floor(M / 2), odd frame carried to next block
 *
 * Pure header-only, no dependencies beyond stdint.h.
 */
#ifndef MINI_AUDIO_LAB_STREAM_H
#define MINI_AUDIO_LAB_STREAM_H

#include <stdint.h>

#if defined(__GNUC__)
#define MA_MAYBE_UNUSED __attribute__((unused))
#else
#define MA_MAYBE_UNUSED
#endif

/* Whole-stream resampler output length for n_in samples at ratio r. */
static inline int64_t stream_total_out(int64_t n_in, double ratio) {
    return (int64_t)((double)n_in / ratio); /* floor of positive value */
}

/* Output length of block `block` (0-based), blocks of block_in input samples. */
static inline int64_t stream_block_out(int64_t block, int64_t block_in,
                                       double ratio) {
    int64_t hi = (int64_t)((double)((block + 1) * block_in) / ratio);
    int64_t lo = (int64_t)((double)(block * block_in) / ratio);
    return hi - lo;
}

/* Total stereo pairs produced from n mono frames (one frame may carry over). */
static inline int64_t stream_total_pairs(int64_t n) { return n / 2; }

/*
 * Stereo pairs produced by one block, given the carry flag from the previous
 * block (1 = one unpaired frame pending). The new carry is written back.
 */
static inline int64_t stream_block_pairs(int64_t n_in_block, int *carry) {
    int64_t frames = (int64_t)(*carry) + n_in_block;
    *carry = (int)(frames & 1);
    return frames / 2;
}

#endif /* MINI_AUDIO_LAB_STREAM_H */
