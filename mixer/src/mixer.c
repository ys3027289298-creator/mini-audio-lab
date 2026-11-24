
/* mixer: float mixing with ramp gain and clipping detection. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "stream.h"

typedef struct { float peak; int clipped; int samples; } MixStat;

static int mix_ramp(const float* a, const float* b, float gain_a, float gain_b,
                    float* out, int n, int ramp_frames, MixStat* st) {
    if (n <= 0) return -1;
    if (ramp_frames < 0 || ramp_frames > n) return -1;
    if (gain_a < 0.0f || gain_b < 0.0f) return -1;
    float peak = 0.0f; int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = ramp_frames > 0 ? (float)i / (float)ramp_frames : 1.0f;
        if (t > 1.0f) t = 1.0f;
        float ga = gain_a * t, gb = gain_b * t;
        float v = a[i] * ga + b[i] * gb;
        if (v > 1.0f || v < -1.0f) { clipped++; if (v > 1.0f) v = 1.0f; else v = -1.0f; }
        out[i] = v;
        if (fabsf(v) > peak) peak = fabsf(v);
    }
    if (st) { st->peak = peak; st->clipped = clipped; st->samples = n; }
    return n;
}

/*
 * Streaming mixer. A mixer is length-neutral: n_out == n_in for every block,
 * so its length is entirely determined by the upstream resampler's unified
 * block lengths. The unified mode keeps one continuous gain ramp across block
 * boundaries instead of restarting the ramp (and the gain) each block.
 * Arithmetic mirrors mix_ramp() exactly, with global position replacing i.
 */
typedef struct {
    float gain_a, gain_b;
    int ramp_frames;
    int64_t pos;
    float peak;
    int clipped;
    int samples;
} MixRampStream;

MA_MAYBE_UNUSED static void mix_stream_init(MixRampStream* st, float gain_a, float gain_b,
                            int ramp_frames) {
    st->gain_a = gain_a;
    st->gain_b = gain_b;
    st->ramp_frames = ramp_frames;
    st->pos = 0;
    st->peak = 0.0f;
    st->clipped = 0;
    st->samples = 0;
}

MA_MAYBE_UNUSED static int mix_stream_block(MixRampStream* st, const float* a, const float* b,
                            float* out, int n) {
    if (n <= 0 || st->gain_a < 0.0f || st->gain_b < 0.0f ||
        st->ramp_frames < 0) return -1;
    for (int i = 0; i < n; i++) {
        float t = st->ramp_frames > 0
                      ? (float)st->pos / (float)st->ramp_frames
                      : 1.0f;
        if (t > 1.0f) t = 1.0f;
        float ga = st->gain_a * t, gb = st->gain_b * t;
        float v = a[i] * ga + b[i] * gb;
        if (v > 1.0f || v < -1.0f) {
            st->clipped++;
            if (v > 1.0f) v = 1.0f; else v = -1.0f;
        }
        out[i] = v;
        if (fabsf(v) > st->peak) st->peak = fabsf(v);
        st->pos++;
        st->samples++;
    }
    return n;
}

#ifndef MA_NO_MAIN
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
#endif /* MA_NO_MAIN */

#ifndef MA_NO_MAIN
int main(void) {
#define N 128
    float a[N], b[N], out[N];
    for (int i = 0; i < N; i++) { a[i] = 0.5f; b[i] = 0.25f; }

    MixStat st;
    check("mix_ok", mix_ramp(a, b, 1.0f, 1.0f, out, N, 0, &st) == N);
    check("mix_sum", fabsf(out[N - 1] - 0.75f) < 1e-5f);
    check("no_clip", st.clipped == 0);

    check("ramp_bad", mix_ramp(a, b, 1.0f, 1.0f, out, N, N + 1, &st) == -1);
    check("neg_gain", mix_ramp(a, b, -1.0f, 1.0f, out, N, 0, &st) == -1);

    float loud[N];
    for (int i = 0; i < N; i++) loud[i] = 0.9f;
    mix_ramp(loud, loud, 1.0f, 1.0f, out, N, 0, &st);
    check("clipped_detected", st.clipped == N && st.peak <= 1.0f + 1e-6f);

    /* ramp must start near zero */
    mix_ramp(a, b, 1.0f, 1.0f, out, N, N, &st);
    check("ramp_start_silent", fabsf(out[0]) < 1e-6f);

    /* streaming mixer must match one whole-buffer ramped mix byte for byte */
    float ref2[N], blk2[N];
    mix_ramp(a, b, 1.0f, 1.0f, ref2, N, N, NULL);
    MixRampStream ms;
    mix_stream_init(&ms, 1.0f, 1.0f, N);
    int half = N / 2;
    check("mix_stream_len",
          mix_stream_block(&ms, a, b, blk2, half) == half &&
          mix_stream_block(&ms, a + half, b + half, blk2 + half, half) == half);
    check("mix_stream_bitexact", memcmp(ref2, blk2, sizeof(blk2)) == 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

#endif /* MA_NO_MAIN */
