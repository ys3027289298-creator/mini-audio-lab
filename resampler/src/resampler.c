
/* resampler: linear resampling with output-length derivation. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "stream.h"

static int resample_linear(const float* in, int n_in, double ratio, float* out, int cap) {
    if (n_in <= 0 || ratio <= 0.0) return -1;
    int n_out = (int)floor((double)n_in / ratio);
    if (n_out <= 0 || n_out > cap) return -1;
    for (int i = 0; i < n_out; i++) {
        double src = (double)i * ratio;
        int i0 = (int)floor(src);
        double t = src - (double)i0;
        int i1 = i0 + 1 < n_in ? i0 + 1 : n_in - 1;
        out[i] = (float)((1.0 - t) * in[i0] + t * in[i1]);
    }
    return n_out;
}

/*
 * Streaming variant. `in` is the whole stream: n_in valid samples plus one
 * duplicated last sample at index n_in (mirroring resample_linear()'s final
 * i1 clamp). Passing the full stream keeps interpolation indexing identical
 * to a single-buffer call, including the interp sample straddling a block
 * boundary. n_emit comes from stream_block_out(), so all blocks together
 * produce exactly stream_total_out(total_in, ratio).
 */
typedef struct {
    double ratio;
    int64_t produced;     /* total output samples emitted */
} ResampleStream;

MA_MAYBE_UNUSED static void resample_stream_init(ResampleStream* st, double ratio) {
    st->ratio = ratio;
    st->produced = 0;
}

MA_MAYBE_UNUSED static int resample_stream_block(ResampleStream* st, const float* in,
                                 int block_in, int block,
                                 float* out, int cap) {
    if (block_in <= 0 || st->ratio <= 0.0) return -1;
    int64_t n_emit = stream_block_out(block, block_in, st->ratio);
    if (n_emit <= 0 || n_emit > cap) return -1;
    for (int64_t i = 0; i < n_emit; i++) {
        /* identical arithmetic order to resample_linear() */
        double src = (double)(st->produced + i) * st->ratio;
        int64_t i0 = (int64_t)floor(src);
        double t = src - (double)i0;
        out[i] = (float)((1.0 - t) * in[i0] + t * in[i0 + 1]);
    }
    st->produced += n_emit;
    return (int)n_emit;
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
#define N 256
    float in[N], out[N * 4];
    for (int i = 0; i < N; i++) in[i] = (float)(i % 7) / 7.0f;
    int n = resample_linear(in, N, 2.0, out, N * 4);
    check("downsample_len", n == N / 2);
    n = resample_linear(in, N, 0.5, out, N * 4);
    check("upsample_len", n == N * 2);
    check("bad_ratio", resample_linear(in, N, 0.0, out, N * 4) == -1);
    check("cap_too_small", resample_linear(in, N, 0.1, out, 4) == -1);

    float dc[N];
    for (int i = 0; i < N; i++) dc[i] = 0.25f;
    n = resample_linear(dc, N, 1.3, out, N * 4);
    check("dc_preserved", fabsf(out[n / 2] - 0.25f) < 1e-5f);

    float ramp[N];
    for (int i = 0; i < N; i++) ramp[i] = (float)i;
    n = resample_linear(ramp, N, 2.0, out, N * 4);
    check("monotone_ramp", out[1] > out[0] && out[3] > out[2]);

    /* streaming blocks must reproduce the whole-buffer reference exactly */
    float ref[N * 4];
    int nref = resample_linear(ramp, N, 1.3, ref, N * 4);
    float full[N + 1];
    memcpy(full, ramp, N * sizeof(float));
    full[N] = ramp[N - 1]; /* duplicated last sample for boundary interp */
    float blk[N * 4];
    ResampleStream rs;
    resample_stream_init(&rs, 1.3);
    int got = 0, bn = 32; /* N is a whole multiple of the block size */
    for (int b = 0; b < N / bn; b++) {
        int m = resample_stream_block(&rs, full, bn, b, blk + got, N * 4 - got);
        check("stream_len", m == (int)stream_block_out(b, bn, 1.3));
        got += m;
    }
    check("stream_total_len", got == nref &&
          got == (int)stream_total_out(N, 1.3));
    check("stream_bitexact", memcmp(blk, ref, (size_t)got * sizeof(float)) == 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

#endif /* MA_NO_MAIN */
