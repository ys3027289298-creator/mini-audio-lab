/*
 * pipeline - streaming demo (resampler -> mixer -> spatial) and timing bench.
 *
 * Build twice from the Makefile:
 *   pipeline_legacy   MA_UNIFY=0  original per-block length derivation
 *   pipeline_unified  MA_UNIFY=1  unified lengths from stream/stream.h
 *
 * Usage:
 *   pipeline demo <outdir>
 *   pipeline bench                       (CSV timing rows to stdout)
 *
 * The old per-stage function signatures are untouched; this file only drives
 * the original functions (legacy mode) or the *_stream_* variants (unified).
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MA_NO_MAIN
#include "resampler.c"
#include "mixer.c"
#include "spatial.c"

#ifndef MA_UNIFY
#define MA_UNIFY 1
#endif

#define SR_IN 8000
#define RATIO 1.3
#define DURATION_S 12
#define BLOCK_S 2
#define N_BLOCKS (DURATION_S / BLOCK_S)
#define TOTAL_IN (SR_IN * DURATION_S)
#define BLOCK_IN (SR_IN * BLOCK_S)

#define MIX_GAIN_A 0.6f
#define MIX_GAIN_B 0.6f
#define RAMP_BLOCK 2048          /* legacy: ramp restarts every block */
#define RAMP_TOTAL 4096          /* unified: one ramp over the stream head */
#define PAN_X (-0.3f)
#define DIST_REF 2.0f
#define DIST_ROLLOFF 1.0f
#define DIST_NOW 4.0f

static void make_input(float* x, int n, int sr) {
    for (int i = 0; i < n; i++) {
        double t = (double)i / sr;
        x[i] = (float)(0.55 * sin(2.0 * M_PI * 220.0 * t)
                     + 0.30 * sin(2.0 * M_PI * 660.0 * t));
    }
}

static void write_f32(const char* path, const float* x, int n) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); exit(2); }
    if (n > 0 && fwrite(x, sizeof(float), (size_t)n, f) != (size_t)n) {
        perror(path); exit(2);
    }
    fclose(f);
}

static int64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1000000000LL + ts.tv_nsec;
}

static int cmp_double(const void* a, const void* b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

static double median_ms(double* secs, int iters) {
    qsort(secs, (size_t)iters, sizeof(double), cmp_double);
    return secs[iters / 2] / 1e6;
}

/* whole-stream references for the byte-exactness checks */
static void ref_resample(const float* in, int n, double ratio, float* out) {
    resample_linear(in, n, ratio, out, n * 4);
}

static void ref_mix(const float* a, const float* b, int n, float* out) {
    mix_ramp(a, b, MIX_GAIN_A, MIX_GAIN_B, out, n, RAMP_TOTAL, NULL);
}

static void ref_spatial(const float* in, int n, Stereo* out, float att) {
    int pairs = n / 2;
    for (int p = 0; p < pairs; p++) {
        spatial_frame(in[2 * p], PAN_X, att, &out[p]);
        Stereo r;
        spatial_frame(in[2 * p + 1], PAN_X, att, &r);
        out[p].left += r.left;
        out[p].right += r.right;
    }
}

static int run_demo(const char* outdir) {
    float* input = malloc((size_t)(TOTAL_IN + 1) * sizeof(float));
    make_input(input, TOTAL_IN, SR_IN);
    input[TOTAL_IN] = input[TOTAL_IN - 1]; /* duplicated last sample */

    int cap = (int)((double)TOTAL_IN / RATIO) + 16;
    float *rs_out = malloc((size_t)cap * sizeof(float));
    float *mix_out = malloc((size_t)cap * sizeof(float));
    Stereo *sp_out = malloc((size_t)cap * sizeof(Stereo));
    float *sp_L = malloc((size_t)cap * sizeof(float));
    float *sp_R = malloc((size_t)cap * sizeof(float));
    float *sp_il = malloc((size_t)cap * 2 * sizeof(float));
    if (!input || !rs_out || !mix_out || !sp_out || !sp_L || !sp_R || !sp_il) {
        fprintf(stderr, "oom\n");
        return 2;
    }

    int blk_rs[N_BLOCKS], blk_mix[N_BLOCKS], blk_sp[N_BLOCKS];
    int n_rs = 0, n_mix = 0, n_sp_pairs = 0, sp_spurious = 0;

    float att = attenuate(DIST_NOW, DIST_REF, DIST_ROLLOFF);

#if MA_UNIFY
    ResampleStream rs; resample_stream_init(&rs, RATIO);
    MixRampStream ms; mix_stream_init(&ms, MIX_GAIN_A, MIX_GAIN_B, RAMP_TOTAL);
    SpatialStream ss; spatial_stream_init(&ss, PAN_X, att);
#else
    /* legacy mode: no cross-block state, each block is self-contained */
#endif

    for (int b = 0; b < N_BLOCKS; b++) {
        const float* inb MA_MAYBE_UNUSED = input + (int64_t)b * BLOCK_IN;

#if MA_UNIFY
        int m = resample_stream_block(&rs, input, BLOCK_IN, b,
                                      rs_out + n_rs, cap - n_rs);
        blk_rs[b] = m;
        n_rs += m;

        int k = mix_stream_block(&ms, rs_out + n_mix, rs_out + n_mix,
                                 mix_out + n_mix, m);
        blk_mix[b] = k;
        n_mix += k;

        int p = spatial_stream_block(&ss, mix_out + (int64_t)n_mix - k, k,
                                     sp_out + n_sp_pairs);
        blk_sp[b] = p;
        n_sp_pairs += p;
#else
        /* legacy: every stage rounds its own block independently */
        int m = resample_linear(inb, BLOCK_IN, RATIO,
                                rs_out + n_rs, cap - n_rs);
        blk_rs[b] = m;
        n_rs += m;

        MixStat st;
        int k = mix_ramp(rs_out + n_mix, rs_out + n_mix, MIX_GAIN_A,
                         MIX_GAIN_B, mix_out + n_mix, m, RAMP_BLOCK, &st);
        blk_mix[b] = k;
        n_mix += k;

        SpatialStream ls; spatial_stream_init(&ls, PAN_X, att);
        int p = spatial_legacy_block(&ls, mix_out + (int64_t)n_mix - k, k,
                                     sp_out + n_sp_pairs);
        blk_sp[b] = p;
        if (k & 1)
            sp_spurious++; /* dangling frame becomes a spurious stereo frame */
        n_sp_pairs += p;
#endif
    }

    for (int p = 0; p < n_sp_pairs; p++) {
        sp_L[p] = sp_out[p].left;
        sp_R[p] = sp_out[p].right;
        sp_il[2 * p] = sp_out[p].left;
        sp_il[2 * p + 1] = sp_out[p].right;
    }

    /* ---- byte-exactness against the appropriate reference ---- */
    float *ref_rs = malloc((size_t)cap * sizeof(float));
    float *ref_mixbuf = malloc((size_t)cap * sizeof(float));
    Stereo *ref_sp = malloc((size_t)cap * sizeof(Stereo));
    int ref_n_rs = (int)stream_total_out(TOTAL_IN, RATIO);
    ref_resample(input, TOTAL_IN, RATIO, ref_rs);
    ref_mix(ref_rs, ref_rs, ref_n_rs, ref_mixbuf);
    ref_spatial(ref_mixbuf, ref_n_rs, ref_sp, att);
    int ref_pairs = (int)stream_total_pairs(ref_n_rs);

    int rs_exact, mix_exact, sp_exact;
#if MA_UNIFY
    rs_exact = n_rs == ref_n_rs &&
               memcmp(rs_out, ref_rs, (size_t)n_rs * sizeof(float)) == 0;
    mix_exact = n_mix == ref_n_rs &&
                memcmp(mix_out, ref_mixbuf, (size_t)n_mix * sizeof(float)) == 0;
    sp_exact = n_sp_pairs == ref_pairs &&
               memcmp(sp_out, ref_sp,
                      (size_t)n_sp_pairs * sizeof(Stereo)) == 0;
#else
    /* legacy reference: replay the same per-block calls into fresh buffers */
    float *l_rs = malloc((size_t)cap * sizeof(float));
    float *l_mix = malloc((size_t)cap * sizeof(float));
    Stereo *l_sp = malloc((size_t)cap * sizeof(Stereo));
    int q1 = 0, q2 = 0, q3 = 0;
    for (int b = 0; b < N_BLOCKS; b++) {
        const float* inb = input + (int64_t)b * BLOCK_IN;
        int m = resample_linear(inb, BLOCK_IN, RATIO, l_rs + q1, cap - q1);
        q1 += m;
        MixStat st;
        int k = mix_ramp(l_rs + q2, l_rs + q2, MIX_GAIN_A, MIX_GAIN_B,
                         l_mix + q2, m, RAMP_BLOCK, &st);
        q2 += k;
        SpatialStream ls; spatial_stream_init(&ls, PAN_X, att);
        q3 += spatial_legacy_block(&ls, l_mix + q2 - k, k, l_sp + q3);
    }
    rs_exact = q1 == n_rs &&
               memcmp(l_rs, rs_out, (size_t)n_rs * sizeof(float)) == 0;
    mix_exact = q2 == n_mix &&
                memcmp(l_mix, mix_out, (size_t)n_mix * sizeof(float)) == 0;
    sp_exact = q3 == n_sp_pairs &&
               memcmp(l_sp, sp_out,
                      (size_t)n_sp_pairs * sizeof(Stereo)) == 0;
    free(l_rs); free(l_mix); free(l_sp);
#endif

    char path[512];
    snprintf(path, sizeof(path), "%s/input.f32", outdir);
    write_f32(path, input, TOTAL_IN);
    snprintf(path, sizeof(path), "%s/resampler.f32", outdir);
    write_f32(path, rs_out, n_rs);
    snprintf(path, sizeof(path), "%s/mixer.f32", outdir);
    write_f32(path, mix_out, n_mix);
    snprintf(path, sizeof(path), "%s/spatial_L.f32", outdir);
    write_f32(path, sp_L, n_sp_pairs);
    snprintf(path, sizeof(path), "%s/spatial_R.f32", outdir);
    write_f32(path, sp_R, n_sp_pairs);
    snprintf(path, sizeof(path), "%s/spatial_interleaved.f32", outdir);
    write_f32(path, sp_il, n_sp_pairs * 2);

    snprintf(path, sizeof(path), "%s/run.json", outdir);
    FILE* j = fopen(path, "w");
    if (!j) { perror(path); return 2; }
    fprintf(j, "{\n");
    fprintf(j, "  \"mode\": \"%s\",\n", MA_UNIFY ? "unified" : "legacy");
    fprintf(j, "  \"sr_in\": %d, \"ratio\": %.6f, \"sr_out\": %.1f,\n",
            SR_IN, RATIO, SR_IN / RATIO);
    fprintf(j, "  \"duration_s\": %d, \"block_s\": %d, \"n_blocks\": %d,\n",
            DURATION_S, BLOCK_S, N_BLOCKS);
    fprintf(j, "  \"total\": {\"input\": %d, \"resampler\": %d, \"mixer\": %d, "
            "\"spatial_pairs\": %d, \"spatial_frames\": %d, "
            "\"spatial_spurious_pairs\": %d},\n",
            TOTAL_IN, n_rs, n_mix, n_sp_pairs, n_sp_pairs * 2, sp_spurious);
    fprintf(j, "  \"block_len\": {\n    \"resampler\": [");
    for (int b = 0; b < N_BLOCKS; b++)
        fprintf(j, "%s%d", b ? "," : "", blk_rs[b]);
    fprintf(j, "],\n    \"mixer\": [");
    for (int b = 0; b < N_BLOCKS; b++)
        fprintf(j, "%s%d", b ? "," : "", blk_mix[b]);
    fprintf(j, "],\n    \"spatial_pairs\": [");
    for (int b = 0; b < N_BLOCKS; b++)
        fprintf(j, "%s%d", b ? "," : "", blk_sp[b]);
    fprintf(j, "]\n  },\n");
    fprintf(j, "  \"unified_formula_total_out\": %lld,\n",
            (long long)stream_total_out(TOTAL_IN, RATIO));
    fprintf(j, "  \"unified_formula_total_pairs\": %d,\n", ref_pairs);
    fprintf(j, "  \"byte_exact\": {\"resampler\": %d, \"mixer\": %d, "
            "\"spatial\": %d},\n", rs_exact, mix_exact, sp_exact);
    fprintf(j, "  \"all_byte_exact\": %d\n}\n",
            rs_exact && mix_exact && sp_exact);
    fclose(j);

    printf("mode=%s input=%d resampler=%d mixer=%d spatial_pairs=%d "
           "spurious_pairs=%d byte_exact=%d/%d/%d\n",
           MA_UNIFY ? "unified" : "legacy", TOTAL_IN, n_rs, n_mix,
           n_sp_pairs, sp_spurious, rs_exact, mix_exact, sp_exact);

    free(input); free(rs_out); free(mix_out); free(sp_out);
    free(sp_L); free(sp_R); free(sp_il);
    free(ref_rs); free(ref_mixbuf); free(ref_sp);
    return (rs_exact && mix_exact && sp_exact) ? 0 : 3;
}

#define BENCH_ITERS 100
#define BENCH_BLOCK 1024

/* volatile sink so -O2 cannot eliminate a stage whose output is unused.
 * noinline checksum over the freshly written buffer: every stage output is
 * observed through a call GCC cannot see through. */
static volatile float g_sink;

__attribute__((noinline))
static void consume(const float* p, int64_t n) {
    float acc = 0.0f;
    for (int64_t i = 0; i < n; i++) acc += p[i];
    g_sink += acc;
}

__attribute__((noinline))
static void consume_stereo(const Stereo* p, int64_t n) {
    float acc = 0.0f;
    for (int64_t i = 0; i < n; i++) acc += p[i].left + p[i].right;
    g_sink += acc;
}

static void bench_stage(const char* stage, int sr_in, int sr_out,
                        int work_in) {
    double ratio = (double)sr_in / sr_out; /* output grows when sr_out > sr_in */
    int n_out = (int)stream_total_out(work_in, ratio);
    float* a = malloc(((size_t)work_in + 1) * sizeof(float));
    float* o = malloc(((size_t)n_out + 8) * sizeof(float));
    Stereo* s = malloc((size_t)n_out * sizeof(Stereo));
    make_input(a, work_in, sr_in);
    a[work_in] = a[work_in - 1];

    double secs[BENCH_ITERS];

    for (int it = 0; it < BENCH_ITERS; it++) {
        int64_t t0 = now_ns();
        if (!strcmp(stage, "resampler")) {
            int blocks = work_in / BENCH_BLOCK;
#if MA_UNIFY
            ResampleStream rs; resample_stream_init(&rs, ratio);
            int got = 0;
            for (int bb = 0; bb < blocks; bb++)
                got += resample_stream_block(&rs, a, BENCH_BLOCK, bb,
                                             o + got, n_out + 8 - got);
#else
            int got = 0;
            for (int bb = 0; bb < blocks; bb++)
                got += resample_linear(a + (int64_t)bb * BENCH_BLOCK,
                                       BENCH_BLOCK, ratio,
                                       o + got, n_out + 8 - got);
#endif
            consume(o, n_out);
        } else if (!strcmp(stage, "mixer")) {
            int blocks = n_out / BENCH_BLOCK;
#if MA_UNIFY
            MixRampStream ms;
            mix_stream_init(&ms, MIX_GAIN_A, MIX_GAIN_B, RAMP_TOTAL);
            int got = 0;
            for (int bb = 0; bb < blocks; bb++)
                got += mix_stream_block(&ms, o + got, o + got,
                                        o + got, BENCH_BLOCK);
#else
            int got = 0;
            MixStat st;
            for (int bb = 0; bb < blocks; bb++) {
                mix_ramp(o + got, o + got, MIX_GAIN_A, MIX_GAIN_B,
                         o + got, BENCH_BLOCK, RAMP_BLOCK, &st);
                got += BENCH_BLOCK;
            }
#endif
            consume(o, n_out);
        } else { /* spatial */
            float att = attenuate(DIST_NOW, DIST_REF, DIST_ROLLOFF);
            int blocks = n_out / BENCH_BLOCK;
#if MA_UNIFY
            SpatialStream ss; spatial_stream_init(&ss, PAN_X, att);
            int got = 0;
            for (int bb = 0; bb < blocks; bb++) {
                got += spatial_stream_block(&ss, o + (int64_t)bb * BENCH_BLOCK,
                                            BENCH_BLOCK, s + got);
            }
#else
            int got = 0;
            for (int bb = 0; bb < blocks; bb++) {
                SpatialStream ls; spatial_stream_init(&ls, PAN_X, att);
                got += spatial_legacy_block(&ls,
                        o + (int64_t)bb * BENCH_BLOCK, BENCH_BLOCK, s + got);
            }
#endif
            consume_stereo(s, n_out / 2);
        }
        secs[it] = (double)(now_ns() - t0);
    }
    double med = median_ms(secs, BENCH_ITERS);
    long long processed = !strcmp(stage, "spatial")
                              ? (long long)n_out * 2
                              : (long long)n_out;
    double msps = (double)processed / (med / 1000.0) / 1e6;
    printf("%s,%s,%d,%d,%d,%d,%.6f,%.3f\n",
           MA_UNIFY ? "unified" : "legacy", stage, sr_in, sr_out,
           n_out, BENCH_ITERS, med, msps);

    free(a); free(o); free(s);
}

static int run_bench(void) {
    printf("mode,stage,sr_in,sr_out,total_samples,iterations,median_ms,"
           "throughput_Msps\n");
    /* 12 s workload; sample rate doubles between the two rows */
    bench_stage("resampler", 8000, 16000, 8000 * DURATION_S);
    bench_stage("resampler", 16000, 32000, 16000 * DURATION_S);
    bench_stage("mixer", 8000, 16000, 8000 * DURATION_S);
    bench_stage("mixer", 16000, 32000, 16000 * DURATION_S);
    bench_stage("spatial", 8000, 16000, 8000 * DURATION_S);
    bench_stage("spatial", 16000, 32000, 16000 * DURATION_S);
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 2 && !strcmp(argv[1], "bench"))
        return run_bench();
    if (argc >= 3 && !strcmp(argv[1], "demo"))
        return run_demo(argv[2]);
    fprintf(stderr, "usage: %s demo <outdir> | bench\n", argv[0]);
    return 2;
}
