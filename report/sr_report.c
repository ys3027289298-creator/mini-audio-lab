/* sr_report: cross-validation report for the unified sample-rate assumption.
 *
 * Links the real DSP code from resampler/mixer/limiter (compiled with
 * -Dstatic= and per-module symbol renames, see Makefile) and emits four
 * datasets:
 *   1. length tables (ms -> frames, resample output lengths)
 *   2. per-stage waveform summaries (n/peak/rms/mean + FNV-1a hash of bytes)
 *   3. timing at 48000 Hz and 44100 Hz
 *   4. old-interface compatibility list (each old API called and verified)
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "sr.h"

/* ---- resampler (from ../resampler/src/resampler.c) ---- */
extern int resample_linear(const float* in, int n_in, double ratio, float* out, int cap);
extern int resample_rates(const float* in, int n_in, float sr_in, float sr_out, float* out, int cap);

/* ---- mixer (from ../mixer/src/mixer.c) ---- */
typedef struct { float peak; int clipped; int samples; } MixStat;
extern int mix_ramp(const float* a, const float* b, float gain_a, float gain_b,
                    float* out, int n, int ramp_frames, MixStat* st);
extern int mix_ramp_ms(const float* a, const float* b, float gain_a, float gain_b,
                       float* out, int n, float ramp_ms, float sr, MixStat* st);

/* ---- limiter (from ../limiter/src/limiter.c) ---- */
typedef struct { float attack_ms, release_ms, sample_rate, env; int engaged; } Limiter;
typedef struct { float peak_in, peak_out; int limited_samples; float max_reduction; } Report;
extern void lim_init(Limiter* l, float attack_ms, float release_ms, float sr);
extern void lim_init_default(Limiter* l, float attack_ms, float release_ms);
extern int lim_process(Limiter* l, const float* in, float* out, int n, float ceiling, Report* rep);

#define CAP 65536
static float buf_a[CAP], buf_b[CAP], buf_out[CAP];

static int failures = 0;
static void check(const char* name, int ok) {
    printf("  %-28s %s\n", name, ok ? "OK" : "FAIL");
    if (!ok) failures++;
}

static uint64_t fnv1a(const float* v, int n) {
    uint64_t h = 1469598103934665603ULL;
    const unsigned char* p = (const unsigned char*)v;
    for (size_t i = 0; i < (size_t)n * sizeof(float); i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}

static void summary(const char* stage, float sr, const float* v, int n) {
    double sum = 0.0, sq = 0.0;
    float pk = 0.0f;
    for (int i = 0; i < n; i++) {
        float a = fabsf(v[i]);
        if (a > pk) pk = a;
        sum += v[i];
        sq += (double)v[i] * v[i];
    }
    printf("  %-9s sr=%5.0f n=%5d peak=%.6f rms=%.6f mean=%+.6f hash=%016llx\n",
           stage, (double)sr, n, (double)pk, n > 0 ? sqrt(sq / n) : 0.0, n > 0 ? sum / n : 0.0,
           (unsigned long long)fnv1a(v, n));
}

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void tone(float* v, int n, float sr, float freq, float amp) {
    for (int i = 0; i < n; i++) v[i] = amp * sinf(2.0f * 3.14159265f * freq * (float)i / sr);
}

int main(void) {
    const float rates[2] = { MAL_SR_DEFAULT, MAL_SR_ALT };

    /* ---- 1. length tables ---- */
    printf("== 1. length table: ms -> frames (mal_ms_to_frames) ==\n");
    printf("  %8s %14s %14s %6s\n", "ms", "frames@48000", "frames@44100", "diff");
    const float durs[] = { 1, 2, 5, 10, 20, 50, 100, 250, 1000 };
    for (int i = 0; i < 9; i++) {
        int f48 = mal_ms_to_frames(durs[i], rates[0]);
        int f44 = mal_ms_to_frames(durs[i], rates[1]);
        printf("  %8.0f %14d %14d %6d\n", (double)durs[i], f48, f44, f48 - f44);
    }
    printf("== 1b. length table: resample output length (resample_rates) ==\n");
    printf("  %8s %16s %16s\n", "n_in", "out 48000>44100", "out 44100>48000");
    const int lens[] = { 128, 256, 1024, 4410, 4800, 11025, 12000 };
    for (int i = 0; i < 7; i++) {
        int n = lens[i];
        for (int j = 0; j < n; j++) buf_a[j] = 0.0f;
        int d = resample_rates(buf_a, n, rates[0], rates[1], buf_out, CAP);
        int u = resample_rates(buf_a, n, rates[1], rates[0], buf_out, CAP);
        printf("  %8d %16d %16d\n", n, d, u);
    }

    /* ---- 2. per-stage waveform summaries at both rates ---- */
    printf("== 2. waveform summary per stage (1 kHz tone, 0.25 s, amp 0.9) ==\n");
    for (int r = 0; r < 2; r++) {
        float sr = rates[r];
        int n = mal_ms_to_frames(250.0f, sr);
        tone(buf_a, n, sr, 1000.0f, 0.9f);

        int m = resample_linear(buf_a, n, 1.0, buf_out, CAP);
        summary("resampler", sr, buf_out, m);

        for (int i = 0; i < n; i++) buf_b[i] = 0.25f;
        MixStat st;
        m = mix_ramp_ms(buf_a, buf_b, 1.0f, 1.0f, buf_out, n, 10.0f, sr, &st);
        summary("mixer", sr, buf_out, m);
        printf("  %-9s sr=%5.0f clipped=%d st.peak=%.6f\n", "mixer", (double)sr, st.clipped, (double)st.peak);

        for (int i = 0; i < n; i++) buf_a[i] *= 2.0f; /* push over ceiling */
        Limiter l; Report rep;
        lim_init(&l, 5.0f, 50.0f, sr);
        m = lim_process(&l, buf_a, buf_out, n, 0.8f, &rep);
        summary("limiter", sr, buf_out, m);
        printf("  %-9s sr=%5.0f limited=%d max_reduction=%.6f peak_out=%.6f\n",
               "limiter", (double)sr, rep.limited_samples, (double)rep.max_reduction, (double)rep.peak_out);
    }
    printf("== 2b. waveform summary: resample 48000 -> 44100 -> 48000 round trip ==\n");
    {
        int n48 = mal_ms_to_frames(250.0f, rates[0]);
        tone(buf_a, n48, rates[0], 1000.0f, 0.9f);
        int n44 = resample_rates(buf_a, n48, rates[0], rates[1], buf_b, CAP);
        summary("rs_48to44", rates[1], buf_b, n44);
        int back = resample_rates(buf_b, n44, rates[1], rates[0], buf_out, CAP);
        summary("rs_44to48", rates[0], buf_out, back);
    }

    /* ---- 3. timing at both sample rates ---- */
    printf("== 3. timing per stage at both rates (0.25 s buffers, 2000 iterations) ==\n");
    printf("  %-9s %6s %8s %10s %12s\n", "stage", "sr", "frames", "total_ms", "ns_per_sample");
    const int iters = 2000;
    for (int r = 0; r < 2; r++) {
        float sr = rates[r];
        int n = mal_ms_to_frames(250.0f, sr);
        tone(buf_a, n, sr, 1000.0f, 0.9f);
        for (int i = 0; i < n; i++) buf_b[i] = 0.25f;

        double t0 = now_s();
        for (int it = 0; it < iters; it++) resample_linear(buf_a, n, 1.0, buf_out, CAP);
        double dt = now_s() - t0;
        printf("  %-9s %6.0f %8d %10.1f %12.1f\n", "resampler", (double)sr, n, dt * 1e3, dt * 1e9 / ((double)n * iters));

        MixStat st;
        t0 = now_s();
        for (int it = 0; it < iters; it++) mix_ramp_ms(buf_a, buf_b, 1.0f, 1.0f, buf_out, n, 10.0f, sr, &st);
        dt = now_s() - t0;
        printf("  %-9s %6.0f %8d %10.1f %12.1f\n", "mixer", (double)sr, n, dt * 1e3, dt * 1e9 / ((double)n * iters));

        Limiter l; Report rep;
        t0 = now_s();
        for (int it = 0; it < iters; it++) { lim_init(&l, 5.0f, 50.0f, sr); lim_process(&l, buf_a, buf_out, n, 0.8f, &rep); }
        dt = now_s() - t0;
        printf("  %-9s %6.0f %8d %10.1f %12.1f\n", "limiter", (double)sr, n, dt * 1e3, dt * 1e9 / ((double)n * iters));
    }

    /* ---- 4. old-interface compatibility list ---- */
    printf("== 4. old-interface compatibility list ==\n");
    printf("  %-10s %-34s %s\n", "module", "old api (unchanged)", "assumption before");
    printf("  %-10s %-34s %s\n", "resampler", "resample_linear(in,n,ratio,out,cap)", "bare ratio, no absolute rate");
    printf("  %-10s %-34s %s\n", "mixer", "mix_ramp(...,ramp_frames,st)", "ramp in frames");
    printf("  %-10s %-34s %s\n", "limiter", "lim_init(...,sr)+lim_process", "ms + caller-supplied sr");
    for (int i = 0; i < 256; i++) buf_a[i] = (float)(i % 7) / 7.0f;
    check("resample_linear ratio=2.0", resample_linear(buf_a, 256, 2.0, buf_out, CAP) == 128);
    check("resample_linear ratio=0.5", resample_linear(buf_a, 256, 0.5, buf_out, CAP) == 512);
    check("resample_linear bad ratio", resample_linear(buf_a, 256, 0.0, buf_out, CAP) == -1);
    for (int i = 0; i < 128; i++) { buf_a[i] = 0.5f; buf_b[i] = 0.25f; }
    MixStat st;
    check("mix_ramp frames=0", mix_ramp(buf_a, buf_b, 1.0f, 1.0f, buf_out, 128, 0, &st) == 128 && st.clipped == 0);
    check("mix_ramp frames=128", mix_ramp(buf_a, buf_b, 1.0f, 1.0f, buf_out, 128, 128, &st) == 128 && fabsf(buf_out[0]) < 1e-6f);
    check("mix_ramp bad ramp", mix_ramp(buf_a, buf_b, 1.0f, 1.0f, buf_out, 128, 129, &st) == -1);
    Limiter l; Report rep;
    lim_init(&l, 5.0f, 50.0f, 48000.0f);
    check("lim_init+lim_process @48000", lim_process(&l, buf_a, buf_out, 128, 0.9f, &rep) == 128 && rep.limited_samples == 0);
    check("lim_process bad ceiling", lim_process(&l, buf_a, buf_out, 128, 0.0f, &rep) == -1);
    printf("  note: default self-test output of every module is byte-identical to\n");
    printf("  pre-unification HEAD; verify with 'make -C report check-compat'.\n");

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
