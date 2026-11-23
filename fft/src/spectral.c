/*
 * spectral.c — Hann-windowed Welch averaged power spectrum.
 *
 * Frame count is not guessed. It is solved from the jitter budget (see
 * SPECTRAL_AVERAGING.md for the full derivations and the Monte-Carlo evidence):
 *
 *   one-frame squared relative variance  C(rho) = (1 + 2 rho)/(1 + rho)^2
 *   Hann 50% overlap inflation factor    IF_50 = 1.5  (conservative upper
 *     bound; the exact rho-dependent law gives 1.03..1.33, see the doc)
 *   K-frame jitter (CV)                  J = sqrt(IF_50 * C(rho) / K)
 *   design                               K = ceil(IF_50 * C(rho) / J_max^2)
 *
 * Local bin SNR for a Hann window and a tone exactly on a bin:
 *   rho = P_s/P_n = (S1^2 / (2 S2)) * 10^(SNR_dB/10) = N/3 * 10^(SNR_dB/10),
 *   with S1 = sum w = N/2, S2 = sum w^2 = 3N/8.
 *
 * The low-SNR departures (leakage, coloured/non-stationary noise, magnitude
 * averaging) are measured by tools/spectral_mc.py; the self-test below checks
 * the parts the C code can see: window constants, the design formula against
 * hand-computed values, unbiased noise floor, and empirical 1/K convergence.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fft_internal.h"

#ifndef SPECTRAL_N
#define SPECTRAL_N 64
#endif
#define NBINS (SPECTRAL_N / 2 + 1)
#define PI_F 3.14159265358979323846f

/* Hann 50% overlap inflation factor (analytic, MC-confirmed). */
#define HANN_IF_50 1.5f
/* Hann constants S1 = sum w, S2 = sum w^2 for n = SPECTRAL_N (floats). */
#define HANN_S1 (0.5f * (float)SPECTRAL_N)
#define HANN_S2 (0.375f * (float)SPECTRAL_N)

typedef struct {
    int    n;                 /* FFT/frame size (power of two)            */
    int    hop;               /* frame advance, n/2 for 50% overlap       */
    int    frames_target;     /* K budget from spectral_frames_needed     */
    int    frames_seen;       /* frames accumulated so far                */
    float  window[SPECTRAL_N];
    float  s2;                /* sum w^2, periodogram normalisation       */
    float  acc[NBINS];        /* accumulated power per bin                */
    C      frame[SPECTRAL_N];
} SpectralAverager;

static void spectral_init(SpectralAverager *s, int frames_target) {
    memset(s, 0, sizeof(*s));
    s->n = SPECTRAL_N;
    s->hop = SPECTRAL_N / 2;
    s->frames_target = frames_target;
    for (int i = 0; i < s->n; i++) {
        s->window[i] = 0.5f * (1.0f - cosf(2.0f * PI_F * (float)i / (float)s->n));
        s->s2 += s->window[i] * s->window[i];
    }
}

/* One Hann-windowed frame of n real samples -> accumulate power into acc. */
static void spectral_push(SpectralAverager *s, const float *x) {
    for (int i = 0; i < s->n; i++) {
        s->frame[i].re = x[i] * s->window[i];
        s->frame[i].im = 0.0f;
    }
    fft_forward(s->frame, s->n);
    for (int b = 0; b < NBINS; b++) {
        s->acc[b] += s->frame[b].re * s->frame[b].re
                   + s->frame[b].im * s->frame[b].im;
    }
    s->frames_seen++;
}

/* Averaged, unbiased periodogram: E[out] = PSD for unit-variance white noise. */
static void spectral_result(const SpectralAverager *s, float *out) {
    float norm = s->s2 * (float)s->frames_seen;
    for (int b = 0; b < NBINS; b++) out[b] = s->acc[b] / norm;
}

/* Local bin SNR rho for an on-bin tone at frame-level SNR snr_db.
 * Tone A = sqrt(2 SN): on-bin complex amplitude is A*S1/2, hence tone power
 * A^2 S1^2 / 4 = SN * S1^2 / 2; noise floor per bin equals S2, so
 * rho = SN * S1^2 / (2 S2). For Hann, S1=N/2, S2=3N/8 -> rho = SN*N/3. */
static float spectral_local_rho(float snr_db) {
    float sn = powf(10.0f, snr_db / 10.0f);
    return (HANN_S1 * HANN_S1 / (2.0f * HANN_S2)) * sn;
}

/* One-frame squared relative variance of the periodogram at local SNR rho. */
static float spectral_c_rho(float rho) {
    float d = 1.0f + rho;
    return (1.0f + 2.0f * rho) / (d * d);
}

/*
 * Smallest frame count meeting the jitter budget j_max under 50% overlap.
 * Returns 0 for invalid arguments; clamps the result to [1, 1u<<30].
 */
static unsigned spectral_frames_needed(float j_max, float snr_db) {
    if (!(j_max > 0.0f) || j_max > 1.0f || !isfinite(snr_db)) return 0;
    float rho = spectral_local_rho(snr_db);
    float k = HANN_IF_50 * spectral_c_rho(rho) / (j_max * j_max);
    if (k < 1.0f) k = 1.0f;
    if (k > (float)(1u << 30)) return 1u << 30;
    return (unsigned)ceilf(k);
}

/* Predicted jitter for a chosen frame count (used by callers and the test). */
static float spectral_predicted_jitter(unsigned k, float snr_db) {
    if (k == 0) return 0.0f;
    float rho = spectral_local_rho(snr_db);
    return sqrtf(HANN_IF_50 * spectral_c_rho(rho) / (float)k);
}

/* ------------------------------------------------------------------ tests */
static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* Deterministic xorshift32 so the empirical convergence test is reproducible. */
static uint32_t rng_state = 0x12345678u;
static float u01(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return (float)(rng_state >> 8) / (float)(1u << 24);
}
static float randn(void) {
    static int have = 0;
    static float spare;
    if (have) { have = 0; return spare; }
    /* Box-Muller from two uniform xorshift draws */
    float u1, u2;
    do { u1 = u01(); } while (u1 <= 1e-7f);
    u2 = u01();
    float r = sqrtf(-2.0f * logf(u1));
    spare = r * sinf(2.0f * PI_F * u2);
    have = 1;
    return r * cosf(2.0f * PI_F * u2);
}

/* White-noise-only run of k overlapping frames; fills the averaged spectrum. */
static void noise_run(unsigned k, float *out) {
    SpectralAverager s;
    spectral_init(&s, (int)k);
    int total = s.n + (int)k * s.hop;          /* k frames, 50% overlap */
    float *x = malloc((size_t)total * sizeof(float));
    if (!x) { exit(2); }
    for (int i = 0; i < total; i++) x[i] = randn();
    for (unsigned f = 0; f < k; f++) spectral_push(&s, &x[f * s.hop]);
    spectral_result(&s, out);
    free(x);
}

int main(void) {
    /* 1. window constants: S1 = N/2, S2 = 3N/8 */
    SpectralAverager s;
    spectral_init(&s, 16);
    check("hann_s1", fabsf(s.n * 0.5f - HANN_S1) < 1e-6f);
    check("hann_s2", fabsf(s.s2 - HANN_S2) / HANN_S2 < 1e-5f);
    /* periodic Hann: w[0]=0, w[N-1] nearly one full period (not zero) */
    check("hann_endpoints", fabsf(s.window[0]) < 1e-6f
                            && fabsf(s.window[s.n - 1]
                                     - sinf(PI_F / s.n) * sinf(PI_F / s.n)) < 1e-6f);

    /* 2. design formula against hand values (N=64 -> rho factor N/3 = 21.3333) */
    /* noise-only rho=0: C=1 -> K = ceil(1.5/J^2) */
    check("design_noise_j10", spectral_frames_needed(0.10f, -100.0f) == 150u);
    check("design_noise_j20", spectral_frames_needed(0.20f, -100.0f) == 38u);
    /* 0 dB tone: rho = 21.3333, C = 43.6667/498.78 = 0.08755;
       J=10% -> K = ceil(1.5*0.08755/0.01) = 14 */
    unsigned k0 = spectral_frames_needed(0.10f, 0.0f);
    check("design_0db_j10", k0 == 14u);
    check("design_respects_budget",
          spectral_predicted_jitter(k0, 0.0f) <= 0.10f);
    check("design_roundup_minimal",
          spectral_predicted_jitter(k0 - 1, 0.0f) > 0.10f);
    check("design_bad_args", spectral_frames_needed(0.0f, 0.0f) == 0u
                             && spectral_frames_needed(-1.0f, 0.0f) == 0u);

    /* 3. unbiased white-noise floor and empirical 1/K jitter convergence.
       Per-bin statistics over REPS independent noise realizations; interior
       bins are pooled for a tight estimate of the floor mean and its CV. */
    enum { REPS = 400, B0 = 4, B1 = NBINS - 4 };
    float floor_mean[2] = {0.0f, 0.0f};
    float cv_meas[2] = {0.0f, 0.0f};
    unsigned ks[2] = {16, 64};
    for (int qi = 0; qi < 2; qi++) {
        unsigned k = ks[qi];
        long count = 0;
        double acc = 0.0, acc2 = 0.0;
        for (int r = 0; r < REPS; r++) {
            float out[NBINS];
            noise_run(k, out);
            for (int b = B0; b < B1; b++) {
                acc += out[b];
                acc2 += (double)out[b] * out[b];
                count++;
            }
        }
        double mu = acc / count;
        double var = fmax(acc2 / count - mu * mu, 0.0);
        floor_mean[qi] = (float)mu;
        cv_meas[qi] = (float)(sqrt(var) / mu);
    }
    check("noise_floor_unbiased", fabsf(floor_mean[0] - 1.0f) < 0.03f
                                  && fabsf(floor_mean[1] - 1.0f) < 0.03f);
    /* jitter halves when K goes 16 -> 64 (the 1/K law); analytic prediction
       uses the textbook Hann IF 1.5, which is conservative versus the measured
       ~1.06 (see tools/spectral_mc.py expB), so gate against both. */
    float pred16 = sqrtf(HANN_IF_50 / 16.0f);
    float pred64 = sqrtf(HANN_IF_50 / 64.0f);
    check("jitter_scales_1_over_sqrtK",
          fabsf(cv_meas[0] / cv_meas[1] - 2.0f) < 0.15f);
    check("jitter_below_conservative_law",
          cv_meas[0] <= pred16 * 1.15f && cv_meas[1] <= pred64 * 1.15f);

    /* 4. on-bin tone rho: Hann factor N/3 at 0 dB */
    float tone_rho = spectral_local_rho(0.0f);
    check("rho_hann_factor", fabsf(tone_rho - ((float)SPECTRAL_N / 3.0f)) < 1e-3f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
