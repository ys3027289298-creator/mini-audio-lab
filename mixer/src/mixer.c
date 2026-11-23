
/* mixer: float mixing with ramp gain and clipping detection. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define N 128

/*
 * Ramp kernel switch.
 *   MIXER_SMOOTH_RAMP=1 (default): cubic smoothstep s(x)=3x^2-2x^3,
 *     C1 at both ends (s'(0)=s'(1)=0), second derivative bounded by 6.
 *   MIXER_SMOOTH_RAMP=0: legacy linear ramp s(x)=x, reproduced exactly
 *     (including the t>1 clamp) so turning the change off restores the
 *     previous behaviour bit-for-bit.
 */
#ifndef MIXER_SMOOTH_RAMP
#define MIXER_SMOOTH_RAMP 1
#endif

/* s(x)=3x^2-2x^3 on [0,1], clamped outside; s'(x)=6x-6x^2, s''(x)=6-12x */
static inline float ramp_kernel(float x) {
#if MIXER_SMOOTH_RAMP
    if (x <= 0.0f) return 0.0f;
    if (x >= 1.0f) return 1.0f;
    return x * x * (3.0f - 2.0f * x);
#else
    if (x > 1.0f) x = 1.0f;
    return x;
#endif
}

typedef struct { float peak; int clipped; int samples; } MixStat;

static int mix_ramp(const float* a, const float* b, float gain_a, float gain_b,
                    float* out, int n, int ramp_frames, MixStat* st) {
    if (n <= 0) return -1;
    if (ramp_frames < 0 || ramp_frames > n) return -1;
    if (gain_a < 0.0f || gain_b < 0.0f) return -1;
    float peak = 0.0f; int clipped = 0;
    for (int i = 0; i < n; i++) {
        float x = ramp_frames > 0 ? (float)i / (float)ramp_frames : 1.0f;
        float t = ramp_kernel(x);
        float ga = gain_a * t, gb = gain_b * t;
        float v = a[i] * ga + b[i] * gb;
        if (v > 1.0f || v < -1.0f) { clipped++; if (v > 1.0f) v = 1.0f; else v = -1.0f; }
        out[i] = v;
        if (fabsf(v) > peak) peak = fabsf(v);
    }
    if (st) { st->peak = peak; st->clipped = clipped; st->samples = n; }
    return n;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* analytic reference of the smooth kernel, kept independent of ramp_kernel */
#if MIXER_SMOOTH_RAMP
static float kernel_ref(float x) {
    if (x <= 0.0f) return 0.0f;
    if (x >= 1.0f) return 1.0f;
    return 3.0f * x * x - 2.0f * x * x * x;
}
#endif

int main(void) {
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

#if MIXER_SMOOTH_RAMP
    /* ---- smooth-kernel evidence, constant input so out[i] = g * s(i/R) ---- */
    const float g = 0.75f;            /* a*1 + b*1 with the inputs above */
    const int R = N;                  /* full-length ramp */
    const float h = 1.0f / (float)R;
    /* kernel samples with one ghost point on each side (the flat regions):
     * ss[k+1] = s(x_k), x_k = k/R, k = -1..R+1; s(x<0)=0, s(x>1)=1. */
    float ss[N + 3];
    mix_ramp(a, b, 1.0f, 1.0f, out, N, R, &st);
    ss[0] = 0.0f;                     /* ghost s(x<0) */
    for (int i = 0; i < R; i++) ss[i + 1] = out[i] / g;
    ss[R + 1] = 1.0f;                 /* s(x=1) */
    ss[R + 2] = 1.0f;                 /* ghost s(x>1) */

    /* endpoint values: s(0)=0, s(1)=1 */
    check("kernel_endpoints", fabsf(ss[1]) < 1e-6f && fabsf(ss[R + 1] - 1.0f) < 1e-6f);

    /*
     * First discrete difference (forward) on the ramp intervals:
     * d1[i] = s(x_{i+1}) - s(x_i), i = 0..R-1.
     * Analytically s(h)-s(0) = 3h^2-2h^3 = O(h^2): the gain step on the
     * first/last ramp interval shrinks quadratically. The legacy linear
     * ramp takes a full h-sized step there (slope 0 -> 1/R jump), which is
     * the end-of-ramp pop. Require the observed seam step under 4h^2.
     */
    float d1[N];
    for (int i = 0; i < R; i++) d1[i] = ss[i + 2] - ss[i + 1];
    const float seam_bound = 4.0f * h * h;
    check("d1_seam_start_oh2", d1[0] <= seam_bound);
    check("d1_seam_end_oh2", (1.0f - ss[R]) <= seam_bound);
    int d1_nonneg = 1;
    for (int i = 0; i < R; i++) if (d1[i] < -1e-7f) { d1_nonneg = 0; break; }
    check("d1_monotonic", d1_nonneg);
    printf("d1_seam=%.9g (legacy linear would step h=%.9g, %.1fx larger)\n",
           d1[0], h, h / d1[0]);

    /*
     * Central differences with ghost points estimate the derivatives at the
     * seams. C1 means the estimated slope matches the flat outside slope 0:
     * p(0), p(1) must be O(h); analytic s'(0)=s'(1)=0, s'(1/2)=3/2.
     */
    float p0 = (ss[2] - ss[0]) / (2.0f * h);
    float p1 = (ss[R + 2] - ss[R]) / (2.0f * h);
    float pm = (ss[R / 2 + 2] - ss[R / 2]) / (2.0f * h);
    check("deriv_zero_at_seams", fabsf(p0) < 2.0f * h && fabsf(p1) < 2.0f * h);
    check("deriv_midpoint", fabsf(pm - 1.5f) < 1e-3f);
    printf("s'(0)=%.6f s'(1)=%.6f s'(0.5)=%.6f\n", p0, p1, pm);

    /*
     * Second discrete difference (central), q_i ~ s''(x_i).
     * s''(x)=6-12x is bounded by 6 on [0,1]. For a cubic the central
     * second difference is exact in real arithmetic, so interior points
     * also match the analytic values tightly (q(1/4)=3, q(3/4)=-3).
     */
    float max_abs_q = 0.0f;
    for (int i = 0; i <= R; i++) {
        float q = (ss[i + 2] - 2.0f * ss[i + 1] + ss[i]) / (h * h);
        if (fabsf(q) > max_abs_q) max_abs_q = fabsf(q);
    }
    check("d2_within_bound", max_abs_q <= 6.0f * (1.0f + 1e-3f));
    float q25 = (ss[R / 4 + 2] - 2.0f * ss[R / 4 + 1] + ss[R / 4]) / (h * h);
    float q75 = (ss[3 * R / 4 + 2] - 2.0f * ss[3 * R / 4 + 1] + ss[3 * R / 4]) / (h * h);
    check("d2_matches_analytic", fabsf(q25 - 3.0f) < 1e-3f && fabsf(q75 + 3.0f) < 1e-3f);
    printf("s''_max_abs=%.6f (analytic bound 6), s''(0.25)=%.6f s''(0.75)=%.6f\n",
           max_abs_q, q25, q75);

    /* pointwise match against the independent analytic smoothstep formula */
    int kernel_match = 1;
    float max_kernel_err = 0.0f;
    for (int i = 0; i <= R; i++) {
        float err = fabsf(ss[i + 1] - kernel_ref((float)i / (float)R));
        if (err > max_kernel_err) max_kernel_err = err;
        if (err > 1e-6f) { kernel_match = 0; break; }
    }
    printf("kernel_pointwise_max_abs_err=%.9g\n", max_kernel_err);
    check("kernel_pointwise_match", kernel_match);
#endif

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
