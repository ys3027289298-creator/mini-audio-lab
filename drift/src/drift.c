
/* drift: cross-correlation clock-drift estimation between two sample clocks.
 *
 * Clock A (reference) and clock B both nominally run at FS Hz, but B runs
 * fast/slow by a few ppm, so over long runs the sample counts disagree.
 * Windowed cross-correlation locates the per-segment offset of B relative
 * to A; a linear fit over the offset sequence gives the drift rate (ppm)
 * and the cumulative sample difference.
 *
 * enable=0 disables estimation: drift_estimate() does nothing and
 * drift_apply() is a bit-identical passthrough (output == current input).
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define FS 48000        /* nominal sample rate, both clocks */
#define SEG 4096        /* correlation segment length */
#define MAXLAG 128      /* offset search window: +-MAXLAG samples */
#define TRUE_PPM 50.0   /* simulated clock-B drift for the self-test */

typedef struct {
    double rate_ppm;    /* estimated drift rate; + means clock B runs fast */
    double cum_samples; /* cumulative sample difference at end of stream */
    int segments;       /* number of offset measurements */
    int valid;          /* 0 when estimation was disabled or failed */
} DriftEstimate;

/* deterministic noise so tests are reproducible */
static unsigned rng_state = 12345u;
static float frand(void) {
    rng_state = rng_state * 1664525u + 1013904223u;
    return (float)((rng_state >> 8) & 0xFFFF) / 32768.0f - 1.0f;
}

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

/* Two-clock stream pair: ref = master noise (clock A timeline); clock B runs
 * fast by (1 + ppm*1e-6), so B's n-th sample lands at master[n/(1+e)] and B
 * accumulates n*e extra samples relative to A. Arrays must hold n+MAXLAG+4. */
static void make_streams(float* ref, float* b, int n, double ppm) {
    for (int i = 0; i < n + MAXLAG + 4; i++) ref[i] = frand();
    double e = ppm * 1e-6;
    for (int i = 0; i < n + MAXLAG + 2; i++) {
        double src = (double)i / (1.0 + e);
        int i0 = (int)floor(src);
        double t = src - (double)i0;
        b[i] = (float)((1.0 - t) * ref[i0] + t * ref[i0 + 1]);
    }
}

/* Offset of b relative to ref at segment start s, via cross-correlation over
 * +-MAXLAG with parabolic refinement of the peak. Positive offset means B is
 * ahead (has run fast) by that many samples. The b segment is gathered into a
 * padded window with out-of-range samples zeroed, so no headroom is needed. */
static double segment_offset(const float* ref, const float* b, int n, int s) {
    float win[SEG + 2 * MAXLAG];
    for (int k = 0; k < SEG + 2 * MAXLAG; k++) {
        int idx = s - MAXLAG + k;
        win[k] = (idx >= 0 && idx < n) ? b[idx] : 0.0f;
    }
    double corr[2 * MAXLAG + 1];
    int best_l = -MAXLAG;
    for (int l = -MAXLAG; l <= MAXLAG; l++) {
        double c = 0.0;
        for (int i = 0; i < SEG; i++) c += (double)ref[s + i] * win[MAXLAG + i + l];
        corr[l + MAXLAG] = c;
        if (c > corr[best_l + MAXLAG]) best_l = l;
    }
    if (best_l > -MAXLAG && best_l < MAXLAG) {
        double y1 = corr[best_l - 1 + MAXLAG];
        double y2 = corr[best_l + MAXLAG];
        double y3 = corr[best_l + 1 + MAXLAG];
        double den = y1 - 2.0 * y2 + y3;
        if (den != 0.0) return best_l + 0.5 * (y1 - y3) / den;
    }
    return best_l;
}

/* Estimate drift from the offset sequence. enable=0: no work, invalid result. */
static int drift_estimate(const float* ref, const float* b, int n,
                          double* offsets, double* times,
                          DriftEstimate* est, int enable) {
    memset(est, 0, sizeof(*est));
    if (!enable) return 0;
    int nseg = n / SEG;
    if (nseg < 2) return -1;
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (int k = 0; k < nseg; k++) {
        int s = k * SEG;
        double t = (double)s / FS;
        double o = segment_offset(ref, b, n, s);
        if (offsets) { offsets[k] = o; times[k] = t; }
        sx += t; sy += o; sxx += t * t; sxy += t * o;
    }
    double denom = (double)nseg * sxx - sx * sx;
    if (denom == 0.0) return -1;
    double slope = ((double)nseg * sxy - sx * sy) / denom; /* samples/second */
    est->rate_ppm = slope / FS * 1e6;
    est->cum_samples = slope * (double)n / FS;
    est->segments = nseg;
    est->valid = 1;
    return 0;
}

/* Map clock-B stream back to the clock-A timeline using the estimated rate.
 * enable=0 (or no valid estimate): bit-identical passthrough. */
static int drift_apply(const float* in, float* out, int n,
                       const DriftEstimate* est, int enable) {
    if (!enable || !est->valid) {
        memcpy(out, in, sizeof(float) * (size_t)n);
        return n;
    }
    double ratio = 1.0 + est->rate_ppm * 1e-6; /* A index -> B index */
    for (int i = 0; i < n; i++) {
        double src = (double)i * ratio;
        int i0 = (int)floor(src);
        double t = src - (double)i0;
        int i1 = i0 + 1 < n ? i0 + 1 : n - 1;
        out[i] = (float)((1.0 - t) * in[i0] + t * in[i1]);
    }
    return n;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* Run one duration tier: build streams, estimate, report. Fills est/offsets. */
static double run_tier(double seconds, DriftEstimate* est, double* offsets,
                       double* times, float** out_ref, float** out_b, int* out_n) {
    int n = (int)(seconds * FS);
    float* ref = malloc(sizeof(float) * (n + MAXLAG + 4));
    float* b = malloc(sizeof(float) * (n + MAXLAG + 4));
    make_streams(ref, b, n, TRUE_PPM);
    double t0 = now_ms();
    int rc = drift_estimate(ref, b, n, offsets, times, est, 1);
    double elapsed = now_ms() - t0;
    check("estimate_rc", rc == 0);
    *out_ref = ref; *out_b = b; *out_n = n;
    return elapsed;
}

int main(void) {
    int maxseg = (int)(10.0 * FS) / SEG;
    double* offsets = malloc(sizeof(double) * maxseg);
    double* times = malloc(sizeof(double) * maxseg);
    DriftEstimate est_short, est_long;
    float* ref_s; float* b_s; float* ref_l; float* b_l; int n_s, n_l;

    /* two duration tiers: 2 s and 10 s */
    double ms_short = run_tier(2.0, &est_short, offsets, times, &ref_s, &b_s, &n_s);
    double ms_long = run_tier(10.0, &est_long, offsets, times, &ref_l, &b_l, &n_l);

    /* deliverable 1: offset sequence table (long tier) */
    printf("== offset sequence (10s tier, seg=%d samples) ==\n", SEG);
    for (int k = 0; k < est_long.segments; k++)
        printf("seg=%d t=%.3f offset=%.3f\n", k, times[k], offsets[k]);

    /* offset sequence must track the monotone drift */
    int mono = 1;
    for (int k = 1; k < est_long.segments; k++)
        if (offsets[k] < offsets[k - 1] - 0.15) mono = 0;
    check("offset_monotonic", mono);

    /* deliverable 2: drift rate + cumulative sample difference, both tiers */
    double true_cum = TRUE_PPM * 1e-6 * n_l;
    printf("drift_rate_ppm_2s=%.2f\n", est_short.rate_ppm);
    printf("drift_rate_ppm_10s=%.2f\n", est_long.rate_ppm);
    printf("cumulative_samples_10s=%.2f true=%.2f\n", est_long.cum_samples, true_cum);
    check("drift_rate_2s", fabs(est_short.rate_ppm - TRUE_PPM) < 5.0);
    check("drift_rate_10s", fabs(est_long.rate_ppm - TRUE_PPM) < 2.0);
    check("cumulative_samples", fabs(est_long.cum_samples - true_cum) < 1.0);
    /* fit predicts the last measured offset: cum scaled to that segment's time */
    double t_last = times[est_long.segments - 1];
    check("cum_matches_last_offset",
          fabs(est_long.cum_samples * t_last * FS / n_l - offsets[est_long.segments - 1]) < 0.5);

    /* deliverable 3: estimation cost for the two duration tiers */
    printf("timing_2s_ms=%.3f\n", ms_short);
    printf("timing_10s_ms=%.3f\n", ms_long);

    /* estimation enabled: correction realigns B onto the A timeline */
    float* corr_out = malloc(sizeof(float) * n_l);
    drift_apply(b_l, corr_out, n_l, &est_long, 1);
    double residual = segment_offset(ref_l, corr_out, n_l, SEG);
    check("corrected_aligned", fabs(residual) < 0.5);

    /* estimation disabled: output identical to current input, no estimate */
    float* pass = malloc(sizeof(float) * n_l);
    DriftEstimate off_est;
    check("disabled_estimate_noop",
          drift_estimate(ref_l, b_l, n_l, NULL, NULL, &off_est, 0) == 0 &&
          off_est.valid == 0);
    drift_apply(b_l, pass, n_l, &est_long, 0);
    check("disabled_passthrough", memcmp(pass, b_l, sizeof(float) * n_l) == 0);

    check("bad_length", drift_estimate(ref_l, b_l, SEG, NULL, NULL, &off_est, 1) == -1);

    free(ref_s); free(b_s); free(ref_l); free(b_l);
    free(corr_out); free(pass); free(offsets); free(times);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
