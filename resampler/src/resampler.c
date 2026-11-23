
/* resampler: linear resampling with output-length derivation.
 *
 * Anti-aliasing lowpass for downsampling (ratio > 1)
 * --------------------------------------------------
 * Plain linear interpolation followed by decimation folds every input
 * component above the new Nyquist (0.5/ratio cycles per input sample)
 * back into the kept band — audible as harsh, metallic ringing.  When
 * ratio > 1 we therefore apply a Kaiser-windowed-sinc FIR lowpass
 * before interpolating.  Design (frequencies normalized to the input
 * sample rate):
 *
 *   new Nyquist        f_nyq = 0.5 / ratio
 *   passband edge      f_p   = 0.4 / ratio   (keep 80% of the new band)
 *   stopband edge      f_s   = 0.5 / ratio   (= f_nyq)
 *   transition width   df    = f_s - f_p = 0.1 / ratio
 *   stopband attenu.   A     = 60 dB
 *   Kaiser order est.  N >= (A - 7.95)/(14.36*df) + 1 = 36.2*ratio + 1,
 *                      rounded up to odd (Type I, integer group delay)
 *   Kaiser beta (A>50) beta  = 0.1102*(A - 8.7) = 5.65
 *   sinc cutoff (-6dB) f_c   = (f_p + f_s)/2 = 0.45 / ratio
 *
 * The impulse response is normalized to unity DC gain.  Measured
 * passband/stopband curves, the aliasing-energy table and timings are
 * in README.md (regenerate with `make measure`).
 *
 * resampler_set_antialias(0) restores the exact pre-filter code path;
 * output is then byte-identical to the original implementation
 * (verified by off_byte_identical below and by `make measure`).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 256
#define AA_ATTEN_DB 60.0

static int g_antialias = 1;

void resampler_set_antialias(int enabled) { g_antialias = enabled != 0; }

/* Zeroth-order modified Bessel function via its power series. */
static double bessel_i0(double x) {
    double y = x * x / 4.0, term = 1.0, sum = 1.0;
    for (int k = 1; k < 60; k++) {
        term *= y / ((double)k * (double)k);
        sum += term;
        if (term < 1e-14 * sum) break;
    }
    return sum;
}

static double kaiser_beta(double atten_db) {
    if (atten_db > 50.0) return 0.1102 * (atten_db - 8.7);
    if (atten_db > 21.0)
        return 0.5842 * pow(atten_db - 21.0, 0.4) + 0.07886 * (atten_db - 21.0);
    return 0.0;
}

/* Odd tap count from the Kaiser order estimate. */
static int kaiser_order(double atten_db, double df) {
    int n = (int)ceil((atten_db - 7.95) / (14.36 * df)) + 1;
    if (n < 3) n = 3;
    return n | 1;
}

/* Windowed-sinc lowpass with unity DC gain; returns the tap count. */
static int design_lowpass(double fc, double df, float* h) {
    int n = kaiser_order(AA_ATTEN_DB, df);
    int m = n / 2;
    double beta = kaiser_beta(AA_ATTEN_DB);
    double denom = bessel_i0(beta);
    double sum = 0.0;
    for (int pass = 0; pass < 2; pass++) {
        for (int i = 0; i < n; i++) {
            int k = i - m;
            double s = k == 0 ? 1.0 : sin(2.0 * M_PI * fc * k) / (2.0 * M_PI * fc * k);
            double r = (double)k / (double)m;
            double w = bessel_i0(beta * sqrt(1.0 - r * r)) / denom;
            double tap = 2.0 * fc * s * w;
            if (pass == 0) sum += tap;
            else h[i] = (float)(tap / sum);
        }
    }
    return n;
}

/* Direct convolution, zero-padded edges. */
static void apply_lowpass(const float* in, int n, const float* h, int ntaps, float* out) {
    int m = ntaps / 2;
    for (int i = 0; i < n; i++) {
        double acc = 0.0;
        for (int k = 0; k < ntaps; k++) {
            int j = i + k - m;
            if (j >= 0 && j < n) acc += (double)h[k] * (double)in[j];
        }
        out[i] = (float)acc;
    }
}

static int resample_linear(const float* in, int n_in, double ratio, float* out, int cap) {
    if (n_in <= 0 || ratio <= 0.0) return -1;
    int n_out = (int)floor((double)n_in / ratio);
    if (n_out <= 0 || n_out > cap) return -1;

    const float* src = in;
    float* filtered = NULL;
    float* taps = NULL;
    if (g_antialias && ratio > 1.0) {
        int ntaps = kaiser_order(AA_ATTEN_DB, 0.1 / ratio);
        taps = (float*)malloc((size_t)ntaps * sizeof(float));
        filtered = (float*)malloc((size_t)n_in * sizeof(float));
        if (taps && filtered) {
            design_lowpass(0.45 / ratio, 0.1 / ratio, taps);
            apply_lowpass(in, n_in, taps, ntaps, filtered);
            src = filtered;
        } else {
            /* allocation failed: fall back to the unfiltered path */
            free(taps);
            free(filtered);
            taps = filtered = NULL;
        }
    }
    for (int i = 0; i < n_out; i++) {
        double spos = (double)i * ratio;
        int i0 = (int)floor(spos);
        double t = spos - (double)i0;
        int i1 = i0 + 1 < n_in ? i0 + 1 : n_in - 1;
        out[i] = (float)((1.0 - t) * src[i0] + t * src[i1]);
    }
    free(taps);
    free(filtered);
    return n_out;
}

/* Verbatim copy of the original (pre-antialias) algorithm, kept as the
 * regression anchor: with the filter disabled, resample_linear must
 * reproduce this byte for byte. */
static int resample_linear_ref(const float* in, int n_in, double ratio, float* out, int cap) {
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

#ifndef RESAMPLER_NO_MAIN

static float peak(const float* v, int n) {
    float p = 0.0f;
    for (int i = 0; i < n; i++) { float a = fabsf(v[i]); if (a > p) p = a; }
    return p;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
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

    /* anti-aliasing: a tone above the new Nyquist (0.25 at ratio 2)
     * must be suppressed when the filter is on */
    float hi[N], y_on[N], y_off[N];
    for (int i = 0; i < N; i++) hi[i] = sinf((float)(2.0 * M_PI * 0.4 * i));
    resampler_set_antialias(1);
    int n_on = resample_linear(hi, N, 2.0, y_on, N);
    resampler_set_antialias(0);
    int n_off = resample_linear(hi, N, 2.0, y_off, N);
    resampler_set_antialias(1);
    double e_on = 0.0, e_off = 0.0;  /* middle half: skip the edge transient */
    for (int i = n_on / 4; i < 3 * n_on / 4; i++) {
        e_on += (double)y_on[i] * y_on[i];
        e_off += (double)y_off[i] * y_off[i];
    }
    check("aa_suppresses_alias", n_on == n_off && e_on < e_off * 1e-4);

    /* filter disabled: byte-identical to the original algorithm */
    float ref[N * 4];
    resampler_set_antialias(0);
    int same = 1;
    for (int r = 0; r < 3; r++) {
        double ratio = r == 0 ? 2.0 : r == 1 ? 1.3 : 0.5;
        int n_a = resample_linear(in, N, ratio, out, N * 4);
        int n_b = resample_linear_ref(in, N, ratio, ref, N * 4);
        if (n_a != n_b || memcmp(out, ref, (size_t)n_a * sizeof(float)) != 0) same = 0;
    }
    resampler_set_antialias(1);
    check("off_byte_identical", same);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

#endif /* RESAMPLER_NO_MAIN */
