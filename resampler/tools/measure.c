
/* Measurement harness for the anti-aliasing lowpass.
 *
 * Produces the materials referenced by resampler/README.md:
 *   - doc/freq_response.csv   measured passband/stopband response (AA on/off)
 *   - stdout                  aliasing-energy table, timing table, and the
 *                             byte-identity check for the disabled path
 *
 * Built and run via `make measure` in the resampler directory.  It
 * includes src/resampler.c directly so the measurements always exercise
 * the shipped code, never a copy.
 */
#define _POSIX_C_SOURCE 199309L
#define RESAMPLER_NO_MAIN
#include "../src/resampler.c"

#include <stdint.h>
#include <time.h>

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e3 + (double)ts.tv_nsec / 1e6;
}

/* Where an input tone at f (cycles/input-sample) lands in the
 * downsampled stream (cycles/output-sample), folded into [0, 0.5].
 * Output sample j reads input position j*ratio, so the tone advances
 * ratio times faster per output sample. */
static double fold_freq(double f, double ratio) {
    double fo = f * ratio;
    fo -= floor(fo);
    if (fo > 0.5) fo = 1.0 - fo;
    return fo;
}

/* Single-bin DFT amplitude of the output at f_out, edges trimmed to
 * skip the filter's zero-padded transient. */
static double tone_amp(const float* y, int n, double f_out) {
    int a = n / 10, b = n - n / 10;
    double re = 0.0, im = 0.0;
    for (int i = a; i < b; i++) {
        double ph = 2.0 * M_PI * f_out * (double)(i - a);
        re += (double)y[i] * cos(ph);
        im -= (double)y[i] * sin(ph);
    }
    return 2.0 * sqrt(re * re + im * im) / (double)(b - a);
}

static double rms_trimmed(const float* y, int n) {
    int a = n / 10, b = n - n / 10;
    double e = 0.0;
    for (int i = a; i < b; i++) e += (double)y[i] * y[i];
    return sqrt(e / (double)(b - a));
}

/* End-to-end response: sweep an input sine across the whole input
 * band, measure the output amplitude at the folded frequency.  Above
 * the new Nyquist this directly measures the aliased residue. */
static void sweep(double ratio, FILE* fp) {
    enum { L = 8192 };
    float* x = (float*)malloc(L * sizeof(float));
    float* y = (float*)malloc(L * sizeof(float));
    if (!x || !y) { fprintf(stderr, "oom\n"); exit(1); }
    for (int k = 0; k <= 240; k++) {
        double f = 0.005 + (0.5 - 0.005) * (double)k / 240.0;
        for (int i = 0; i < L; i++) x[i] = sinf((float)(2.0 * M_PI * f * i));
        resampler_set_antialias(1);
        int n = resample_linear(x, L, ratio, y, L);
        double a_on = tone_amp(y, n, fold_freq(f, ratio));
        resampler_set_antialias(0);
        n = resample_linear(x, L, ratio, y, L);
        double a_off = tone_amp(y, n, fold_freq(f, ratio));
        fprintf(fp, "%.1f,%.5f,%.3f,%.3f\n", ratio, f,
                20.0 * log10(a_on + 1e-9), 20.0 * log10(a_off + 1e-9));
    }
    free(x);
    free(y);
}

/* Aliasing energy: input lives entirely in the stopband
 * [0.5/ratio, 0.5], so every bit of output energy is folded residue. */
static void alias_energy(double ratio) {
    enum { L = 16384, K = 24 };
    static float x[L], y[L];
    double f_lo = 0.5 / ratio + 0.004, f_hi = 0.495;
    for (int i = 0; i < L; i++) {
        double v = 0.0;
        for (int k = 0; k < K; k++) {
            double f = f_lo + (f_hi - f_lo) * (double)k / (double)(K - 1);
            v += sin(2.0 * M_PI * f * i + 0.7 * k);
        }
        x[i] = (float)(v / K);
    }
    double rms_in = rms_trimmed(x, L);
    resampler_set_antialias(0);
    int n = resample_linear(x, L, ratio, y, L);
    double rms_off = rms_trimmed(y, n);
    resampler_set_antialias(1);
    n = resample_linear(x, L, ratio, y, L);
    double rms_on = rms_trimmed(y, n);
    printf("| %.1f | %.4f | %.4f | %.6f | %.1f |\n",
           ratio, rms_in, rms_off, rms_on, 20.0 * log10(rms_off / rms_on));
}

/* Timing: 1M-sample buffer, best of 5 runs, filter on vs off. */
static void timing(double ratio, double* ms_on, double* ms_off) {
    enum { L = 1 << 20, REPS = 5 };
    float* x = (float*)malloc(L * sizeof(float));
    float* y = (float*)malloc(L * sizeof(float));
    if (!x || !y) { fprintf(stderr, "oom\n"); exit(1); }
    uint32_t s = 12345u;
    for (int i = 0; i < L; i++) {
        s = s * 1664525u + 1013904223u;
        x[i] = (float)(s >> 8) / 8388608.0f - 1.0f;
    }
    for (int aa = 1; aa >= 0; aa--) {
        resampler_set_antialias(aa);
        double best = 1e30;
        for (int r = 0; r < REPS; r++) {
            double t0 = now_ms();
            resample_linear(x, L, ratio, y, L);
            double dt = now_ms() - t0;
            if (dt < best) best = dt;
        }
        if (aa) *ms_on = best; else *ms_off = best;
    }
    resampler_set_antialias(1);
    free(x);
    free(y);
}

/* Disabled path vs the reference copy of the original algorithm. */
static int byte_identical(void) {
    enum { L = 4096 };
    static float x[L], a[L * 2], b[L * 2];
    uint32_t s = 777u;
    for (int i = 0; i < L; i++) {
        s = s * 1664525u + 1013904223u;
        x[i] = (float)(s >> 8) / 8388608.0f - 1.0f;
    }
    resampler_set_antialias(0);
    const double ratios[] = {2.0, 4.0, 1.3, 0.5};
    for (int r = 0; r < 4; r++) {
        int na = resample_linear(x, L, ratios[r], a, L * 2);
        int nb = resample_linear_ref(x, L, ratios[r], b, L * 2);
        if (na != nb || memcmp(a, b, (size_t)na * sizeof(float)) != 0) return 0;
    }
    resampler_set_antialias(1);
    return 1;
}

int main(void) {
    FILE* fp = fopen("doc/freq_response.csv", "w");
    if (!fp) {
        fprintf(stderr, "cannot write doc/freq_response.csv (run `make measure`)\n");
        return 1;
    }
    fprintf(fp, "ratio,freq,db_on,db_off\n");
    sweep(2.0, fp);
    sweep(4.0, fp);
    fclose(fp);
    printf("wrote doc/freq_response.csv\n\n");

    printf("## 混叠能量（输入全部位于阻带，输出即回卷残留）\n\n");
    printf("| 比率 | 输入 RMS | 输出 RMS（关滤波） | 输出 RMS（开滤波） | 抑制 (dB) |\n");
    printf("|---|---|---|---|---|\n");
    alias_energy(2.0);
    alias_energy(4.0);

    printf("\n## 耗时（1M 采样，5 次取最优，毫秒）\n\n");
    printf("| 比率 | 关滤波 | 开滤波 | 倍数 |\n");
    printf("|---|---|---|---|\n");
    for (int r = 0; r < 2; r++) {
        double ratio = r == 0 ? 2.0 : 4.0, on, off;
        timing(ratio, &on, &off);
        printf("| %.1f | %.2f | %.2f | %.1fx |\n", ratio, off, on, on / off);
    }

    printf("\n## 兼容性\n\n");
    printf("关闭滤波后输出与旧实现逐字节一致: %s\n",
           byte_identical() ? "PASS" : "FAIL");
    return 0;
}
