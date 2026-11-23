
/* resampler: linear resampling with output-length derivation. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define N 256

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

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
