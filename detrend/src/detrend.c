
/* detrend: polynomial detrend via orthonormal discrete Legendre basis,
   windowed-sinc FIR high-pass, double-precision radix-2 FFT. */
#include "detrend.h"
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 在均匀网格 x_i = -1 + 2i/(n-1) 上构造 0..p 阶离散 Legendre 基,
   三项递推生成后做两遍修正 Gram-Schmidt 精确正交归一。 */
static void build_basis(double* B, int n, int p) {
    for (int k = 0; k <= p; k++) {
        double* row = B + (size_t)k * n;
        for (int i = 0; i < n; i++) {
            double x = -1.0 + 2.0 * i / (n - 1);
            if (k == 0) row[i] = 1.0;
            else if (k == 1) row[i] = x;
            else {
                const double* p1 = B + (size_t)(k - 1) * n;
                const double* p2 = B + (size_t)(k - 2) * n;
                row[i] = ((2.0 * k - 1.0) * x * p1[i] - (k - 1.0) * p2[i]) / k;
            }
        }
    }
    for (int pass = 0; pass < 2; pass++) {
        for (int k = 0; k <= p; k++) {
            double* row = B + (size_t)k * n;
            for (int j = 0; j < k; j++) {
                const double* prev = B + (size_t)j * n;
                double c = 0.0;
                for (int i = 0; i < n; i++) c += row[i] * prev[i];
                for (int i = 0; i < n; i++) row[i] -= c * prev[i];
            }
            double nrm = 0.0;
            for (int i = 0; i < n; i++) nrm += row[i] * row[i];
            nrm = sqrt(nrm);
            for (int i = 0; i < n; i++) row[i] /= nrm;
        }
    }
}

void detrend_poly(const double* x, int n, int p, double* out) {
    double* B = malloc((size_t)(p + 1) * n * sizeof(double));
    build_basis(B, n, p);
    for (int i = 0; i < n; i++) out[i] = x[i];
    for (int k = 0; k <= p; k++) {
        const double* row = B + (size_t)k * n;
        double c = 0.0;
        for (int i = 0; i < n; i++) c += x[i] * row[i];
        for (int i = 0; i < n; i++) out[i] -= c * row[i];
    }
    free(B);
}

int fir_highpass_design(double fc, double fs, double trans, double* h, int maxtaps) {
    int m = (int)(3.3 * fs / trans) | 1;   /* Hamming 过渡带 ~3.3 fs/M,强制奇数 */
    if (m > maxtaps) m = maxtaps | 1;
    if (m < 3) m = 3;
    int c = (m - 1) / 2;
    double fn = fc / fs;
    for (int i = 0; i < m; i++) {
        int k = i - c;
        double lp = (k == 0) ? 2.0 * fn : sin(2.0 * M_PI * fn * k) / (M_PI * k);
        double w = 0.54 - 0.46 * cos(2.0 * M_PI * i / (m - 1));
        h[i] = (k == 0 ? 1.0 : 0.0) - lp * w;   /* 谱反演:HP = delta - LP */
    }
    return m;
}

void fir_apply_centered(const double* x, int n, const double* h, int m, double* y) {
    int c = (m - 1) / 2;
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int k = 0; k < m; k++) {
            int j = i + k - c;
            if (j >= 0 && j < n) s += h[k] * x[j];
        }
        y[i] = s;
    }
}

void fft_double(double* re, double* im, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double t;
            t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * M_PI / len * (inverse ? 1.0 : -1.0);
        double wr = cos(ang), wi = sin(ang);
        for (int i = 0; i < n; i += len) {
            double cwr = 1.0, cwi = 0.0;
            for (int k = 0; k < len / 2; k++) {
                double ur = re[i + k], ui = im[i + k];
                double vr = re[i + k + len / 2] * cwr - im[i + k + len / 2] * cwi;
                double vi = re[i + k + len / 2] * cwi + im[i + k + len / 2] * cwr;
                re[i + k] = ur + vr;             im[i + k] = ui + vi;
                re[i + k + len / 2] = ur - vr;   im[i + k + len / 2] = ui - vi;
                double nwr = cwr * wr - cwi * wi;
                double nwi = cwr * wi + cwi * wr;
                cwr = nwr; cwi = nwi;
            }
        }
    }
    if (inverse) {
        for (int i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
    }
}

#ifdef DETREND_SELFTEST
#include <stdio.h>

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    enum { N = 4096 };
    double* x = malloc(N * sizeof(double));
    double* r = malloc(N * sizeof(double));

    /* 3 次多项式 + 小扰动,p=3 应精确去除多项式部分 */
    for (int i = 0; i < N; i++) {
        double t = -1.0 + 2.0 * i / (N - 1);
        x[i] = 5.0 * t * t * t - 2.0 * t * t + t - 0.3;
    }
    detrend_poly(x, N, 3, r);
    double maxr = 0.0;
    for (int i = 0; i < N; i++) if (fabs(r[i]) > maxr) maxr = fabs(r[i]);
    check("poly3_exact", maxr < 1e-8);

    /* p=0 等价于去均值 */
    for (int i = 0; i < N; i++) x[i] = 7.0 + 0.001 * (i % 3);
    detrend_poly(x, N, 0, r);
    double mean = 0.0;
    for (int i = 0; i < N; i++) mean += r[i];
    check("order0_demean", fabs(mean / N) < 1e-9);

    /* 高通:DC 衰减 > 60 dB,4 倍截止处正弦通过 */
    double h[2048];
    int m = fir_highpass_design(100.0, 8000.0, 50.0, h, 2048);
    for (int i = 0; i < N; i++) x[i] = 1.0;
    fir_apply_centered(x, N, h, m, r);
    double dc = 0.0;
    for (int i = N / 4; i < 3 * N / 4; i++) if (fabs(r[i]) > dc) dc = fabs(r[i]);
    check("hp_dc_60db", dc < 1e-3);
    for (int i = 0; i < N; i++) x[i] = sin(2.0 * M_PI * 400.0 * i / 8000.0);
    fir_apply_centered(x, N, h, m, r);
    double amp = 0.0;
    for (int i = N / 4; i < 3 * N / 4; i++) if (fabs(r[i]) > amp) amp = fabs(r[i]);
    check("hp_passband", fabs(amp - 1.0) < 0.02);

    /* FFT 往返 */
    double re[64], im[64], orig[64];
    for (int i = 0; i < 64; i++) { re[i] = (double)(i % 7) - 3.0; im[i] = 0.0; orig[i] = re[i]; }
    fft_double(re, im, 64, 0);
    fft_double(re, im, 64, 1);
    double maxd = 0.0;
    for (int i = 0; i < 64; i++) if (fabs(re[i] - orig[i]) > maxd) maxd = fabs(re[i] - orig[i]);
    check("fft_roundtrip", maxd < 1e-9);

    free(x); free(r);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
