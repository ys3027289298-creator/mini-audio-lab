/* peakpick: windowed-spectrum two-tone peak picking self-test.
 *
 * 复现 Python 分析引擎（tools/common.py）的确定性试验，验证三条结论：
 *   1) 旁瓣误检随窗型变化：-40dB 报告门限、ROI±8bin 内单音假峰数
 *      rect=14, hann=2, hamming=0, blackman=0；
 *   2) 两相邻音合并：blackman L=256 时 df=300Hz(1.60bin) 只剩 1 个局部极大
 *      （漏检），df=400Hz(2.13bin) 分成 2 个峰但带旁瓣假峰；
 *   3) 加大窗长后同时消除漏检与误检：blackman L=2048、df=400Hz 全 clean。
 * 窗公式必须与 tools/common.py 逐式一致（对称型，分母 L-1）。
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FS 48000.0
#define GATE_DB 40.0
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double re, im; } C;

static void fft(C* a, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = -2.0 * M_PI / len;
        double wr = cos(ang), wi = sin(ang);
        for (int i = 0; i < n; i += len) {
            double cwr = 1.0, cwi = 0.0;
            for (int k = 0; k < len / 2; k++) {
                C u = a[i + k];
                C v;
                C* q = &a[i + k + len / 2];
                v.re = q->re * cwr - q->im * cwi;
                v.im = q->re * cwi + q->im * cwr;
                a[i + k].re = u.re + v.re;
                a[i + k].im = u.im + v.im;
                q->re = u.re - v.re;
                q->im = u.im - v.im;
                double nr = cwr * wr - cwi * wi;
                cwi = cwr * wi + cwi * wr;
                cwr = nr;
            }
        }
    }
}

enum { WIN_RECT, WIN_HANN, WIN_HAMMING, WIN_BLACKMAN };

static double winval(int type, int n, int L) {
    double a = 2.0 * M_PI * n / (L - 1);
    switch (type) {
        case WIN_RECT:    return 1.0;
        case WIN_HANN:    return 0.50 - 0.50 * cos(a);
        case WIN_HAMMING: return 0.54 - 0.46 * cos(a);
        default:          return 0.42 - 0.50 * cos(a) + 0.08 * cos(2.0 * a);
    }
}

/* x=cos(2*pi f0 n/fs)+cos(2*pi f1 n/fs + phi)，加窗后零填充到 nfft，
 * 输出前 nfft/2+1 个幅值谱。f1==f0 时等价单音。*/
static void spectrum(int type, int L, double f0, double f1, double phi,
                     int nfft, double* mag) {
    C* a = (C*)calloc((size_t)nfft, sizeof(C));
    for (int n = 0; n < L; n++) {
        double v = cos(2.0 * M_PI * f0 * n / FS)
                 + cos(2.0 * M_PI * f1 * n / FS + phi);
        a[n].re = v * winval(type, n, L);
    }
    fft(a, nfft);
    mag[0] = fabs(a[0].re);
    for (int k = 1; k <= nfft / 2; k++)
        mag[k] = hypot(a[k].re, a[k].im);
    free(a);
}

/* 与 Python pick_peaks 同规格：严格局部极大且 >= 主峰*10^(-GATE/20)。
 * 返回峰个数，峰下标写入 idx（容量 cap）。*/
static int pick(const double* mag, int nbin, int* idx, int cap) {
    double mx = 0.0;
    for (int k = 0; k < nbin; k++) if (mag[k] > mx) mx = mag[k];
    double floor_ = mx * pow(10.0, -GATE_DB / 20.0);
    int cnt = 0;
    for (int k = 1; k < nbin - 1; k++) {
        if (mag[k] > mag[k - 1] && mag[k] > mag[k + 1] && mag[k] >= floor_) {
            if (cnt < cap) idx[cnt] = k;
            cnt++;
        }
    }
    return cnt;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    const char* names[4] = {"rect", "hann", "hamming", "blackman"};
    /* Python 引擎实测（L=512, NFFT=2048, ROI±8bin, -40dB 门限）：14,2,0,0 */
    int expect_fp[4] = {14, 2, 0, 0};

    /* ---- 试验 1：单音 3000Hz、L=512，ROI 内旁瓣假峰数 ----
     * 与 Python 同口径：ROI±8 个分析窗 bin，配对容差 ±0.5 分析窗 bin。
     * raw_* 行打印原始计数，供 tools/crosscheck.py 与 Python 引擎逐项比对。*/
    enum { N1 = 2048 };
    static double m1[N1 / 2 + 1];
    static int i1[256];
    double df1 = FS / N1;                       /* FFT bin 宽 23.4375Hz */
    double abin1 = FS / 512.0;                  /* 分析窗 bin 宽 93.75Hz */
    for (int t = 0; t < 4; t++) {
        spectrum(t, 512, 3000.0, 3000.0, 0.0, N1, m1);
        int cnt = pick(m1, N1 / 2 + 1, i1, 256);
        int fp = 0;
        for (int i = 0; i < cnt; i++) {
            double hz = i1[i] * df1;
            if (fabs(hz - 3000.0) > 8.0 * abin1) continue;
            if (fabs(hz - 3000.0) > 0.5 * abin1) fp++;
        }
        char label[64];
        snprintf(label, sizeof label, "sidelobe_fp_%s", names[t]);
        printf("raw_%s=%d\n", label, fp);
        check(label, fp == expect_fp[t]);
    }

    /* ---- 试验 2：blackman L=256，双音合并（df=300Hz, 1.60 bin）---- */
    enum { N2 = 1024 };
    static double m2[N2 / 2 + 1];
    static int i2[256];
    spectrum(WIN_BLACKMAN, 256, 3000.0, 3300.0, 0.0, N2, m2);
    int cnt_merged = pick(m2, N2 / 2 + 1, i2, 256);
    printf("raw_blackman_L256_df300_peakcount=%d\n", cnt_merged);
    check("blackman_L256_df300_merged_single_peak", cnt_merged == 1);

    /* ---- 试验 3：同窗长 df=400Hz（2.13 bin）分成两峰，但有旁瓣假峰 ---- */
    spectrum(WIN_BLACKMAN, 256, 3000.0, 3400.0, 0.0, N2, m2);
    int cnt_split = pick(m2, N2 / 2 + 1, i2, 256);
    printf("raw_blackman_L256_df400_peakcount=%d\n", cnt_split);
    check("blackman_L256_df400_two_peaks", cnt_split >= 2);

    /* ---- 试验 4：blackman L=2048、df=400Hz，ROI 内零误检且两峰都在 ---- */
    enum { N4 = 8192 };
    static double m4[N4 / 2 + 1];
    static int i4[256];
    spectrum(WIN_BLACKMAN, 2048, 3000.0, 3400.0, 0.0, N4, m4);
    int cnt4 = pick(m4, N4 / 2 + 1, i4, 256);
    double df4 = FS / N4;
    double bin4 = FS / 2048.0;
    int b0 = (int)lround(3000.0 / df4), b1 = (int)lround(3400.0 / df4);
    int hit0 = 0, hit1 = 0, fp4 = 0, cnt4_in_roi = 0;
    for (int i = 0; i < cnt4; i++) {
        double hz = i4[i] * df4;
        if (hz < 3000.0 - 8 * bin4 || hz > 3400.0 + 8 * bin4) continue;
        cnt4_in_roi++;
        if (abs(i4[i] - b0) <= 0) hit0 = 1;       /* 峰落在真实音最近 bin */
        else if (abs(i4[i] - b1) <= 0) hit1 = 1;
        else fp4++;
    }
    check("blackman_L2048_two_true_peaks", hit0 && hit1);
    check("blackman_L2048_zero_fp", fp4 == 0);
    printf("raw_blackman_L2048_peakcount=%d\n", cnt4_in_roi);
    printf("raw_blackman_L2048_fp=%d\n", fp4);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
