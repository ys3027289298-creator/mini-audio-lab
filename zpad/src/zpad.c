/* zpad: 补零插值误差衰减律标定 + 最小倍率推导 + 检测器接入验证。
 *
 * 推导两个量:
 *   量1 插值误差衰减率: eps(R, Delta) ≈ C0 · Delta^-m · R^-k
 *       (R = 补零倍率, Delta = 双音间距/bin, 误差为抛物线插值对 DTFT 真实峰位的偏差)
 *   量2 达标最小倍率:   R_min = min{ 2^j : C_use · R^-k <= EPS_TARGET }
 * 检测器以 R_min 补零; 阈值由 Dirichlet 核副瓣电平与幅度指标推导, 不拍脑袋。
 *
 * 输出五段: 1.倍率-误差曲线 2.衰减率拟合与 R_min 3.阈值推导+误判计数
 *           4.两档倍率耗时 5.衰减率失效的间距反例
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define N    256
#define RMAX 64
#define MMAX (N * RMAX)
#define PI   3.14159265358979323846

#define EPS_TARGET 0.01   /* 目标插值误差 (bin) */
#define D_MIN      4.0    /* 设计最小间距 (bin) */
#define A_MIN      0.5    /* 指标: 最小音幅 */
#define A_MAX      1.0    /* 指标: 最大音幅 */
#define F_MIN      16.0   /* 指标: 最低音频 (bin), 决定镜像泄漏距离 */
#define F_MAX      110.0
#define N_TRIAL    5000

typedef struct { double f, a, ph; } Tone;

static double complex buf[MMAX];
static double   magv[MMAX];
static unsigned char sup[MMAX];
static double complex xs[N];
static double xr[N];

static void fft(double complex *a, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { double complex t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double complex wl = cexp(-I * 2.0 * PI / len);
        for (int i = 0; i < n; i += len) {
            double complex w = 1.0;
            for (int k = 0; k < len / 2; k++) {
                double complex u = a[i + k], v = a[i + k + len / 2] * w;
                a[i + k] = u + v; a[i + k + len / 2] = u - v;
                w *= wl;
            }
        }
    }
}

static uint64_t rs = 0x9e3779b97f4a7c15ull;
static double rnd(void) { /* xorshift64* -> [0,1) */
    rs ^= rs >> 12; rs ^= rs << 25; rs ^= rs >> 27;
    return (double)(rs * 0x2545F4914F6CDD1Dull >> 11) * (1.0 / 9007199254740992.0);
}

static void gen_complex(const Tone *ts, int nt) {
    for (int n = 0; n < N; n++) {
        double complex s = 0;
        for (int j = 0; j < nt; j++)
            s += ts[j].a * cexp(I * (2.0 * PI * ts[j].f * n / N + ts[j].ph));
        xs[n] = s;
    }
}

static void gen_real(const Tone *ts, int nt) {
    for (int n = 0; n < N; n++) {
        double s = 0;
        for (int j = 0; j < nt; j++)
            s += ts[j].a * cos(2.0 * PI * ts[j].f * n / N + ts[j].ph);
        xr[n] = s;
    }
}

/* 解析 DTFT (Dirichlet 核求和) 及其一二阶导数, 用于牛顿法精化真实峰位 */
static void dtft(const Tone *ts, int nt, double nu,
                 double complex *X, double complex *Xp, double complex *Xpp) {
    *X = 0; *Xp = 0; *Xpp = 0;
    for (int j = 0; j < nt; j++) {
        double complex c = ts[j].a * cexp(I * ts[j].ph);
        double u = nu - ts[j].f;
        double complex D = 0, D1 = 0, D2 = 0;
        for (int n = 0; n < N; n++) {
            double complex w = cexp(-I * 2.0 * PI * u * n / N);
            double complex dn = -I * 2.0 * PI * n / N;
            D += w; D1 += dn * w; D2 += dn * dn * w;
        }
        *X += c * D; *Xp += c * D1; *Xpp += c * D2;
    }
}

static double refine_peak(const Tone *ts, int nt, double nu0) {
    double nu = nu0;
    for (int it = 0; it < 60; it++) {
        double complex X, Xp, Xpp;
        dtft(ts, nt, nu, &X, &Xp, &Xpp);
        double g1 = 2.0 * creal(conj(X) * Xp);
        double g2 = 2.0 * (cabs(Xp) * cabs(Xp) + creal(conj(X) * Xpp));
        double step = g1 / g2;
        nu -= step;
        if (fabs(step) < 1e-13) break;
    }
    return nu;
}

/* 补零 R 倍 -> FFT -> 窗内最大点抛物线插值, 返回估计峰位 (bin) */
static double est_freq(int R, double flo, double fhi) {
    int M = N * R;
    for (int i = 0; i < N; i++) buf[i] = xs[i];
    for (int i = N; i < M; i++) buf[i] = 0;
    fft(buf, M);
    int klo = (int)ceil(flo * R), khi = (int)(fhi * R);
    int k0 = klo;
    double best = -1;
    for (int k = klo; k <= khi; k++) {
        double m = cabs(buf[k]);
        if (m > best) { best = m; k0 = k; }
    }
    double d = 0;
    if (k0 > klo && k0 < khi) {
        double ym = cabs(buf[k0 - 1]), y0 = cabs(buf[k0]), yp = cabs(buf[k0 + 1]);
        double den = ym - 2.0 * y0 + yp;
        if (fabs(den) > 1e-30) d = 0.5 * (ym - yp) / den;
    }
    return (k0 + d) / R;
}

/* 检测器: R 倍补零, 迭代取最大 + 主峰半宽(±2 bin)抑制, 阈值 T 截断 */
static int detect(int R, double T, double *fout, int maxout) {
    int M = N * R;
    for (int i = 0; i < N; i++) buf[i] = xr[i];
    for (int i = N; i < M; i++) buf[i] = 0;
    fft(buf, M);
    for (int i = 0; i < M; i++) magv[i] = cabs(buf[i]);
    int klo = 4 * R, khi = 124 * R;
    memset(sup + klo - 1, 0, (size_t)(khi - klo + 3));
    int rad = (int)(2.0 * R + 0.5); /* 主瓣零-零宽 2 bin, 抑制其半宽 */
    int nd = 0;
    while (nd < maxout) {
        int k0 = -1;
        double best = 0;
        for (int k = klo; k <= khi; k++)
            if (!sup[k] && magv[k] > best) { best = magv[k]; k0 = k; }
        if (k0 < 0 || best <= T) break;
        double d = 0;
        if (k0 > klo && k0 < khi) {
            double ym = magv[k0 - 1], y0 = magv[k0], yp = magv[k0 + 1];
            double den = ym - 2.0 * y0 + yp;
            if (fabs(den) > 1e-30) d = 0.5 * (ym - yp) / den;
        }
        fout[nd++] = (k0 + d) / R;
        int a = k0 - rad < klo ? klo : k0 - rad;
        int b = k0 + rad > khi ? khi : k0 + rad;
        for (int k = a; k <= b; k++) sup[k] = 1;
    }
    return nd;
}

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_d(const void *a, const void *b) {
    double x = *(const double *)a - *(const double *)b;
    return (x > 0) - (x < 0);
}

int main(void) {
    int failures = 0;
    static const int Rs[] = { 1, 2, 4, 8, 16, 32, 64 };
    static const double Ds[] = { 0.0, 16.0, 8.0, 4.0 }; /* 0 = 单音 */
    static const double offs[] = { 0.13, 0.47, 0.83 };
    enum { NR = 7, ND = 4, NPH = 12 };
    double err[ND][NR];

    /* ---- 1. 倍率-误差曲线: eps = max_{相位,频偏} |f_hat(R) - nu_DTFT| ---- */
    for (int d = 0; d < ND; d++)
        for (int r = 0; r < NR; r++) err[d][r] = 0;
    for (int d = 0; d < ND; d++) {
        for (int o = 0; o < 3; o++) {
            for (int p = 0; p < NPH; p++) {
                Tone ts[2];
                int nt;
                ts[0].f = 60.0 + offs[o]; ts[0].a = (d == 0) ? 1.0 : A_MIN; ts[0].ph = 0;
                if (d == 0) { nt = 1; }
                else {
                    ts[1].f = ts[0].f + Ds[d]; ts[1].a = A_MAX;
                    ts[1].ph = 2.0 * PI * p / NPH;
                    nt = 2;
                }
                gen_complex(ts, nt);
                double nustar = refine_peak(ts, nt, ts[0].f);
                for (int r = 0; r < NR; r++) {
                    double fh = est_freq(Rs[r], ts[0].f - 1.5, ts[0].f + 1.5);
                    double e = fabs(fh - nustar);
                    if (e > err[d][r]) err[d][r] = e;
                }
            }
        }
    }
    printf("== 1. 倍率-误差曲线 (max over 12 相位 x 3 频偏, 单位 bin) ==\n");
    printf("%4s %12s %12s %12s %12s\n", "R", "single", "d=16", "d=8", "d=4");
    for (int r = 0; r < NR; r++)
        printf("%4d %12.3e %12.3e %12.3e %12.3e\n",
               Rs[r], err[0][r], err[1][r], err[2][r], err[3][r]);

    /* ---- 2. 衰减率拟合: k 来自单音, C(Delta) 拟合 Delta^-m, 推 R_min ---- */
    double sx = 0, sy = 0, sxx = 0, sxy = 0; int np = 0;
    for (int r = 1; r < NR; r++) { /* R>=2 进入渐近区 */
        double lx = log((double)Rs[r]), ly = log(err[0][r] + 1e-300);
        sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly; np++;
    }
    double k = -(np * sxy - sx * sy) / (np * sxx - sx * sx);
    double Cs = 0;
    for (int r = 0; r < NR; r++) {
        double c = err[0][r] * pow((double)Rs[r], k);
        if (c > Cs) Cs = c;
    }
    double Cd[3], Dd[3] = { 16.0, 8.0, 4.0 };
    for (int d = 1; d < ND; d++) {
        Cd[d - 1] = 0;
        for (int r = 0; r < NR; r++) {
            double c = err[d][r] * pow((double)Rs[r], k);
            if (c > Cd[d - 1]) Cd[d - 1] = c;
        }
    }
    sx = sy = sxx = sxy = 0; np = 0;
    for (int i = 0; i < 3; i++) {
        double lx = log(Dd[i]), ly = log(Cd[i]);
        sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly; np++;
    }
    double m = -(np * sxy - sx * sy) / (np * sxx - sx * sx);
    double C0 = exp((sy + m * sx) / np);
    double rho = 0; /* 律的最大残差比 */
    for (int d = 1; d < ND; d++)
        for (int r = 0; r < NR; r++) {
            double pred = C0 * pow(Ds[d], -m) * pow((double)Rs[r], -k);
            double q = err[d][r] / pred;
            if (q > rho) rho = q;
        }
    double C_two = C0 * pow(D_MIN, -m) * rho * 1.25; /* 间距外推至 D_MIN + 25% 余量 */
    double C_use = C_two > Cs * 1.25 ? C_two : Cs * 1.25;
    int R_min = 0;
    for (int R = 1; R <= RMAX; R <<= 1)
        if (C_use * pow((double)R, -k) <= EPS_TARGET) { R_min = R; break; }
    printf("\n== 2. 衰减率与最小倍率 ==\n");
    printf("k=%.3f  Cs(single)=%.4f\n", k, Cs);
    printf("C(16)=%.4f  C(8)=%.4f  C(4)=%.4f\n", Cd[0], Cd[1], Cd[2]);
    printf("m=%.3f  C0=%.4f  残差比 rho=%.3f\n", m, C0, rho);
    printf("律: eps(R,Delta) <= %.4f · Delta^(-%.3f) · R^(-%.3f)  (含余量 C_use=%.4f)\n",
           C0, m, k, C_use);
    printf("目标 eps<=%.3f bin @ Delta_min=%.1f -> R_min=%d\n", EPS_TARGET, D_MIN, R_min);
    if (R_min == 0) { printf("FATAL: R_min 超出 64\n"); failures++; R_min = RMAX; }

    /* ---- 3. 阈值推导 + 误判计数 ---- */
    double s_side = 0; /* |sinc| 在 |x|>=2 bin 的最大值 (第一可漏副瓣) */
    for (double x = 2.0; x <= 128.0; x += 1e-4) {
        double v = fabs(sin(PI * x) / (PI * x));
        if (v > s_side) s_side = v;
    }
    double s_img = 1.0 / (PI * (2.0 * F_MIN - 2.0));       /* 镜像裙边上界 */
    double d_bias = (3.0 / (PI * PI)) * (A_MAX / A_MIN) / D_MIN + EPS_TARGET;
    double s_nb = 1.0 / (PI * (D_MIN - d_bias));           /* 近邻主瓣裙边@真峰 */
    double B_rej = 2.0 * (s_side + s_img) * A_MAX;         /* 副瓣巧合上界 (N/2 单位) */
    double B_acc = A_MIN - A_MAX * s_nb - A_MIN * s_img;   /* 真峰幅度下界 */
    double Tn = sqrt(B_rej * B_acc);                       /* 几何均值阈值 */
    double T = Tn * N / 2.0;
    printf("\n== 3. 阈值推导与误判计数 ==\n");
    printf("s_side(>=2bin)=%.5f  s_img=%.5f  d_bias=%.4f  s_nb=%.5f\n",
           s_side, s_img, d_bias, s_nb);
    printf("B_rej=%.4f  B_acc=%.4f  (单位 N/2, 需 B_rej < B_acc)\n", B_rej, B_acc);
    printf("T = sqrt(B_rej·B_acc) · N/2 = %.4f · N/2 = %.2f\n", Tn, T);
    if (!(B_rej < B_acc)) { printf("FATAL: 阈值窗口不存在\n"); failures++; }

    int misjud = 0, boundviol = 0;
    double maxferr = 0;
    double fdet[8];
    for (int t = 0; t < N_TRIAL; t++) {
        Tone ts[2];
        int nt;
        if (t % 20 == 0) { nt = 0; }
        else {
            nt = (rnd() < 0.5) ? 1 : 2;
            ts[0].f = F_MIN + rnd() * (F_MAX - F_MIN);
            ts[0].a = A_MIN + rnd() * (A_MAX - A_MIN);
            ts[0].ph = rnd() * 2 * PI;
            if (nt == 2) {
                double room = ts[0].f - F_MIN;
                if (F_MAX - ts[0].f < room) room = F_MAX - ts[0].f;
                if (room < D_MIN) { nt = 1; }
                else {
                    double dlt = D_MIN + rnd() * (room < 40.0 ? room - D_MIN : 40.0 - D_MIN);
                    double dir = (rnd() < 0.5) ? 1.0 : -1.0;
                    double f2 = ts[0].f + dir * dlt;
                    if (f2 < F_MIN || f2 > F_MAX) f2 = ts[0].f - dir * dlt;
                    ts[1].f = f2;
                    ts[1].a = A_MIN + rnd() * (A_MAX - A_MIN);
                    ts[1].ph = rnd() * 2 * PI;
                }
            }
        }
        gen_real(ts, nt);
        int nd = detect(R_min, T, fdet, 8);
        if (nd != nt) { misjud++; continue; }
        if (nt > 0) {
            double ft[2];
            ft[0] = ts[0].f;
            ft[1] = (nt == 2) ? ts[1].f : 0;
            qsort(fdet, (size_t)nd, sizeof(double), cmp_d);
            qsort(ft, (size_t)nt, sizeof(double), cmp_d);
            for (int i = 0; i < nt; i++) {
                double e = fabs(fdet[i] - ft[i]);
                if (e > maxferr) maxferr = e;
                double bound = EPS_TARGET + 0.02;
                if (nt == 2) {
                    double a_this = (fabs(ft[i] - ts[0].f) < 1e-9) ? ts[0].a : ts[1].a;
                    double a_oth = (fabs(ft[i] - ts[0].f) < 1e-9) ? ts[1].a : ts[0].a;
                    double dd = fabs(ts[1].f - ts[0].f);
                    bound = EPS_TARGET + 1.5 * (3.0 / (PI * PI)) * (a_oth / a_this) / dd + 0.012;
                }
                if (e > bound) boundviol++;
            }
        }
    }
    printf("trials=%d  misjudgments=%d  max_freq_err=%.4f bin  bound_violations=%d\n",
           N_TRIAL, misjud, maxferr, boundviol);
    if (misjud != 0) { printf("FATAL: 误判不为零\n"); failures++; }
    if (boundviol != 0) { printf("FATAL: 频率误差越出推导界\n"); failures++; }

    /* ---- 4. 两档倍率耗时 ---- */
    Tone td[2] = { { 60.3, 1.0, 0.0 }, { 66.7, 0.7, 1.0 } };
    gen_real(td, 2);
    int R2 = R_min * 4 <= RMAX ? R_min * 4 : RMAX;
    int reps = 2000;
    volatile double sink = 0;
    double t0 = now_s();
    for (int i = 0; i < reps; i++) { sink += (double)detect(R_min, T, fdet, 8) + fdet[0]; }
    double t1 = now_s();
    for (int i = 0; i < reps; i++) { sink += (double)detect(R2, T, fdet, 8) + fdet[0]; }
    double t2 = now_s();
    double us1 = (t1 - t0) / reps * 1e6, us2 = (t2 - t1) / reps * 1e6;
    printf("\n== 4. 两档倍率耗时 (检测全流程, %d 次平均) ==\n", reps);
    printf("R=%-3d %8.2f us/frame\nR=%-3d %8.2f us/frame  (x%.2f)\n",
           R_min, us1, R2, us2, us2 / us1);

    /* ---- 5. 反例: 间距低于主瓣宽, 衰减率失效 ---- */
    printf("\n== 5. 衰减率失效反例: Delta=0.75 bin (< 主瓣零-零宽 2 bin) ==\n");
    printf("%4s %14s %14s\n", "R", "eps_true(bin)", "律外推(bin)");
    static const int Rc[] = { 4, 8, 16, 32, 64 };
    for (int r = 0; r < 5; r++) {
        double emax = 0;
        for (int o = 0; o < 3; o++)
            for (int p = 0; p < NPH; p++) {
                Tone tc[2] = {
                    { 60.0 + offs[o], 1.0, 0.0 },
                    { 60.0 + offs[o] + 0.75, 1.0, 2.0 * PI * p / NPH }
                };
                gen_complex(tc, 2);
                double fh = est_freq(Rc[r], tc[0].f - 1.5, tc[0].f + 1.5);
                double e = fabs(fh - tc[0].f);
                if (e > emax) emax = e;
            }
        double pred = C0 * pow(0.75, -m) * pow((double)Rc[r], -k);
        printf("%4d %14.3e %14.3e\n", Rc[r], emax, pred);
    }
    Tone tm[2] = { { 60.13, 1.0, 0.0 }, { 60.88, 1.0, 0.6 } };
    gen_real(tm, 2);
    int ndm = detect(R_min, T, fdet, 8);
    printf("检测器在 Delta=0.75 上报告 %d 个音 (真实 2 个): 峰合并, 间距指标 D_MIN=%.1f 是衰减律的适用下界\n",
           ndm, D_MIN);

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
