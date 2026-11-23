
/* bench: 去趋势阶数 p、高通截止 fc 与低频残留的耦合实验。
   输出:控制台表格 + ../artifacts/detrend_*.csv */
#define _POSIX_C_SOURCE 199309L
#include "detrend.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FS    49152.0   /* 使 K*FB/FS = 40 恰为整数 bin */
#define N     16384     /* 分析段长 */
#define K     8192      /* 能量评估窗(取中段,避开边界瞬态) */
#define FB    240.0     /* 低频评估带 [0, FB) */
#define KB    40        /* 带内单边 bin 数 = K*FB/FS */
#define NREP  64        /* 噪声实现平均次数 */
#define PMAX  6
#define NFC   9

static const double FC[NFC] = {5, 10, 20, 40, 80, 120, 160, 240, 360};

static unsigned long long rng_s = 0x9e3779b97f4a7c15ull;
static double frand(void) {
    rng_s ^= rng_s >> 12; rng_s ^= rng_s << 25; rng_s ^= rng_s >> 27;
    return (rng_s * 0x2545F4914F6CDD1Dull >> 11) * (1.0 / 9007199254740992.0);
}
static double grand(void) { /* Box-Muller */
    double u = frand() + 1e-300, v = frand();
    return sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
}

/* 趋势 = 3 次多项式 + 非多项式慢漂(半波正弦,~1.5 Hz) */
static double trend(double x) {
    return 20.0 * (x * x * x - 0.5 * x * x + 0.25 * x - 0.1)
         + 3.0 * sin(M_PI * (x + 1.0));
}

static double trans_of(double fc) { double t = 0.5 * fc; return t < 15.0 ? 15.0 : t; }

static double hann[K];
static double W2;

/* 取中段 K 点加 Hann 窗做 FFT */
static void eval_spectrum(const double* sig, double* re, double* im) {
    int off = (N - K) / 2;
    for (int i = 0; i < K; i++) { re[i] = sig[off + i] * hann[i]; im[i] = 0.0; }
    fft_double(re, im, K, 0);
}

/* 带内均方(每样本):E = 2/(K*W2) * sum_{k0..k1} S(k) G(k),0 号 bin 半权重 */
static double band_ms_range(const double* re, const double* im, const double* G,
                            int k0, int k1) {
    double s = 0.0;
    for (int k = k0; k <= k1; k++) {
        double wgt = (k == 0) ? 0.5 : 1.0;
        s += wgt * (re[k] * re[k] + im[k] * im[k]) * (G ? G[k] : 1.0);
    }
    return 2.0 * s / (K * W2);
}
static double band_ms(const double* re, const double* im, const double* G) {
    return band_ms_range(re, im, G, 0, KB);
}

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int main(void) {
    static double x[N], w[N], t[N], r[N], y[N];
    static double re[K], im[K];
    static double G[NFC][KB + 1];       /* 各截止的功率响应(评估 bin 上) */
    static double hbuf[16384];
    static double gre[32768], gim[32768];

    for (int i = 0; i < K; i++) {
        hann[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / (K - 1));
        W2 += hann[i] * hann[i];
    }
    for (int i = 0; i < N; i++)
        t[i] = trend(-1.0 + 2.0 * i / (N - 1));

    /* 各 fc 的 FIR 功率响应:32768 点 FFT,bin 4k 与评估 bin k 同频 */
    for (int f = 0; f < NFC; f++) {
        int m = fir_highpass_design(FC[f], FS, trans_of(FC[f]), hbuf, 16384);
        for (int i = 0; i < 32768; i++) { gre[i] = i < m ? hbuf[i] : 0.0; gim[i] = 0.0; }
        fft_double(gre, gim, 32768, 0);
        for (int k = 0; k <= KB; k++)
            G[f][k] = gre[4 * k] * gre[4 * k] + gim[4 * k] * gim[4 * k];
    }

    static double E_d[PMAX + 1], E_casc[PMAX + 1][NFC];
    static double E_nd[PMAX + 1][NFC], E_nd_d[PMAX + 1];
    static double E_tr[PMAX + 1];
    double E_w = 0.0, E_wsub = 0.0, E_ndsub = 0.0;
    static double E_wf[NFC];
    double norm_prev = 1e300;
    int norm_mono = 1;

    /* 趋势 leftover(确定性,单次) */
    for (int p = 0; p <= PMAX; p++) {
        detrend_poly(t, N, p, r);
        eval_spectrum(r, re, im);
        E_tr[p] = band_ms(re, im, NULL);
        double nrm = 0.0;
        for (int i = 0; i < N; i++) nrm += r[i] * r[i];
        if (nrm > norm_prev * 1.0001) norm_mono = 0;
        norm_prev = nrm;
    }
    printf("sanity_leftover_norm_nonincreasing=%d\n", norm_mono);

    for (int rep = 0; rep < NREP; rep++) {
        for (int i = 0; i < N; i++) { w[i] = grand(); x[i] = t[i] + w[i]; }

        eval_spectrum(w, re, im);              /* 纯噪声基准 */
        E_w += band_ms(re, im, NULL);
        E_wsub += band_ms_range(re, im, NULL, 1, 4);   /* 最低 4 bin(0-24 Hz) */
        for (int f = 0; f < NFC; f++) E_wf[f] += band_ms(re, im, G[f]);

        for (int p = 0; p <= PMAX; p++) {
            detrend_poly(x, N, p, r);          /* 信号+噪声:去趋势 */
            eval_spectrum(r, re, im);
            E_d[p] += band_ms(re, im, NULL);
            for (int f = 0; f < NFC; f++) E_casc[p][f] += band_ms(re, im, G[f]);

            detrend_poly(w, N, p, r);          /* 纯噪声:去趋势(底噪) */
            eval_spectrum(r, re, im);
            E_nd_d[p] += band_ms(re, im, NULL);
            if (p == PMAX) E_ndsub += band_ms_range(re, im, NULL, 1, 4);
            for (int f = 0; f < NFC; f++) E_nd[p][f] += band_ms(re, im, G[f]);
        }
    }
    E_wsub /= NREP; E_ndsub /= NREP;
    for (int p = 0; p <= PMAX; p++) {
        E_d[p] /= NREP; E_nd_d[p] /= NREP;
        for (int f = 0; f < NFC; f++) { E_casc[p][f] /= NREP; E_nd[p][f] /= NREP; }
    }
    E_w /= NREP;
    for (int f = 0; f < NFC; f++) E_wf[f] /= NREP;

    const double D_B = 2.0 * KB;               /* 带内噪声自由度 */
    const double NOISE_MS = 2.0 * FB / FS;     /* 白噪声带内均方(理论) */

    /* ---- 1. 各阶残留表 ---- */
    printf("== residual table (band [0,%g Hz], ms/sample, avg of %d reps) ==\n", FB, NREP);
    printf("p, E_meas, E_trend_left, E_noise_meas, E_noise_model, E_theory, meas/theory\n");
    FILE* fp = fopen("../artifacts/detrend_residual.csv", "w");
    fprintf(fp, "p,E_meas,E_trend_left,E_noise_meas,E_noise_model,E_theory,meas_over_theory\n");
    for (int p = 0; p <= PMAX; p++) {
        double noise_model = NOISE_MS * (1.0 - (p + 1) / D_B);
        double eth = E_tr[p] + E_nd_d[p];
        printf("%d, %.6e, %.6e, %.6e, %.6e, %.6e, %.3f\n",
               p, E_d[p], E_tr[p], E_nd_d[p], noise_model, eth, E_d[p] / eth);
        fprintf(fp, "%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.6f\n",
                p, E_d[p], E_tr[p], E_nd_d[p], noise_model, eth, E_d[p] / eth);
    }
    fclose(fp);
    printf("noise band-ms: raw %.6e, detrended(p=%d) %.6e; sub-band[6-24Hz]: raw %.6e, detrended %.6e\n",
           E_w, PMAX, E_nd_d[PMAX], E_wsub, E_ndsub);

    /* ---- 选阶:首个进入平台(<= 1.5x 最高阶残留)的 p ---- */
    int pstar = PMAX;
    for (int p = 0; p <= PMAX; p++)
        if (E_d[p] <= 1.5 * E_d[PMAX]) { pstar = p; break; }
    printf("p* = %d (first within 1.5x of p=%d plateau)\n", pstar, PMAX);

    /* ---- 2. 耦合项:E_meas vs 可分离预测 E_d(p)*G(fc) ---- */
    printf("\n== coupling C(p,fc) = E_meas / (E_d(p)*G(fc)), dB ==\n");
    printf("p\\fc");
    for (int f = 0; f < NFC; f++) printf(", %g", FC[f]);
    printf("\n");
    fp = fopen("../artifacts/detrend_coupling.csv", "w");
    fprintf(fp, "p,fc,E_meas,G_fc,E_sep,coupling_dB\n");
    for (int p = 0; p <= PMAX; p++) {
        printf("%d", p);
        for (int f = 0; f < NFC; f++) {
            double gfc = E_wf[f] / E_w;                 /* 实测带内能量保持率 */
            double esep = E_d[p] * gfc;
            double cdb = 10.0 * log10(E_casc[p][f] / esep);
            printf(", %+.2f", cdb);
            fprintf(fp, "%d,%g,%.8e,%.8e,%.8e,%.4f\n",
                    p, FC[f], E_casc[p][f], gfc, esep, cdb);
        }
        printf("\n");
    }
    fclose(fp);

    /* ---- 选截止:趋势 leftover 没入噪声(串联残留 <= 1.5x 串联底噪) ---- */
    int fstar = NFC - 1;
    for (int f = 0; f < NFC; f++)
        if (E_casc[pstar][f] <= 1.5 * E_nd[pstar][f]) { fstar = f; break; }
    printf("fc* = %g Hz (smallest fc with cascade residual <= 1.5x cascade floor)\n",
           FC[fstar]);

    /* ---- 3. 串联底噪 vs 单模块 ---- */
    printf("\n== noise floor: cascade(detrend p=%d + HP) vs HP alone ==\n", pstar);
    printf("fc, floor_single, floor_cascade, floor_model, casc/single dB\n");
    fp = fopen("../artifacts/detrend_floor.csv", "w");
    fprintf(fp, "fc,floor_single,floor_cascade,floor_model,casc_over_single_dB\n");
    for (int f = 0; f < NFC; f++) {
        double model = E_wf[f] * (1.0 - (pstar + 1) / D_B);
        double ratio_db = 10.0 * log10(E_nd[pstar][f] / E_wf[f]);
        printf("%g, %.6e, %.6e, %.6e, %+.3f\n",
               FC[f], E_wf[f], E_nd[pstar][f], model, ratio_db);
        fprintf(fp, "%g,%.8e,%.8e,%.8e,%.4f\n",
                FC[f], E_wf[f], E_nd[pstar][f], model, ratio_db);
    }
    fclose(fp);

    /* ---- 4. 两档截止耗时 ---- */
    printf("\n== timing: detrend(p=%d) + FIR HP, N=%d ==\n", pstar, N);
    fp = fopen("../artifacts/detrend_timing.csv", "w");
    fprintf(fp, "fc,taps,reps,ms_per_run\n");
    printf("fc, taps, reps, ms_per_run\n");
    double tfc[2] = {40.0, 240.0};
    for (int c = 0; c < 2; c++) {
        int m = fir_highpass_design(tfc[c], FS, trans_of(tfc[c]), hbuf, 16384);
        int reps = (int)(1.5e9 / ((double)N * m));
        if (reps < 5) reps = 5;
        if (reps > 2000) reps = 2000;
        double t0 = now_s();
        for (int rep = 0; rep < reps; rep++) {
            detrend_poly(x, N, pstar, r);
            fir_apply_centered(r, N, hbuf, m, y);
        }
        double ms = 1e3 * (now_s() - t0) / reps;
        printf("%g, %d, %d, %.3f\n", tfc[c], m, reps, ms);
        fprintf(fp, "%g,%d,%d,%.4f\n", tfc[c], m, reps, ms);
    }
    fclose(fp);
    return 0;
}
