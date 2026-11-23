/* early_reflection: pre-delayed early reflection taps with fractional delays.
 *
 * 抽头时刻 <-> 距离:
 *   t_i    = PRE_MS/1000 + dist_i / SOS              [s]
 *   D_i    = t_i * fs                                [samples]
 *   n_i    = floor(D_i),  mu_i = D_i - n_i           [整数 / 分数部分]
 *   dist_i = ((n_i + mu_i)/fs - PRE_MS/1000) * SOS   [m]
 *
 * 分数延迟 (Lagrange, 奇数阶 N): n = (N+1)/2, 抽头 j = -(n-1)..n
 *   h_j(mu) = prod_{i != j} (mu - i)/(j - i)
 *   y[k]    = sum_j h_j * x[k - n_i - j]              (延迟 D_i 的近似)
 *
 * 两个恒等式:
 *   sum_j h_j = 1                (直流增益为 1)
 *   sum_j j*h_j = mu             (一阶矩)
 *   所以抽头质心 = n_i + mu = D_i, 与解析延迟零偏差.
 *
 * 插值误差上界 (输入频谱限制在 fmax 以内, 幅度不超过 amp):
 *   |E(e^jw)| = w^(N+1)/(N+1)! * prod_j |mu - j|,  w = 2*pi*fmax/fs
 *   err <= amp * w^(N+1)/(N+1)! * prod_j |mu - j|
 *   prod_j |mu-j| 在 mu = 1/2 处取最大:
 *     N=1: 0.25, N=3: 0.5625, N=5: 3.515625
 *   相位延迟偏差 |tau_p - mu| <= err_bound / w.
 *
 * 延迟 D_i <= 0 (反射和直达声挤在同一时刻) 以及反射抽头落在偏移 0 上,
 * 均报无效.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define FS        48000.0   /* sample rate, Hz */
#define SOS       343.0     /* speed of sound, m/s */
#define PRE_MS    5.0       /* pre-delay, ms */
#define MAX_REFL  4
#define MAX_ORDER 5
#define MAX_TAPS  ((MAX_ORDER + 1) * MAX_REFL)
#define PI        3.14159265358979323846

typedef struct { int delay; double frac; } Tap;      /* D = delay + frac */

typedef struct {
    int n;
    int order;
    int pre_samples;
    Tap tap[MAX_REFL];
} TapTable;

typedef struct {
    int order;
    int n_taps;
    int first_offset;
    int offset[MAX_TAPS];
    float weight[MAX_TAPS];
} Ir;

static int valid_order(int order) {
    return order >= 1 && order <= MAX_ORDER && (order % 2) == 1;
}

/* Fractional-delay Lagrange coefficients, taps j = -(n-1)..n. */
static int frac_taps(int order, double mu, int* js, double* h) {
    if (!valid_order(order) || js == NULL || h == NULL) return -1;
    if (!(mu >= 0.0 && mu < 1.0)) return -1;
    int n = (order + 1) / 2;
    int lo = -(n - 1), hi = n;
    for (int j = lo; j <= hi; j++) {
        double w = 1.0;
        for (int i = lo; i <= hi; i++)
            if (i != j) w *= (mu - (double)i) / (double)(j - i);
        js[j - lo] = j;
        h[j - lo] = w;
    }
    return order + 1;
}

static double factorial(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; i++) f *= (double)i;
    return f;
}

/* max over mu of prod_j |mu - j|, attained at mu = 1/2 */
static double lagrange_pmax(int order) {
    if (!valid_order(order)) return -1.0;
    int n = (order + 1) / 2;
    double p = 1.0;
    for (int j = -(n - 1); j <= n; j++) p *= fabs(0.5 - (double)j);
    return p;
}

/* Worst-case interpolation error for amplitude amp within fmax. */
static double frac_error_bound(int order, double fmax, double amp) {
    if (!valid_order(order) || fmax < 0.0 || fmax > FS * 0.5) return -1.0;
    if (!(amp >= 0.0)) return -1.0;
    double w = 2.0 * PI * fmax / FS;
    return amp * pow(w, (double)(order + 1)) / factorial(order + 1)
           * lagrange_pmax(order);
}

/* distance -> tap time (integer + fractional sample delay) */
static int tap_from_distance(double dist, double pre_ms, double fs, int order, Tap* t) {
    if (t == NULL || !isfinite(dist) || dist < 0.0) return -1;
    if (!isfinite(pre_ms) || pre_ms < 0.0 || !isfinite(fs) || fs <= 0.0) return -1;
    if (!valid_order(order)) return -1;
    double D = (pre_ms * 0.001 + dist / SOS) * fs;
    if (!(D > 0.0)) return -1;                  /* 延迟为零: 反射与直达声重合 */
    double fl = floor(D);
    if (fl > 2.0e9) return -1;
    t->delay = (int)fl;
    t->frac = D - fl;
    return 0;
}

/* tap time -> distance (inverse of tap_from_distance) */
static double distance_from_tap(const Tap* t, double pre_ms, double fs) {
    if (t == NULL || !isfinite(pre_ms) || !isfinite(fs) || fs <= 0.0) return -1.0;
    return (((double)t->delay + t->frac) / fs - pre_ms * 0.001) * SOS;
}

static int tap_table_build(const float* dist, const float* gain, int nref, int order,
                           double pre_ms, double fs, TapTable* tt) {
    if (dist == NULL || gain == NULL || tt == NULL) return -1;
    if (nref <= 0 || nref > MAX_REFL || !valid_order(order)) return -1;
    if (!isfinite(pre_ms) || pre_ms < 0.0 || !isfinite(fs) || fs <= 0.0) return -1;
    memset(tt, 0, sizeof(*tt));
    tt->order = order;
    tt->pre_samples = (int)floor(pre_ms * 0.001 * fs + 0.5);
    int half = (order + 1) / 2;
    for (int r = 0; r < nref; r++) {
        if (gain[r] < 0.0f || !isfinite(gain[r])) return -1;
        if (tap_from_distance((double)dist[r], pre_ms, fs, order, &tt->tap[r]) != 0)
            return -1;
        /* earliest fractional tap must sit after the direct sound */
        if (tt->tap[r].delay - half + 1 <= 0) return -1;
    }
    tt->n = nref;
    return 0;
}

/* Expand each reflection into order+1 weighted taps, sorted by offset. */
static int ir_build(const TapTable* tt, const float* gain, Ir* ir) {
    if (tt == NULL || gain == NULL || ir == NULL) return -1;
    if (tt->n <= 0 || tt->n > MAX_REFL || !valid_order(tt->order)) return -1;
    int per = tt->order + 1;
    if (per * tt->n > MAX_TAPS) return -1;

    memset(ir, 0, sizeof(*ir));
    ir->order = tt->order;
    int k = 0, first = 0x7fffffff;
    for (int r = 0; r < tt->n; r++) {
        if (gain[r] < 0.0f || !isfinite(gain[r])) return -1;
        int js[MAX_ORDER + 1];
        double h[MAX_ORDER + 1];
        if (frac_taps(tt->order, tt->tap[r].frac, js, h) < 0) return -1;
        for (int m = 0; m < per; m++) {
            int off = tt->tap[r].delay + js[m];
            if (off < 0) return -1;             /* 预延迟不足: 抽头落到负样本 */
            if (off == 0) return -1;            /* 反射落在直达声上 = 零延迟 */
            ir->offset[k] = off;
            ir->weight[k] = (float)((double)gain[r] * h[m]);
            if (off < first) first = off;
            k++;
        }
    }
    for (int a = 1; a < k; a++) {               /* insertion sort by offset */
        int off = ir->offset[a];
        float w = ir->weight[a];
        int b = a - 1;
        while (b >= 0 && ir->offset[b] > off) {
            ir->offset[b + 1] = ir->offset[b];
            ir->weight[b + 1] = ir->weight[b];
            b--;
        }
        ir->offset[b + 1] = off;
        ir->weight[b + 1] = w;
    }
    ir->n_taps = k;
    ir->first_offset = first;
    return 0;
}

/* Audio buffer for the reflection IR, or -1 if it does not fit. */
static int ir_synth(const Ir* ir, float* out, int cap) {
    if (ir == NULL || out == NULL || ir->n_taps <= 0) return -1;
    int len = ir->offset[ir->n_taps - 1] + 1;
    if (len > cap) return -1;
    memset(out, 0, (size_t)len * sizeof(float));
    for (int i = 0; i < ir->n_taps; i++) out[ir->offset[i]] += ir->weight[i];
    return len;
}

/* dry + reflections; wet[0] is the direct sound, reflections delayed. */
static int ir_apply(const Ir* ir, const float* dry, float* wet, int n) {
    if (ir == NULL || dry == NULL || wet == NULL || n <= 0) return -1;
    for (int k = 0; k < ir->n_taps; k++)
        if (ir->offset[k] <= 0) return -1;
    for (int i = 0; i < n; i++) wet[i] = dry[i];
    for (int k = 0; k < ir->n_taps; k++) {
        int off = ir->offset[k];
        for (int i = off; i < n; i++) wet[i] += ir->weight[k] * dry[i - off];
    }
    return n;
}

/* gain-weighted tap centroid within [lo, hi] */
static double ir_centroid(const Ir* ir, int lo, int hi) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < ir->n_taps; i++) {
        if (ir->offset[i] < lo || ir->offset[i] > hi) continue;
        num += (double)ir->offset[i] * ir->weight[i];
        den += (double)ir->weight[i];
    }
    return den != 0.0 ? num / den : -1.0;
}

/* 500 Hz + 3 kHz + 7.9 kHz probe, total amplitude 1.0, bandlimited to 7.9 kHz */
static double probe(double t) {
    return 0.5 * sin(2.0 * PI * 500.0 * t)
         + 0.3 * sin(2.0 * PI * 3000.0 * t + 0.4)
         + 0.2 * sin(2.0 * PI * 7900.0 * t + 1.1);
}

static double measured_interp_error(int order, double D) {
    int js[MAX_ORDER + 1];
    double h[MAX_ORDER + 1];
    int per = frac_taps(order, D - floor(D), js, h);
    if (per < 0) return -1.0;
    int d = (int)floor(D);
    double worst = 0.0;
    for (int k = d + MAX_ORDER + 4; k < d + 4096; k++) {
        double y = 0.0;
        for (int m = 0; m < per; m++)
            y += h[m] * probe((double)(k - d - js[m]) / FS);
        double ref = probe(((double)k - D) / FS);
        double e = fabs(y - ref);
        if (e > worst) worst = e;
    }
    return worst;
}

/* Worst-case phase-delay deviation (samples) implied by the error bound.
 * If |R - e^{-j w mu}| <= e, then |arg R + w mu| <= asin(e / (1 - e)). */
static double phase_error_bound(int order, double freq, double amp) {
    if (freq <= 0.0 || freq >= FS * 0.5) return -1.0;
    double w = 2.0 * PI * freq / FS;
    double e = frac_error_bound(order, freq, amp);
    if (e < 0.0 || e >= 0.5) return -1.0;
    return asin(e / (1.0 - e)) / w;
}

/* Phase-delay deviation from mu at freq, in samples. */
static double phase_delay_dev(int order, double mu, double freq) {
    int js[MAX_ORDER + 1];
    double h[MAX_ORDER + 1];
    int per = frac_taps(order, mu, js, h);
    if (per < 0 || freq <= 0.0 || freq >= FS * 0.5) return -1.0;
    double w = 2.0 * PI * freq / FS;
    double re = 0.0, im = 0.0;
    for (int m = 0; m < per; m++) {
        re += h[m] * cos(w * (mu - (double)js[m]));
        im += h[m] * sin(w * (mu - (double)js[m]));
    }
    return -atan2(im, re) / w;      /* 已扣掉 mu, 直接是偏差 */
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void print_table(const TapTable* tt, const Ir* ir, const float* dist) {
    int half = (tt->order + 1) / 2;
    for (int r = 0; r < tt->n; r++) {
        const Tap* t = &tt->tap[r];
        double D = (double)t->delay + t->frac;
        int lo = t->delay - half + 1, hi = t->delay + half;
        double cent = ir_centroid(ir, lo, hi);
        printf("tap%d_dist_m=%.6f\n", r, dist[r]);
        printf("tap%d_delay_ms=%.6f\n", r, D / FS * 1000.0);
        printf("tap%d_delay_samples=%.6f\n", r, D);
        printf("tap%d_int=%d\n", r, t->delay);
        printf("tap%d_frac=%.6f\n", r, t->frac);
        printf("tap%d_dist_back_m=%.6f\n", r, distance_from_tap(t, PRE_MS, FS));
        printf("tap%d_centroid_samples=%.6f\n", r, cent);
        printf("tap%d_centroid_dev=%.3e\n", r, cent - D);
        for (int i = 0; i < ir->n_taps; i++) {
            if (ir->offset[i] < lo || ir->offset[i] > hi) continue;
            printf("tap%d_ir%d=%d:%.9f\n", r, ir->offset[i] - lo,
                   ir->offset[i], ir->weight[i]);
        }
    }
}

int main(void) {
    const float dist[3] = { 2.5f, 4.1f, 6.3f };
    const float gain[3] = { 0.50f, 0.35f, 0.25f };
    const int nref = 3;
    const int order = 3;

    TapTable tt;
    Ir ir;
    if (tap_table_build(dist, gain, nref, order, PRE_MS, FS, &tt) != 0 ||
        ir_build(&tt, gain, &ir) != 0) {
        check("build", 0);
        printf("all_ok=0\n");
        return 1;
    }
    print_table(&tt, &ir, dist);
    printf("ir_len=%d\n", ir.offset[ir.n_taps - 1] + 1);
    printf("ir_first_offset=%d\n", ir.first_offset);
    printf("predelay_samples=%.0f\n", PRE_MS * 0.001 * FS);

    /* --- tap time <-> distance --- */
    double back_err = 0.0;
    for (int r = 0; r < tt.n; r++) {
        double back = distance_from_tap(&tt.tap[r], PRE_MS, FS);
        double e = fabs(back - (double)dist[r]);
        if (e > back_err) back_err = e;
    }
    check("dist_roundtrip", back_err < 1e-9);

    double D0 = (PRE_MS * 0.001 + 2.5 / SOS) * FS;
    check("tap_formula", fabs(((double)tt.tap[0].delay + tt.tap[0].frac) - D0) < 1e-9);
    check("predelay_gap", D0 >= PRE_MS * 0.001 * FS && ir.first_offset >= 240);

    Tap t;
    TapTable bad;
    float zero = 0.0f, tiny = 0.001f, one = 1.0f;
    check("zero_delay_invalid", tap_from_distance(0.0, 0.0, FS, 3, &t) == -1);
    check("neg_dist_invalid", tap_from_distance(-1.0, PRE_MS, FS, 3, &t) == -1);
    check("even_order_invalid", tap_from_distance(2.0, PRE_MS, FS, 2, &t) == -1);
    check("zero_dist_build_invalid",
          tap_table_build(&zero, &one, 1, 1, 0.0, FS, &bad) == -1);
    check("tap_on_direct_invalid",
          tap_table_build(&tiny, &one, 1, 1, 0.0, FS, &bad) == -1);
    check("order5_build", tap_table_build(dist, gain, nref, 5, PRE_MS, FS, &bad) == 0);

    /* --- interpolation weights: DC gain, first moment --- */
    int js[MAX_ORDER + 1];
    double h[MAX_ORDER + 1];
    check("mu_range_invalid", frac_taps(3, 1.0, js, h) == -1);
    check("zero_delay_frac_invalid", frac_taps(3, -1e-19, js, h) == -1);
    double max_dc = 0.0, max_mom = 0.0;
    for (int o = 1; o <= MAX_ORDER; o += 2) {
        for (int k = 0; k < 100; k++) {
            double mu = (double)k / 100.0;
            int per = frac_taps(o, mu, js, h);
            double s = 0.0, m = 0.0;
            for (int i = 0; i < per; i++) { s += h[i]; m += (double)js[i] * h[i]; }
            if (fabs(s - 1.0) > max_dc) max_dc = fabs(s - 1.0);
            if (fabs(m - mu) > max_mom) max_mom = fabs(m - mu);
        }
    }
    check("unit_dc_all_orders", max_dc < 1e-9);
    check("centroid_mu_all_orders", max_mom < 1e-9);

    check("pmax_o1", fabs(lagrange_pmax(1) - 0.25) < 1e-12);
    check("pmax_o3", fabs(lagrange_pmax(3) - 0.5625) < 1e-12);
    check("pmax_o5", fabs(lagrange_pmax(5) - 3.515625) < 1e-12);
    check("bad_order_pmax", lagrange_pmax(2) < 0.0);

    /* --- IR structure --- */
    check("ir_tap_count", ir.n_taps == (order + 1) * nref);
    check("ir_first_offset", ir.first_offset == tt.tap[0].delay - 1);
    double cent_dev = 0.0;
    for (int r = 0; r < tt.n; r++) {
        int lo = tt.tap[r].delay - (order + 1) / 2 + 1;
        int hi = tt.tap[r].delay + (order + 1) / 2;
        double D = (double)tt.tap[r].delay + tt.tap[r].frac;
        double e = fabs(ir_centroid(&ir, lo, hi) - D);
        if (e > cent_dev) cent_dev = e;
    }
    check("ir_centroid_matches", cent_dev < 1e-3);

    /* order-1 IR must reproduce the analytic linear weights */
    TapTable t1;
    Ir ir1;
    tap_table_build(dist, gain, 1, 1, PRE_MS, FS, &t1);
    ir_build(&t1, gain, &ir1);
    double mu0 = t1.tap[0].frac;
    check("linear_weight0",
          fabsf(ir1.weight[0] - (float)((double)gain[0] * (1.0 - mu0))) < 1e-6f);
    check("linear_weight1",
          fabsf(ir1.weight[1] - (float)((double)gain[0] * mu0)) < 1e-6f);

    /* --- reflection IR buffer and dry mix --- */
    static float irsig[4096], dry[4096], wet[4096];
    int len = ir_synth(&ir, irsig, 4096);
    check("synth_len", len == ir.offset[ir.n_taps - 1] + 1);
    check("synth_silent_before_taps",
          irsig[0] == 0.0f && irsig[ir.first_offset - 1] == 0.0f &&
          fabsf(irsig[ir.first_offset]) > 1e-6f);
    check("synth_cap_invalid", ir_synth(&ir, irsig, 64) == -1);
    check("apply_bad", ir_apply(NULL, dry, wet, 4096) == -1);

    memset(dry, 0, sizeof(dry));
    dry[0] = 1.0f;
    int rn = ir_apply(&ir, dry, wet, 4096);
    int first_nz = -1;
    for (int i = 1; i < 4096; i++)
        if (fabsf(wet[i]) > 1e-7f) { first_nz = i; break; }
    float sum = 0.0f;
    for (int i = 0; i < 4096; i++) sum += wet[i];
    check("render_ok", rn == 4096 && wet[0] == 1.0f);
    check("render_predelay_gap", first_nz == ir.first_offset && first_nz >= 240);
    check("render_gain_sum", fabsf(sum - 2.1f) < 1e-4f);

    /* --- error bounds vs measured error, per interpolation order --- */
    double prev_bound = 1e9;
    for (int o = 1; o <= MAX_ORDER; o += 2) {
        char name[64];
        double bound = frac_error_bound(o, 7900.0, 1.0);
        double err = measured_interp_error(o, 100.3);
        printf("interp_bound_o%d=%.9f\n", o, bound);
        printf("interp_err_o%d=%.9f\n", o, err);
        printf("interp_ratio_o%d=%.6f\n", o, err / bound);
        snprintf(name, sizeof(name), "interp_o%d_within_bound", o);
        check(name, err >= 0.0 && err < bound);
        snprintf(name, sizeof(name), "bound_shrinks_o%d", o);
        check(name, bound < prev_bound);
        prev_bound = bound;
        double dev = phase_delay_dev(o, 0.3, 4000.0);
        double lim = phase_error_bound(o, 4000.0, 1.0);
        printf("phase_dev_o%d_4k=%.3e\n", o, dev);
        printf("phase_lim_o%d_4k=%.3e\n", o, lim);
        snprintf(name, sizeof(name), "phase_o%d_within_bound", o);
        check(name, fabs(dev) <= lim);
    }
    check("bad_order_bound", frac_error_bound(2, 4000.0, 1.0) < 0.0);
    check("fmax_over_nyquist", frac_error_bound(1, FS * 0.5 + 1.0, 1.0) < 0.0);
    check("bound_amp_linear",
          fabs(frac_error_bound(3, 4000.0, 0.5) - 0.5 * frac_error_bound(3, 4000.0, 1.0))
          < 1e-12);
    check("bound_monotone_fmax",
          frac_error_bound(3, 2000.0, 1.0) < frac_error_bound(3, 4000.0, 1.0));

    /* published bound table: fmax in {1,4,8} kHz, amp = 1 */
    static const double flist[3] = { 1000.0, 4000.0, 8000.0 };
    for (int fi = 0; fi < 3; fi++)
        for (int o = 1; o <= MAX_ORDER; o += 2)
            printf("err_bound_o%d_f%d=%.6e\n", o, (int)flist[fi],
                   frac_error_bound(o, flist[fi], 1.0));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
