/* gainstage: unified dB-domain gain budget for the five-stage bus
 * (resample -> eq -> limiter -> reverb -> pan) with a full-chain peak audit.
 *
 * The bug this fixes: every stage used to multiply its own gain with no
 * global accounting, so the bus level drifted upward and clipped even with
 * the mixer faders at zero. Here every stage's trim comes from one budget
 * table in the dB domain; the table is queryable and the chain is audited
 * with a drum+vocal signal, peaks recorded before/after each stage.
 *
 * The limiter is included verbatim from limiter/src/limiter.c: its trigger
 * behavior is untouched, only the gain structure ahead of it is budgeted.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define main limiter_selftest_main /* keep the limiter self-test out of our main */
#include "../../limiter/src/limiter.c"
#undef main

#define SR 48000.0f
#define N 4096
#define TWO_PI 6.28318530718f

static float db_to_lin(float db) { return powf(10.0f, db / 20.0f); }
static float lin_to_db(float x) { return x > 1e-9f ? 20.0f * log10f(x) : -180.0f; }

/* ---------------- gain budget table (dB domain) ---------------- */

#define STAGES 5
#define BUDGET_DB 8.0f /* total headroom the five stages together may spend */

typedef struct {
    const char* name;
    float trim_db;       /* gain applied at this stage (never > 0) */
    float cumulative_db; /* sum of trims from chain input to this stage's output */
    float consumed_db;   /* budget spent up to and including this stage */
    float remaining_db;  /* budget left after this stage */
} StageGain;

/* per-stage trim allocation; sum must equal -BUDGET_DB */
static const char* k_stage_name[STAGES] = { "resample", "eq", "limiter", "reverb", "pan" };
static const float k_stage_trim_db[STAGES] = { -1.0f, -4.0f, 0.0f, -2.0f, -1.0f };

static void build_table(StageGain t[STAGES]) {
    float cum = 0.0f;
    for (int i = 0; i < STAGES; i++) {
        cum += k_stage_trim_db[i];
        t[i].name = k_stage_name[i];
        t[i].trim_db = k_stage_trim_db[i];
        t[i].cumulative_db = cum;
        t[i].consumed_db = -cum;
        t[i].remaining_db = BUDGET_DB + cum;
    }
}

/* queryable lookup: returns 1 and fills *out when the stage exists */
static int gainstage_find(const char* name, StageGain* out) {
    StageGain t[STAGES];
    build_table(t);
    for (int i = 0; i < STAGES; i++) {
        if (strcmp(t[i].name, name) == 0) {
            if (out) *out = t[i];
            return 1;
        }
    }
    return 0;
}

static void print_table(void) {
    StageGain t[STAGES];
    build_table(t);
    printf("gain_budget_db=%.1f\n", BUDGET_DB);
    printf("%-9s %8s %13s %12s %13s\n", "stage", "trim_dB", "cumulative_dB", "consumed_dB", "remaining_dB");
    for (int i = 0; i < STAGES; i++)
        printf("%-9s %8.1f %13.1f %12.1f %13.1f\n",
               t[i].name, t[i].trim_db, t[i].cumulative_db, t[i].consumed_db, t[i].remaining_db);
}

/* ---------------- the five stages ---------------- */

/* 1. resample: linear interpolation (ratio 1.0 here; length preserved) */
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

/* 2. eq: RBJ peaking biquad; the +4 dB design boost is what the stage trim pays for */
typedef struct { float b0, b1, b2, a1, a2, x1, x2, y1, y2; } Biquad;

static void eq_init(Biquad* q, float f0, float q_bw, float boost_db) {
    float A = powf(10.0f, boost_db / 40.0f);
    float w0 = TWO_PI * f0 / SR;
    float alpha = sinf(w0) / (2.0f * q_bw);
    float a0 = 1.0f + alpha / A;
    q->b0 = (1.0f + alpha * A) / a0;
    q->b1 = (-2.0f * cosf(w0)) / a0;
    q->b2 = (1.0f - alpha * A) / a0;
    q->a1 = (-2.0f * cosf(w0)) / a0;
    q->a2 = (1.0f - alpha / A) / a0;
    q->x1 = q->x2 = q->y1 = q->y2 = 0.0f;
}

static float eq_tick(Biquad* q, float x) {
    float y = q->b0 * x + q->b1 * q->x1 + q->b2 * q->x2 - q->a1 * q->y1 - q->a2 * q->y2;
    q->x2 = q->x1; q->x1 = x;
    q->y2 = q->y1; q->y1 = y;
    return y;
}

/* 4. reverb: single feedback comb, dry + wet */
#define COMB_D 1776 /* ~37 ms at 48 kHz */
typedef struct { float buf[COMB_D]; int pos; float feedback; float wet; } Comb;

static void comb_init(Comb* c, float feedback, float wet) {
    memset(c->buf, 0, sizeof(c->buf));
    c->pos = 0;
    c->feedback = feedback;
    c->wet = wet;
}

static float comb_tick(Comb* c, float x) {
    float v = x + c->feedback * c->buf[c->pos];
    c->buf[c->pos] = v;
    c->pos = (c->pos + 1) % COMB_D;
    return x + c->wet * v;
}

/* 5. pan: constant-power, returns peak over both channels */
static void pan_gains(float pan, float* left, float* right) {
    float a = (pan + 1.0f) * 0.25f * 3.14159265f;
    *left = cosf(a);
    *right = sinf(a);
}

/* ---------------- chain run with peak probes ---------------- */

typedef struct { float in_peak[STAGES]; float out_peak[STAGES]; float final_peak; } ChainAudit;

static float peak_of(const float* v, int n) {
    float p = 0.0f;
    for (int i = 0; i < n; i++) { float a = fabsf(v[i]); if (a > p) p = a; }
    return p;
}

/* use_trim=1 applies the budget table; use_trim=0 is the old unmanaged behavior */
static ChainAudit run_chain(const float* src, int n, int use_trim) {
    static float s1[N], s2[N], s3[N], s4[N];
    ChainAudit au;
    memset(&au, 0, sizeof(au));

    StageGain t[STAGES];
    build_table(t);
    float g[STAGES];
    for (int i = 0; i < STAGES; i++) g[i] = use_trim ? db_to_lin(t[i].trim_db) : 1.0f;

    /* 1. resample */
    int n1 = resample_linear(src, n, 1.0, s1, N);
    au.in_peak[0] = peak_of(src, n);
    for (int i = 0; i < n1; i++) s1[i] *= g[0];
    au.out_peak[0] = peak_of(s1, n1);

    /* 2. eq: +4 dB peaking boost at 800 Hz, then trim */
    Biquad q;
    eq_init(&q, 800.0f, 1.0f, 4.0f);
    au.in_peak[1] = au.out_peak[0];
    for (int i = 0; i < n1; i++) s2[i] = eq_tick(&q, s1[i]) * g[1];
    au.out_peak[1] = peak_of(s2, n1);

    /* 3. limiter: verbatim module, ceiling 0.9, trim (0 dB by allocation) */
    Limiter lim;
    lim_init(&lim, 5.0f, 50.0f, SR);
    Report rep;
    au.in_peak[2] = au.out_peak[1];
    lim_process(&lim, s2, s3, n1, 0.9f, &rep);
    for (int i = 0; i < n1; i++) s3[i] *= g[2];
    au.out_peak[2] = peak_of(s3, n1);

    /* 4. reverb: comb, wet 0.5, then trim */
    Comb c;
    comb_init(&c, 0.35f, use_trim ? 0.5f : 1.0f);
    au.in_peak[3] = au.out_peak[2];
    for (int i = 0; i < n1; i++) s4[i] = comb_tick(&c, s3[i]) * g[3];
    au.out_peak[3] = peak_of(s4, n1);

    /* 5. pan: constant-power to stereo, peak over both channels, then trim */
    float gl, gr;
    pan_gains(0.3f, &gl, &gr);
    au.in_peak[4] = au.out_peak[3];
    float p = 0.0f;
    for (int i = 0; i < n1; i++) {
        float l = s4[i] * gl * g[4], r = s4[i] * gr * g[4];
        if (fabsf(l) > p) p = fabsf(l);
        if (fabsf(r) > p) p = fabsf(r);
    }
    au.out_peak[4] = p;
    au.final_peak = p;
    return au;
}

/* ---------------- test material: drums + vocal ---------------- */

static void make_material(float* out, int n) {
    for (int i = 0; i < n; i++) {
        float t = (float)i / SR;
        /* kick every 0.3 s: 55 Hz sine with fast decay */
        float drum = 0.0f;
        float local = fmodf(t, 0.3f);
        drum = sinf(TWO_PI * 55.0f * local) * expf(-local * 25.0f);
        /* vocal: 220 Hz + harmonics with slow amplitude modulation */
        float vib = 0.7f + 0.3f * sinf(TWO_PI * 4.0f * t);
        float vocal = (0.6f * sinf(TWO_PI * 220.0f * t)
                     + 0.3f * sinf(TWO_PI * 440.0f * t)
                     + 0.15f * sinf(TWO_PI * 880.0f * t)) * vib;
        out[i] = 0.8f * drum + 0.5f * vocal;
    }
    float p = peak_of(out, n);
    for (int i = 0; i < n; i++) out[i] = out[i] / p * 0.8f; /* chain head at -1.9 dBFS */
}

/* ---------------- self test ---------------- */

static int gs_failures = 0;
static void gs_check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) gs_failures++;
}

int main(void) {
    print_table();

    /* table integrity */
    StageGain t[STAGES];
    build_table(t);
    float cum = 0.0f;
    int sums_ok = 1, no_boost = 1, remain_ok = 1;
    for (int i = 0; i < STAGES; i++) {
        cum += k_stage_trim_db[i];
        if (fabsf(t[i].cumulative_db - cum) > 1e-6f) sums_ok = 0;
        if (t[i].trim_db > 0.0f) no_boost = 0;
        if (fabsf(t[i].consumed_db + t[i].remaining_db - BUDGET_DB) > 1e-6f) remain_ok = 0;
    }
    gs_check("table_cumulative", sums_ok);
    gs_check("table_no_boost", no_boost);
    gs_check("table_budget_closed", remain_ok && fabsf(t[STAGES - 1].remaining_db) < 1e-6f);

    /* queryable */
    StageGain row;
    gs_check("find_eq", gainstage_find("eq", &row) && fabsf(row.trim_db - (-4.0f)) < 1e-6f);
    gs_check("find_pan_remaining", gainstage_find("pan", &row) && fabsf(row.remaining_db) < 1e-6f);
    gs_check("find_unknown", !gainstage_find("chorus", &row));

    /* full chain with drum+vocal material */
    static float src[N];
    make_material(src, N);
    float head_peak = peak_of(src, N);

    ChainAudit au = run_chain(src, N, 1);
    printf("%-9s %12s %12s\n", "stage", "in_dBFS", "out_dBFS");
    for (int i = 0; i < STAGES; i++)
        printf("%-9s %12.2f %12.2f\n", t[i].name, lin_to_db(au.in_peak[i]), lin_to_db(au.out_peak[i]));
    printf("chain_head_dBFS=%.2f\n", lin_to_db(head_peak));
    printf("chain_tail_dBFS=%.2f\n", lin_to_db(au.final_peak));

    gs_check("tail_le_head", au.final_peak <= head_peak + 1e-4f);
    gs_check("tail_under_full_scale", au.final_peak <= 1.0f);

    /* contrast: the old unmanaged gains still run away and clip */
    ChainAudit naive = run_chain(src, N, 0);
    printf("naive_tail_dBFS=%.2f\n", lin_to_db(naive.final_peak));
    gs_check("naive_would_clip", naive.final_peak > 1.0f);
    gs_check("budget_beats_naive", au.final_peak < naive.final_peak);

    /* limiter trigger behavior intact: same loud-block stimulus as its own self test */
    Limiter lim;
    float loud[BLK], lout[BLK];
    for (int i = 0; i < BLK; i++) loud[i] = 2.0f;
    lim_init(&lim, 5.0f, 50.0f, SR);
    Report rep;
    lim_process(&lim, loud, lout, BLK, 0.8f, &rep);
    gs_check("limiter_still_triggers", rep.limited_samples > 0 && rep.max_reduction > 0.0f);
    gs_check("limiter_still_caps", rep.peak_out <= 0.8f + 1e-4f);

    printf("all_ok=%d\n", gs_failures == 0);
    return gs_failures == 0 ? 0 : 1;
}
