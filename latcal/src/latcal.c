
/* latcal: end-to-end latency calibration of two audio chains (BGM vs judge
 * sound), cross-chain compensation, and buffer-length auto adjustment.
 *
 * Each chain is modelled as block buffering + fixed algorithmic delay:
 *   BGM chain   : double-buffered stream, block delay = 2*B, fx = 192 smp
 *   judge chain : one-block synth path, block delay = B,   fx =  32 smp
 * so the true latency gap grows with the buffer length B and a single fixed
 * compensation cannot hold across buffer settings.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define SR              48000.0   /* sample rate */
#define MAX_DELAY       4096      /* samples, per delay line */
#define CAPTURE_BLOCKS  64        /* coarse probe observation window */
#define DETECT          0.5f      /* impulse peak threshold */
#define IMPULSE_TRIALS  3
#define SPREAD_MAX      2         /* samples: repeated trials must agree */
#define RESIDUAL_MAX_MS 1.0       /* post-compensation hard limit */
#define N_BUF           3
#define FX_BG           192
#define FX_JUD          32

static const int BUF_SIZES[N_BUF] = {64, 256, 1024};

/* ---------- generic delay line ---------- */

typedef struct { float* buf; int cap; int pos; } DelayLine;

static int dl_init(DelayLine* d, int cap) {
    d->buf = (float*)calloc((size_t)cap, sizeof(float));
    if (!d->buf) return -1;
    d->cap = cap;
    d->pos = 0;
    return 0;
}

static void dl_free(DelayLine* d) {
    free(d->buf);
    d->buf = NULL;
}

static void dl_reset(DelayLine* d) {
    if (d->buf) memset(d->buf, 0, (size_t)d->cap * sizeof(float));
    d->pos = 0;
}

/* push x, read the sample `delay` positions behind; 0 <= delay < cap */
static float dl_push(DelayLine* d, float x, int delay) {
    d->buf[d->pos] = x;
    int idx = d->pos - delay;
    if (idx < 0) idx += d->cap;
    float y = d->buf[idx];
    if (++d->pos == d->cap) d->pos = 0;
    return y;
}

/* ---------- audio chain ---------- */

typedef struct {
    int block_frames;          /* B */
    int block_factor;          /* 2 = double buffer (BGM), 1 = judge */
    int fx_delay;              /* fixed algorithmic delay */
    int blk_delay;             /* block_factor * B */
    int comp_delay;            /* compensation attached to this chain */
    DelayLine blk, fx, comp;
} Chain;

static int chain_init(Chain* c, int B, int factor, int fx) {
    if (B <= 0 || factor <= 0 || fx < 0) return -1;
    memset(c, 0, sizeof(*c));
    c->block_frames = B;
    c->block_factor = factor;
    c->fx_delay = fx;
    c->blk_delay = factor * B;
    if (c->blk_delay + fx >= MAX_DELAY) return -1;
    if (dl_init(&c->blk, MAX_DELAY) ||
        dl_init(&c->fx, MAX_DELAY) ||
        dl_init(&c->comp, MAX_DELAY)) {
        dl_free(&c->blk);
        dl_free(&c->fx);
        dl_free(&c->comp);
        return -1;
    }
    return 0;
}

static void chain_reset(Chain* c) {
    dl_reset(&c->blk);
    dl_reset(&c->fx);
    dl_reset(&c->comp);
}

static void chain_free(Chain* c) {
    dl_free(&c->blk);
    dl_free(&c->fx);
    dl_free(&c->comp);
}

static float chain_process(Chain* c, float x) {
    float y = dl_push(&c->blk, x, c->blk_delay);
    y = dl_push(&c->fx, y, c->fx_delay);
    if (c->comp_delay > 0)
        y = dl_push(&c->comp, y, c->comp_delay);
    return y;
}

/* ---------- latency probes ---------- */

/* Reliable method: inject a unit impulse and locate the output peak.
 * inject=0 simulates a missed trigger and must never report a latency.
 * Returns latency in samples, -1 when no peak is detected, -2 when the
 * peak fails the gain/purity checks. */
static int probe_impulse_once(Chain* c, int inject) {
    chain_reset(c);
    int peak_idx = -1;
    float peak = 0.0f, other = 0.0f;
    for (int i = 0; i < MAX_DELAY * 2; i++) {
        float y = chain_process(c, (inject && i == 0) ? 1.0f : 0.0f);
        float a = fabsf(y);
        if (a > peak) { peak = a; peak_idx = i; }
        else if (a > other) other = a;
    }
    if (peak_idx < 0 || peak < DETECT) return -1;
    if (fabsf(peak - 1.0f) > 1e-4f) return -2;   /* chain must be unity gain */
    if (other > 1e-5f) return -2;                /* no echoes/ringing */
    return peak_idx;
}

/* Combine repeated trials: reject when the spread is too wide, otherwise
 * take the nearest-integer mean (equals the median for tight trials). */
static int trials_decide(const int* t, int n, int* out) {
    if (n <= 0) return -1;
    int lo = t[0], hi = t[0];
    long sum = 0;
    for (int i = 0; i < n; i++) {
        if (t[i] < lo) lo = t[i];
        if (t[i] > hi) hi = t[i];
        sum += t[i];
    }
    if (hi - lo > SPREAD_MAX) return -1;
    *out = (int)((sum + n / 2) / n);
    return 0;
}

static int probe_impulse(Chain* c, int trials) {
    if (trials <= 0 || trials > 8) return -1;
    int t[8];
    for (int i = 0; i < trials; i++) {
        t[i] = probe_impulse_once(c, 1);
        if (t[i] < 0) return -1;
    }
    int v;
    return trials_decide(t, trials, &v) == 0 ? v : -3;
}

/* Faulty method: the probe can only timestamp at block-callback boundaries,
 * so every reading is rounded up to a whole block. Quantisation error grows
 * with B — fine at small buffers, wrong by several ms at large ones. */
static int probe_coarse(Chain* c) {
    chain_reset(c);
    for (int blk = 1; blk <= CAPTURE_BLOCKS; blk++) {
        int seen = 0;
        for (int i = 0; i < c->block_frames; i++) {
            float y = chain_process(c, (blk == 1 && i == 0) ? 1.0f : 0.0f);
            if (fabsf(y) > DETECT) seen = 1;
        }
        if (seen) return blk * c->block_frames;   /* callback boundary */
    }
    return -1;
}

/* ---------- calibration ---------- */

typedef struct {
    int block_frames;
    int lat_bg_smp, lat_jud_smp;        /* accepted measurements */
    int coarse_bg_smp, coarse_jud_smp;  /* first (faulty) method readings */
    int comp_smp;                       /* delay attached to judge chain */
    int residual_smp;                   /* measured arrival gap after comp */
    int attempts;                       /* 1 = coarse passed, 2 = re-measured */
    int ok;
} CalPoint;

typedef struct {
    double sr;
    int n;
    CalPoint pts[N_BUF];
} CalTable;

static double to_ms(int smp, double sr) { return (double)smp * 1000.0 / sr; }

/* Re-probe both chains end-to-end with compensation already attached and
 * return the absolute arrival-time gap in samples. */
static int measure_residual(Chain* bg, Chain* jud) {
    int a = probe_impulse(bg, 1);
    int b = probe_impulse(jud, 1);
    if (a < 0 || b < 0) return -1;
    return abs(a - b);
}

/* Calibrate one buffer length. Return codes:
 *   0 ok; -1 probe failure; -2 judge chain slower than BGM (model broken);
 *  -3 reliable probe failed on re-measure;
 *  -4 residual still >= 1 ms after re-measure: refuse, do not add more comp. */
static int cal_point(Chain* bg, Chain* jud, int B, CalPoint* p) {
    memset(p, 0, sizeof(*p));
    p->block_frames = B;

    /* attempt 1: cheap block-quantised probe */
    int l_bg = probe_coarse(bg);
    int l_jud = probe_coarse(jud);
    if (l_bg < 0 || l_jud < 0) return -1;
    p->coarse_bg_smp = l_bg;
    p->coarse_jud_smp = l_jud;
    p->attempts = 1;

    int delta = l_bg - l_jud;
    if (delta < 0) return -2;
    jud->comp_delay = delta;            /* overwrite, never accumulate */

    int res = measure_residual(bg, jud);
    if (res < 0) return -1;

    if (to_ms(res, SR) >= RESIDUAL_MAX_MS) {
        /* The measurement method is at fault: discard its result entirely
         * and re-measure with the reliable impulse probe. */
        jud->comp_delay = 0;
        int r_bg = probe_impulse(bg, IMPULSE_TRIALS);
        int r_jud = probe_impulse(jud, IMPULSE_TRIALS);
        if (r_bg < 0 || r_jud < 0) return -3;
        delta = r_bg - r_jud;
        if (delta < 0) return -2;
        jud->comp_delay = delta;        /* fresh value, not coarse + extra */
        res = measure_residual(bg, jud);
        if (res < 0) return -1;
        if (to_ms(res, SR) >= RESIDUAL_MAX_MS) return -4;
        l_bg = r_bg;
        l_jud = r_jud;
        p->attempts = 2;
    }

    p->lat_bg_smp = l_bg;
    p->lat_jud_smp = l_jud;
    p->comp_smp = delta;
    p->residual_smp = res;
    p->ok = 1;
    return 0;
}

/* Compensation in samples for an arbitrary buffer length, interpolated
 * between calibration points. Out of range means recalibrate: -1. */
static int comp_auto_smp(const CalTable* t, int B) {
    if (t->n < 2 || B < t->pts[0].block_frames ||
        B > t->pts[t->n - 1].block_frames)
        return -1;
    for (int i = 0; i < t->n - 1; i++) {
        const CalPoint* a = &t->pts[i];
        const CalPoint* b = &t->pts[i + 1];
        if (B >= a->block_frames && B <= b->block_frames) {
            if (B == a->block_frames) return a->comp_smp;
            if (B == b->block_frames) return b->comp_smp;
            double frac = (double)(B - a->block_frames) /
                          (double)(b->block_frames - a->block_frames);
            double v = (double)a->comp_smp +
                       frac * (double)(b->comp_smp - a->comp_smp);
            return (int)lround(v);
        }
    }
    return -1;
}

/* ---------- self test ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    CalTable tbl;
    tbl.sr = SR;
    tbl.n = 0;

    printf("== end-to-end latency, sr=%.0f Hz ==\n", SR);
    for (int i = 0; i < N_BUF; i++) {
        int B = BUF_SIZES[i];
        Chain bg, jud;
        check("chain_init", chain_init(&bg, B, 2, FX_BG) == 0 &&
                            chain_init(&jud, B, 1, FX_JUD) == 0);

        int rc = cal_point(&bg, &jud, B, &tbl.pts[i]);
        check("cal_ok", rc == 0 && tbl.pts[i].ok);
        tbl.n++;

        CalPoint* p = &tbl.pts[i];
        printf("buf=%-5d bg=%4d smp/%7.3f ms  judge=%4d smp/%7.3f ms  "
               "comp=%4d smp/%7.3f ms  residual=%d smp/%.3f ms  attempts=%d\n",
               B, p->lat_bg_smp, to_ms(p->lat_bg_smp, SR),
               p->lat_jud_smp, to_ms(p->lat_jud_smp, SR),
               p->comp_smp, to_ms(p->comp_smp, SR),
               p->residual_smp, to_ms(p->residual_smp, SR), p->attempts);

        /* six measured groups (two chains x three buffers) must match the
         * analytic chain model; drop compensation before probing */
        jud.comp_delay = 0;
        check("bg_matches_model", probe_impulse(&bg, IMPULSE_TRIALS) == 2 * B + FX_BG);
        check("jud_matches_model", probe_impulse(&jud, IMPULSE_TRIALS) == B + FX_JUD);

        /* after a re-measurement compensation is the true gap, never
         * coarse-plus-an-extra-kick; the coarse-only point is allowed any
         * value as long as its residual stays under 1 ms (checked above) */
        if (p->attempts == 2)
            check("comp_exact", p->comp_smp == B + (FX_BG - FX_JUD));
        check("residual_under_1ms", to_ms(p->residual_smp, SR) < RESIDUAL_MAX_MS);

        chain_free(&bg);
        chain_free(&jud);
    }

    /* faulty probe exposed: its quantisation error grows with the block size */
    {
        Chain bg64, jud64, bg1k, jud1k;
        chain_init(&bg64, 64, 2, FX_BG);
        chain_init(&jud64, 64, 1, FX_JUD);
        chain_init(&bg1k, 1024, 2, FX_BG);
        chain_init(&jud1k, 1024, 1, FX_JUD);

        int c64 = probe_coarse(&bg64) - probe_coarse(&jud64);
        int c1k = probe_coarse(&bg1k) - probe_coarse(&jud1k);
        check("coarse_small_within_1ms",
              abs(c64 - (64 + FX_BG - FX_JUD)) * 1000.0 / SR < RESIDUAL_MAX_MS);
        check("coarse_large_wrong",
              abs(c1k - (1024 + FX_BG - FX_JUD)) * 1000.0 / SR >= RESIDUAL_MAX_MS);
        check("small_buffer_no_remeasure", tbl.pts[0].attempts == 1);
        check("remeasured_large_buffer", tbl.pts[2].attempts == 2 &&
              tbl.pts[2].residual_smp == 0);

        /* missed trigger must not invent a latency */
        check("no_peak_rejected", probe_impulse_once(&jud1k, 0) == -1);

        chain_free(&bg64);
        chain_free(&jud64);
        chain_free(&bg1k);
        chain_free(&jud1k);
    }

    /* jittery trials rejected instead of averaged into a bogus value */
    {
        int tight[IMPULSE_TRIALS] = {320, 321, 320};
        int wide[IMPULSE_TRIALS] = {320, 325, 318};
        int v;
        check("trials_tight", trials_decide(tight, IMPULSE_TRIALS, &v) == 0 &&
              v == 320);
        check("trials_wide_rejected", trials_decide(wide, IMPULSE_TRIALS, &v) == -1);
    }

    /* auto adjustment: compensation follows buffer length by interpolation,
     * and the interpolated value actually lands within 1 ms end-to-end */
    {
        int c512 = comp_auto_smp(&tbl, 512);
        check("auto_interp", c512 == 512 + FX_BG - FX_JUD);
        Chain bg, jud;
        chain_init(&bg, 512, 2, FX_BG);
        chain_init(&jud, 512, 1, FX_JUD);
        jud.comp_delay = c512;
        int res = measure_residual(&bg, &jud);
        check("auto_residual_under_1ms", res >= 0 && to_ms(res, SR) < RESIDUAL_MAX_MS);
        check("auto_endpoints", comp_auto_smp(&tbl, 64) == tbl.pts[0].comp_smp &&
              comp_auto_smp(&tbl, 1024) == tbl.pts[2].comp_smp);
        check("auto_out_of_range", comp_auto_smp(&tbl, 32) == -1 &&
              comp_auto_smp(&tbl, 2048) == -1);
        chain_free(&bg);
        chain_free(&jud);
    }

    /* bad construction is rejected */
    {
        Chain bad;
        check("bad_chain", chain_init(&bad, 0, 2, FX_BG) == -1);
        check("huge_chain", chain_init(&bad, MAX_DELAY, 2, FX_BG) == -1);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
