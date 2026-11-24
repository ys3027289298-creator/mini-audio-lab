
/* limiter: peak limiter with attack/release envelope and gain reduction report.
 *
 * Release is program-dependent: a level envelope splits the release rate in
 * two - fast while the signal is loud, slow while it is quiet. The fast/slow
 * blend weight is smoothed, so when the level crosses the split point the
 * gain hands off continuously along the reduction direction instead of
 * stepping. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define BLK 128

/* level envelope and blend smoothing time constants */
#define LEVEL_ATTACK_MS 1.0f
#define LEVEL_RELEASE_MS 200.0f
#define BLEND_MS 10.0f

typedef struct {
    float attack_ms;
    float release_fast_ms; /* release while the level envelope is hot */
    float release_slow_ms; /* release while the level envelope is quiet */
    float split;           /* level threshold between fast and slow release */
    float sample_rate;
    float env;             /* gain envelope */
    float level;           /* level envelope, selects the release rate */
    float w;               /* smoothed fast/slow blend weight, 1 = fast */
    int engaged;
} Limiter;

typedef struct {
    float peak_in;
    float peak_out;
    int limited_samples;
    float max_reduction; /* deepest gain dip in the block, 1 - min(env) */
    float max_gain;      /* largest gain envelope value in the block */
    float end_gain;      /* gain envelope after the last sample */
} Report;

static void lim_init(Limiter* l, float attack_ms, float release_fast_ms,
                     float release_slow_ms, float split, float sr) {
    memset(l, 0, sizeof(*l));
    l->attack_ms = attack_ms;
    l->release_fast_ms = release_fast_ms;
    l->release_slow_ms = release_slow_ms;
    l->split = split;
    l->sample_rate = sr;
    l->env = 1.0f;
}

/* one-pole coefficients from time constants */
static float coef(float ms, float sr) {
    if (ms <= 0.0f || sr <= 0.0f) return 0.0f;
    return expf(-1.0f / (ms * 0.001f * sr));
}

static int lim_process(Limiter* l, const float* in, float* out, int n, float ceiling, Report* rep) {
    if (l == NULL || in == NULL || out == NULL || n <= 0) return -1;
    if (ceiling <= 0.0f || ceiling > 1.0f) return -1;
    if (l->attack_ms <= 0.0f || l->release_fast_ms <= 0.0f || l->release_slow_ms <= 0.0f) return -1;
    if (l->split <= 0.0f || l->split > 1.0f) return -1;

    float ca = coef(l->attack_ms, l->sample_rate);
    float cf = coef(l->release_fast_ms, l->sample_rate);
    float cs = coef(l->release_slow_ms, l->sample_rate);
    float cw = coef(BLEND_MS, l->sample_rate);
    float cla = coef(LEVEL_ATTACK_MS, l->sample_rate);
    float clr = coef(LEVEL_RELEASE_MS, l->sample_rate);
    float pin = 0.0f, pout = 0.0f, red = 0.0f, gmax = 0.0f;
    int lim = 0;

    for (int i = 0; i < n; i++) {
        float x = in[i];
        float ax = fabsf(x);
        if (ax > pin) pin = ax;
        /* level envelope: fast attack, slow decay */
        float cl = ax > l->level ? cla : clr;
        l->level = cl * l->level + (1.0f - cl) * ax;
        float target = ax > ceiling ? ceiling / ax : 1.0f;
        if (target < l->env) {
            l->env = ca * l->env + (1.0f - ca) * target;
        } else {
            /* hot level releases fast, quiet level releases slow; the
             * smoothed weight glides between the two so the gain never
             * steps at the handoff */
            float wt = l->level >= l->split ? 1.0f : 0.0f;
            l->w = cw * l->w + (1.0f - cw) * wt;
            float cr = l->w * cf + (1.0f - l->w) * cs;
            l->env = cr * l->env + (1.0f - cr) * target;
        }
        if (l->env > gmax) gmax = l->env;
        if (l->env < 1.0f) { lim++; float r = 1.0f - l->env; if (r > red) red = r; }
        float y = x * l->env;
        if (y > ceiling) y = ceiling;
        if (y < -ceiling) y = -ceiling;
        out[i] = y;
        if (fabsf(y) > pout) pout = fabsf(y);
    }
    if (rep) {
        rep->peak_in = pin; rep->peak_out = pout;
        rep->limited_samples = lim; rep->max_reduction = red;
        rep->max_gain = gmax; rep->end_gain = l->env;
    }
    l->engaged = lim > 0;
    return n;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Limiter a, b;
    float in[BLK], out[BLK];
    Report rep, ra, rb;
    int i;

    for (i = 0; i < BLK; i++) in[i] = 0.5f;
    lim_init(&a, 5.0f, 60.0f, 400.0f, 0.25f, 48000.0f);
    check("quiet_passthrough", lim_process(&a, in, out, BLK, 0.9f, &rep) == BLK && rep.limited_samples == 0);
    check("quiet_peak", fabsf(rep.peak_out - 0.5f) < 1e-4f);

    for (i = 0; i < BLK; i++) in[i] = 2.0f;
    lim_init(&a, 5.0f, 60.0f, 400.0f, 0.25f, 48000.0f);
    lim_process(&a, in, out, BLK, 0.8f, &rep);
    check("loud_limited", rep.limited_samples > 0);
    check("below_ceiling", rep.peak_out <= 0.8f + 1e-4f);
    check("reduction_positive", rep.max_reduction > 0.0f);

    /* hot level envelope picks the fast release */
    lim_init(&a, 5.0f, 60.0f, 400.0f, 0.25f, 48000.0f);
    lim_init(&b, 5.0f, 400.0f, 400.0f, 0.25f, 48000.0f);
    lim_process(&a, in, out, BLK, 0.8f, &rep);   /* in[] still 2.0 */
    lim_process(&b, in, out, BLK, 0.8f, &rep);
    for (i = 0; i < BLK; i++) in[i] = 0.0f;
    lim_process(&a, in, out, BLK, 0.8f, &ra);
    lim_process(&b, in, out, BLK, 0.8f, &rb);
    check("fast_release_when_hot", ra.end_gain > rb.end_gain);

    /* quiet level envelope picks the slow release */
    lim_init(&a, 5.0f, 200.0f, 800.0f, 0.25f, 48000.0f);
    lim_init(&b, 5.0f, 200.0f, 200.0f, 0.25f, 48000.0f);
    for (i = 0; i < BLK; i++) in[i] = 2.0f;
    lim_process(&a, in, out, BLK, 0.8f, &rep);
    lim_process(&b, in, out, BLK, 0.8f, &rep);
    for (i = 0; i < BLK; i++) in[i] = 0.0f;
    for (int k = 0; k < 300; k++) { /* ~0.8 s: level decays below split, env still recovering */
        lim_process(&a, in, out, BLK, 0.8f, &ra);
        lim_process(&b, in, out, BLK, 0.8f, &rb);
    }
    check("slow_release_when_quiet", ra.end_gain < rb.end_gain);

    /* pure impulse: gain dips and recovers, never overshoots above 1 */
    for (i = 0; i < BLK; i++) in[i] = 0.0f;
    in[0] = 1.0f;
    lim_init(&a, 5.0f, 60.0f, 400.0f, 0.25f, 48000.0f);
    lim_process(&a, in, out, BLK, 0.8f, &rep);
    check("impulse_dips", rep.max_reduction > 0.0f);
    check("impulse_no_overshoot", rep.max_gain <= 1.0f + 1e-6f);

    check("bad_ceiling", lim_process(&a, in, out, BLK, 0.0f, &rep) == -1);
    lim_init(&a, 0.0f, 60.0f, 400.0f, 0.25f, 48000.0f);
    check("bad_attack", lim_process(&a, in, out, BLK, 0.8f, &rep) == -1);
    lim_init(&a, 5.0f, 60.0f, 0.0f, 0.25f, 48000.0f);
    check("bad_release", lim_process(&a, in, out, BLK, 0.8f, &rep) == -1);
    lim_init(&a, 5.0f, 60.0f, 400.0f, 0.0f, 48000.0f);
    check("bad_split", lim_process(&a, in, out, BLK, 0.8f, &rep) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
