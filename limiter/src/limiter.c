/* limiter: peak limiter with attack/release envelope and gain reduction report.
 *
 * The envelope is a one-pole recursion driven by the per-sample target gain
 * (ceiling / |x|); attack and release are selected by whether the target lies
 * below or above the current envelope.  For a constant target the recursion
 * has the closed form
 *
 *     env[n] = target + (env[0] - target) * exp(-n / (tau * sr))
 *
 * `make test` prints four acceptance artifacts:
 *   [envelope_curve]        env(t) over a step-in / release-out signal
 *   [gain_reduction_table]  per-window peaks and max gain reduction in dB
 *   [time_constants]        measured vs analytic settling times (耗时)
 *   [bypass_identity]       bypassed output is byte-identical to the input
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define BLK 128

typedef struct {
    float attack_ms;
    float release_ms;
    float sample_rate;
    float env;
    int engaged;
    int bypass;     /* 1: pass samples through untouched (skip limiting) */
} Limiter;

typedef struct { float peak_in; float peak_out; int limited_samples; float max_reduction; } Report;

static void lim_init(Limiter* l, float attack_ms, float release_ms, float sr) {
    memset(l, 0, sizeof(*l));
    l->attack_ms = attack_ms;
    l->release_ms = release_ms;
    l->sample_rate = sr;
    l->env = 1.0f;
}

/* one-pole coefficient from a time constant: env reaches 1-1/e after ms */
static float coef(float ms, float sr) {
    if (ms <= 0.0f || sr <= 0.0f) return 0.0f;
    return expf(-1.0f / (ms * 0.001f * sr));
}

/* one envelope update; the same recursion used by lim_process */
static float env_update(float env, float target, float ca, float cr) {
    float c = target < env ? ca : cr;
    return c * env + (1.0f - c) * target;
}

/* lim_process plus an optional per-sample envelope trace (may be NULL). */
static int lim_process_ex(Limiter* l, const float* in, float* out, int n,
                          float ceiling, Report* rep, float* env_trace) {
    if (l == NULL || in == NULL || out == NULL || n <= 0) return -1;
    if (!(ceiling > 0.0f && ceiling <= 1.0f)) return -1;
    if (!(l->attack_ms > 0.0f) || !(l->release_ms > 0.0f) || !(l->sample_rate > 0.0f)) return -1;

    float pin = 0.0f, pout = 0.0f, red = 0.0f;
    int lim = 0;

    if (l->bypass) {
        /* skip limiting entirely: no gain multiply, no clipping, no state update */
        if (out != in) memmove(out, in, (size_t)n * sizeof(float));
        for (int i = 0; i < n; i++) {
            float a = fabsf(in[i]);
            if (a > pin) pin = a;
            if (env_trace) env_trace[i] = 1.0f;
        }
        if (rep) {
            rep->peak_in = pin;
            rep->peak_out = pin;
            rep->limited_samples = 0;
            rep->max_reduction = 0.0f;
        }
        l->engaged = 0;
        return n;
    }

    float ca = coef(l->attack_ms, l->sample_rate);
    float cr = coef(l->release_ms, l->sample_rate);

    for (int i = 0; i < n; i++) {
        float x = in[i];
        float ax = fabsf(x);
        if (ax > pin) pin = ax;
        float target = ax > ceiling ? ceiling / ax : 1.0f;
        l->env = env_update(l->env, target, ca, cr);
        if (l->env < 1.0f) {
            lim++;
            float r = 1.0f - l->env;
            if (r > red) red = r;
        }
        float y = x * l->env;
        if (y > ceiling) y = ceiling;
        if (y < -ceiling) y = -ceiling;
        out[i] = y;
        if (env_trace) env_trace[i] = l->env;
        if (fabsf(y) > pout) pout = fabsf(y);
    }
    if (rep) {
        rep->peak_in = pin;
        rep->peak_out = pout;
        rep->limited_samples = lim;
        rep->max_reduction = red;
    }
    l->engaged = lim > 0;
    return n;
}

static int lim_process(Limiter* l, const float* in, float* out, int n, float ceiling, Report* rep) {
    return lim_process_ex(l, in, out, n, ceiling, rep, NULL);
}

/* ---- measurement / reporting ------------------------------------------- */

#define SR          48000.0f
#define NSAMP       19200                        /* 400 ms */
#define ATTACK_MS   5.0f
#define RELEASE_MS  50.0f
#define CEILING     0.8f
#define STEP_AT     480                          /* 10 ms */
#define RELEASE_AT  2400                         /* 50 ms */

static float g_in[NSAMP];
static float g_out[NSAMP];
static float g_env[NSAMP];
static float g_env_analytic[NSAMP];

static void fill_step_signal(float* x, int n) {
    for (int i = 0; i < n; i++) {
        float amp = (i >= STEP_AT && i < RELEASE_AT) ? 2.0f : 0.2f;
        x[i] = (i & 1) ? -amp : amp;             /* square wave, |x| is the step */
    }
}

static void print_curve_point(const float* env, const float* ana, int i) {
    double db = 20.0 * log10((double)env[i]);
    if (db > -0.0005) db = 0.0;                  /* avoid "-0.000" */
    double err = fabs((double)env[i] - (double)ana[i]);
    printf("%8.2f,%.6f,%.3f,%.6f,%.2e\n", 1000.0 * i / SR, env[i], db, ana[i], err);
}

static void report_envelope_curve(const float* env, const float* ana) {
    static const int times_ms[] = {
        0, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 25,
        30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 120, 140,
        160, 200, 250, 300, 350, 399
    };
    printf("[envelope_curve] sr=%.0f attack_ms=%.1f release_ms=%.1f ceiling=%.1f "
           "signal=square 0.2 -> 2.0 @10ms -> 0.2 @50ms\n",
           SR, ATTACK_MS, RELEASE_MS, CEILING);
    printf("t_ms,env,gain_db,analytic_env,error\n");
    for (size_t k = 0; k < sizeof(times_ms) / sizeof(times_ms[0]); k++) {
        int i = (int)(times_ms[k] * SR / 1000.0f + 0.5f);
        if (i >= NSAMP) i = NSAMP - 1;
        print_curve_point(env, ana, i);
    }
}

static void report_gain_reduction(const float* in, const float* out, const float* env) {
    const int win = (int)(0.020f * SR);          /* 20 ms windows */
    printf("[gain_reduction_table] window_ms=20\n");
    printf("t_start_ms,t_end_ms,in_peak,out_peak,env_min,gr_max_db,limited_samples\n");
    for (int s = 0; s < NSAMP; s += win) {
        int e = s + win;
        if (e > NSAMP) e = NSAMP;
        float ip = 0.0f, op = 0.0f, emin = 1.0f;
        int lim = 0;
        for (int i = s; i < e; i++) {
            float a = fabsf(in[i]);
            float b = fabsf(out[i]);
            if (a > ip) ip = a;
            if (b > op) op = b;
            if (env[i] < emin) emin = env[i];
            if (env[i] < 1.0f - 1e-6f) lim++;
        }
        double gr = emin < 1.0f - 1e-6f ? -20.0 * log10((double)emin) : 0.0;
        printf("%.2f,%.2f,%.4f,%.4f,%.6f,%.3f,%d\n",
               1000.0 * s / SR, 1000.0 * e / SR, ip, op, emin, gr, lim);
    }
}

/* first sample time (ms) at which env has covered `frac` of the change
 * env0 -> target; returns -1 if it never does inside the segment */
static double settle_ms(const float* env, int start, int count,
                        float env0, float target, double frac) {
    double level = (double)env0 + ((double)target - (double)env0) * frac;
    for (int i = 0; i < count; i++) {
        float e = env[start + i];
        int reached = target < env0 ? (e <= level) : (e >= level);
        if (reached) return 1000.0 * (i + 1) / SR;
    }
    return -1.0;
}

static double analytic_settle_ms(double tau_ms, double frac) {
    return -tau_ms * log(1.0 - frac);
}

static void report_time_constants(const float* env, double* t63_a, double* t99_a,
                                  double* t63_r, double* t99_r) {
    const double f63 = 1.0 - 1.0 / 2.718281828459045;
    const double f90 = 0.90;
    const double f99 = 0.99;
    float rel_env0 = env[RELEASE_AT - 1];

    double a63 = settle_ms(env, STEP_AT, RELEASE_AT - STEP_AT, 1.0f, 0.4f, f63);
    double a90 = settle_ms(env, STEP_AT, RELEASE_AT - STEP_AT, 1.0f, 0.4f, f90);
    double a99 = settle_ms(env, STEP_AT, RELEASE_AT - STEP_AT, 1.0f, 0.4f, f99);
    double r63 = settle_ms(env, RELEASE_AT, NSAMP - RELEASE_AT, rel_env0, 1.0f, f63);
    double r90 = settle_ms(env, RELEASE_AT, NSAMP - RELEASE_AT, rel_env0, 1.0f, f90);
    double r99 = settle_ms(env, RELEASE_AT, NSAMP - RELEASE_AT, rel_env0, 1.0f, f99);

    *t63_a = a63;
    *t99_a = a99;
    *t63_r = r63;
    *t99_r = r99;

    printf("[time_constants] settling 耗时, measured from the envelope trace\n");
    printf("constant,tau_ms,t63_meas_ms,t63_analytic_ms,t90_meas_ms,t90_analytic_ms,"
           "t99_meas_ms,t99_analytic_ms\n");
    printf("attack,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
           (double)ATTACK_MS, a63, analytic_settle_ms(ATTACK_MS, f63),
           a90, analytic_settle_ms(ATTACK_MS, f90),
           a99, analytic_settle_ms(ATTACK_MS, f99));
    printf("release,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
           (double)RELEASE_MS, r63, analytic_settle_ms(RELEASE_MS, f63),
           r90, analytic_settle_ms(RELEASE_MS, f90),
           r99, analytic_settle_ms(RELEASE_MS, f99));
}

static unsigned long long fnv1a(const void* data, size_t bytes) {
    const unsigned char* p = (const unsigned char*)data;
    unsigned long long h = 1469598103934665603ULL;
    for (size_t i = 0; i < bytes; i++) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

/* report max |env - closed form| for the constant-target attack / release
 * segments of the step signal */
static void analytic_error(const float* env, double* err_attack, double* err_release) {
    double ea = 0.0, er = 0.0;
    double tau_a = (double)ATTACK_MS * 0.001 * (double)SR;
    double tau_r = (double)RELEASE_MS * 0.001 * (double)SR;
    for (int i = STEP_AT; i < RELEASE_AT; i++) {
        double n = i - STEP_AT + 1;
        double ana = 0.4 + (1.0 - 0.4) * exp(-n / tau_a);
        double e = fabs((double)env[i] - ana);
        if (e > ea) ea = e;
    }
    double env0 = (double)env[RELEASE_AT - 1];
    for (int i = RELEASE_AT; i < NSAMP; i++) {
        double n = i - RELEASE_AT + 1;
        double ana = 1.0 + (env0 - 1.0) * exp(-n / tau_r);
        double e = fabs((double)env[i] - ana);
        if (e > er) er = e;
    }
    *err_attack = ea;
    *err_release = er;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Limiter l;
    Report rep;

    /* ---- legacy self tests ---- */
    float blk_in[BLK], blk_out[BLK];
    for (int i = 0; i < BLK; i++) blk_in[i] = 0.5f;
    lim_init(&l, 5.0f, 50.0f, 48000.0f);
    check("quiet_passthrough", lim_process(&l, blk_in, blk_out, BLK, 0.9f, &rep) == BLK && rep.limited_samples == 0);
    check("quiet_peak", fabsf(rep.peak_out - 0.5f) < 1e-4f);

    for (int i = 0; i < BLK; i++) blk_in[i] = 2.0f;
    lim_init(&l, 5.0f, 50.0f, 48000.0f);
    lim_process(&l, blk_in, blk_out, BLK, 0.8f, &rep);
    check("loud_limited", rep.limited_samples > 0);
    check("below_ceiling", rep.peak_out <= 0.8f + 1e-4f);
    check("reduction_positive", rep.max_reduction > 0.0f);

    check("bad_ceiling", lim_process(&l, blk_in, blk_out, BLK, 0.0f, &rep) == -1);
    lim_init(&l, 0.0f, 50.0f, 48000.0f);
    check("bad_attack", lim_process(&l, blk_in, blk_out, BLK, 0.8f, &rep) == -1);

    /* ---- step response: envelope trace, tables, analytic comparison ---- */
    fill_step_signal(g_in, NSAMP);
    lim_init(&l, ATTACK_MS, RELEASE_MS, SR);
    lim_process_ex(&l, g_in, g_out, NSAMP, CEILING, &rep, g_env);

    double err_a = 0.0, err_r = 0.0;
    analytic_error(g_env, &err_a, &err_r);

    /* monotonicity: attack never rises, release never falls */
    int mono = 1;
    for (int i = STEP_AT; i < RELEASE_AT; i++) if (g_env[i] > g_env[i - 1] + 1e-7f) mono = 0;
    for (int i = RELEASE_AT; i < NSAMP; i++) if (g_env[i] < g_env[i - 1] - 1e-7f) mono = 0;
    check("envelope_monotonic", mono);
    check("envelope_attack_analytic", err_a < 1e-4);
    check("envelope_release_analytic", err_r < 1e-4);

    /* ---- bypass: bytes must match the input exactly ---- */
    static float g_bypass[NSAMP];
    for (int i = 0; i < NSAMP; i++) g_bypass[i] = 123.0f;
    lim_init(&l, ATTACK_MS, RELEASE_MS, SR);
    l.bypass = 1;
    int rc = lim_process(&l, g_in, g_bypass, NSAMP, CEILING, &rep);
    int bypass_ok = rc == NSAMP && memcmp(g_in, g_bypass, sizeof(g_in)) == 0;

    /* quiet input through the *active* path must also be bit-exact */
    static float quiet_in[NSAMP], quiet_out[NSAMP];
    for (int i = 0; i < NSAMP; i++) quiet_in[i] = (i & 1) ? -0.2f : 0.2f;
    lim_init(&l, ATTACK_MS, RELEASE_MS, SR);
    int rcq = lim_process(&l, quiet_in, quiet_out, NSAMP, CEILING, &rep);
    int quiet_ok = rcq == NSAMP && memcmp(quiet_in, quiet_out, sizeof(quiet_in)) == 0;

    /* ---- print the four acceptance artifacts ---- */
    double t63_a = 0.0, t99_a = 0.0, t63_r = 0.0, t99_r = 0.0;

    /* closed form: constant target 0.4 during attack, 1.0 during release */
    double tau_a = (double)ATTACK_MS * 0.001 * (double)SR;
    double tau_r = (double)RELEASE_MS * 0.001 * (double)SR;
    double rel0 = (double)g_env[RELEASE_AT - 1];
    for (int i = 0; i < NSAMP; i++) {
        if (i < STEP_AT) g_env_analytic[i] = 1.0f;
        else if (i < RELEASE_AT)
            g_env_analytic[i] = (float)(0.4 + 0.6 * exp(-(i - STEP_AT + 1) / tau_a));
        else
            g_env_analytic[i] = (float)(1.0 + (rel0 - 1.0) * exp(-(i - RELEASE_AT + 1) / tau_r));
    }

    report_envelope_curve(g_env, g_env_analytic);
    report_gain_reduction(g_in, g_out, g_env);
    report_time_constants(g_env, &t63_a, &t99_a, &t63_r, &t99_r);

    printf("[bypass_identity]\n");
    printf("samples=%d,bytes=%zu,fnv1a_in=%016llx,fnv1a_bypass=%016llx,memcmp=%s\n",
           NSAMP, sizeof(g_in),
           fnv1a(g_in, sizeof(g_in)), fnv1a(g_bypass, sizeof(g_in)),
           bypass_ok ? "identical" : "different");
    printf("quiet_active_memcmp=%s\n", quiet_ok ? "identical" : "different");

    /* ---- checks that consume the measured artifacts ---- */
    check("attack_t63_matches", fabs(t63_a - analytic_settle_ms(ATTACK_MS, 1.0 - 1.0 / 2.718281828459045)) < 0.05);
    check("attack_t99_matches", fabs(t99_a - analytic_settle_ms(ATTACK_MS, 0.99)) < 0.10);
    check("release_t63_matches", fabs(t63_r - analytic_settle_ms(RELEASE_MS, 1.0 - 1.0 / 2.718281828459045)) < 0.05);
    check("release_t99_matches", fabs(t99_r - analytic_settle_ms(RELEASE_MS, 0.99)) < 0.10);
    check("bypass_byte_identical", bypass_ok);
    check("quiet_byte_identical", quiet_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
