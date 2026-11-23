/* meter: frame-averaged power detector; variance vs latency trade-off bench. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define MAXN 256
#define FRAME 256
#define SR 48000.0
#define TF_MS ((double)FRAME / SR * 1000.0)   /* 5.3333... ms */
#define PI 3.14159265358979323846

#define DB_PER_LN  (10.0 / log(10.0))                    /* 4.34294...   */
#define K_AGE  (DB_PER_LN * DB_PER_LN / SR * 1000.0)     /* 0.39294 dB^2*ms */
#define K_WIN  (2.0 * K_AGE)                              /* 0.78588 dB^2*ms */

typedef struct {
    int frame_len;
    int n_frames;
    int idx;
    int fill;
    double sum;
    double hist[MAXN];
} Meter;

static int meter_init(Meter* m, int frame_len, int n_frames) {
    if (m == NULL || frame_len <= 0 || n_frames <= 0 || n_frames > MAXN) return -1;
    memset(m, 0, sizeof(*m));
    m->frame_len = frame_len;
    m->n_frames = n_frames;
    return 0;
}

static double frame_power(const float* x, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += (double)x[i] * (double)x[i];
    return s / (double)n;
}

/* push one frame's mean square, return the mean over the last n_frames
 * (partial window during the first n_frames-1 pushes). */
static double meter_push(Meter* m, double p) {
    if (m->fill == m->n_frames) m->sum -= m->hist[m->idx];
    else m->fill++;
    m->hist[m->idx] = p;
    m->sum += p;
    m->idx = (m->idx + 1) % m->n_frames;
    return m->sum / (double)m->fill;
}

static double to_db(double p) { return 10.0 * log10(p + 1e-30); }

/* ---------- self test ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

typedef struct { uint64_t s; } Rng;
static void rng_seed(Rng* r, uint64_t s) { r->s = s ? s : 1; }
static uint64_t rng_next(Rng* r) {
    r->s ^= r->s << 13; r->s ^= r->s >> 7; r->s ^= r->s << 17;
    return r->s;
}
static double rng_u(Rng* r) { return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0); }
static double rng_n(Rng* r) {
    double u = rng_u(r), v = rng_u(r);
    if (u < 1e-12) u = 1e-12;
    return sqrt(-2.0 * log(u)) * cos(2.0 * PI * v);
}

/* lowest k with k/N >= frac, closed form for the boxcar step response */
static int ceil_frames(int N, double frac) { return (int)ceil((double)N * frac); }

/* first frame index after a step where output is within 0.5 dB of final */
static int step_settle_frames(int N) {
    Meter m = {0};
    if (meter_init(&m, FRAME, N) != 0) return -1;
    for (int i = 0; i < 100; i++) meter_push(&m, 1.0);
    for (int k = 1; k <= 1000; k++) {
        double y = meter_push(&m, 100.0);
        if (to_db(y) >= to_db(100.0) - 0.5) return k;
    }
    return -1;
}

/* centroid of the boxcar impulse response, in frames */
static double group_delay_frames(int N) {
    Meter m = {0};
    if (meter_init(&m, FRAME, N) != 0) return -1.0;
    for (int i = 0; i < N - 1; i++) meter_push(&m, 0.0);
    double num = 0.0, den = 0.0;
    for (int k = 0; k < N; k++) {
        double y = meter_push(&m, k == 0 ? 1.0 : 0.0);
        num += (double)k * y;
        den += y;
    }
    return num / den;
}

/* mean event-to-report age: impulse centroid plus the half-frame average
 * age of an event whose position inside the frame is uniform. */
static double mean_age_frames(int N) { return group_delay_frames(N) + 0.5; }

static void selftest(void) {
    Meter m = {0};
    float dc[FRAME];

    check("init_bad_n", meter_init(&m, FRAME, 0) == -1);
    check("init_bad_frame", meter_init(&m, 0, 8) == -1);
    check("init_too_many", meter_init(&m, FRAME, MAXN + 1) == -1);

    for (int i = 0; i < FRAME; i++) dc[i] = 0.5f;
    check("frame_power_dc", fabs(frame_power(dc, FRAME) - 0.25) < 1e-12);

    check("init_ok", meter_init(&m, FRAME, 1) == 0);
    for (int i = 0; i < 10; i++) meter_push(&m, 0.25);
    check("n1_passthrough", fabs(meter_push(&m, 0.5) - 0.5) < 1e-12);

    check("init_ok8", meter_init(&m, FRAME, 8) == 0);
    for (int i = 0; i < 8; i++) meter_push(&m, 0.25);
    check("dc_hold", fabs(meter_push(&m, 0.25) - 0.25) < 1e-12);

    /* closed form for the step: t_settle = ceil(0.8902*N) frames */
    check("step_n1", step_settle_frames(1) == ceil_frames(1, 0.8902));
    check("step_n8", step_settle_frames(8) == ceil_frames(8, 0.8902));
    check("step_n64", step_settle_frames(64) == ceil_frames(64, 0.8902));

    check("gd_n1", fabs(group_delay_frames(1) - 0.0) < 1e-12);
    check("gd_n8", fabs(group_delay_frames(8) - 3.5) < 1e-12);
    check("gd_n64", fabs(group_delay_frames(64) - 31.5) < 1e-12);
    check("age_n8", fabs(mean_age_frames(8) - 4.0) < 1e-12);
    check("age_n64", fabs(mean_age_frames(64) - 32.0) < 1e-12);

    /* dB pipeline: lognormal frame power, sigma_log = 0.30,
     * sd of 10*log10(P) must be DB_PER_LN*0.30 = 1.303 dB. */
    {
        enum { M1 = 1 << 20 };
        Rng r; rng_seed(&r, 42);
        Meter mm = {0};
        meter_init(&mm, FRAME, 1);
        double acc = 0.0, acc2 = 0.0;
        for (int i = 0; i < M1; i++) {
            double d = to_db(exp(0.30 * rng_n(&r)));
            acc += d; acc2 += d * d;
        }
        double mu = acc / M1;
        double sd = sqrt((acc2 - M1 * mu * mu) / (M1 - 1));
        check("db_scale_lognormal", fabs(sd - DB_PER_LN * 0.30) < 0.01 * DB_PER_LN * 0.30);
    }

    /* 1/sqrt(N) law on iid frame powers, sampled on disjoint windows.
     * sigma_log = 0.02 keeps the log-domain second-order term below 0.01%. */
    {
        enum { M = 1 << 21, NSEED = 4 };
        const double SIG = 0.02;
        double sd1 = 0.0;
        for (int N = 1; N <= 64; N *= 8) {
            double sum = 0.0;
            for (int s = 0; s < NSEED; s++) {
                Rng r; rng_seed(&r, 1000 + s);
                Meter mm = {0};
                if (meter_init(&mm, FRAME, N) != 0) { check("init_bench", 0); return; }
                double acc = 0.0, acc2 = 0.0;
                long c = 0;
                for (int i = 0; i < M; i++) {
                    double y = meter_push(&mm, exp(SIG * rng_n(&r)));
                    if ((i + 1) % N == 0) {
                        double d = to_db(y);
                        acc += d; acc2 += d * d; c++;
                    }
                }
                double mu = acc / c;
                sum += sqrt((acc2 - c * mu * mu) / (c - 1));
            }
            double sd = sum / NSEED;
            if (N == 1) sd1 = sd;
            char name[64];
            snprintf(name, sizeof(name), "var_law_n%d", N);
            /* 1/sqrt(N) predicts sd*sqrt(N) == sd1 (2-sigma of the estimator) */
            check(name, fabs(sd * sqrt((double)N) - sd1) < 0.006 * sd1);
        }
    }

    printf("all_ok=%d\n", failures == 0);
}

/* ---------- bench ---------- */

#define SECS 600.0
#define STEPS 400
#define WARM 400
#define NSEED 4
#define SEED0 0x1234567ULL

typedef struct { double b0, b1, b2; } Pink;
static double pink_step(Pink* p, double w) {
    p->b0 = 0.99765 * p->b0 + w * 0.0990460;
    p->b1 = 0.96300 * p->b1 + w * 0.2965164;
    p->b2 = 0.57000 * p->b2 + w * 1.0526913;
    return p->b0 + p->b1 + p->b2 + w * 0.1848;
}

/* Steady-state sd of the dB output, estimated on disjoint windows.
 * Waveshaper: white = iid gaussian; pink = Kellet 1/f (corner ~18 Hz). */
static int variance_run(int pink, int N, uint64_t seed, double* sd_db) {
    double scale = 1.0;
    Rng r; Pink pf = {0.0, 0.0, 0.0};

    if (pink) {
        rng_seed(&r, seed);
        double s = 0.0;
        int n = 0;
        for (int i = 0; i < (int)SR; i++) {
            double v = pink_step(&pf, rng_n(&r));
            if (i >= 4800) { s += v * v; n++; }
        }
        if (!(s > 0.0)) return -1;
        scale = 0.25 / sqrt(s / n);
        pf.b0 = pf.b1 = pf.b2 = 0.0;
    }
    rng_seed(&r, seed);

    Meter m = {0};
    if (meter_init(&m, FRAME, N) != 0) return -1;
    int nf = (int)(SECS * SR / FRAME);
    float buf[FRAME];
    double acc = 0.0, acc2 = 0.0;
    long c = 0;
    for (int f = 0; f < nf; f++) {
        for (int i = 0; i < FRAME; i++) {
            double w = pink ? pink_step(&pf, rng_n(&r)) : rng_n(&r);
            buf[i] = (float)(w * scale);
        }
        double y = meter_push(&m, frame_power(buf, FRAME));
        if (f >= WARM && (f + 1) % N == 0) {
            double d = to_db(y);
            acc += d; acc2 += d * d; c++;
        }
    }
    double mu = acc / c;
    *sd_db = sqrt((acc2 - c * mu * mu) / (double)(c - 1));
    return 0;
}

/* Exponential smoother on frame power, for the cross-check that the
 * exchange relation is not a boxcar artefact. */
static int exp_variance_run(int N, double* sd_db, double* age_ms) {
    const double alpha = 1.0 / (double)N;
    Rng r; rng_seed(&r, SEED0 + 99);
    int nf = (int)(SECS * SR / FRAME);
    float buf[FRAME];
    double y = 0.0;
    long c = 0;
    double acc = 0.0, acc2 = 0.0;
    for (int f = 0; f < nf; f++) {
        for (int i = 0; i < FRAME; i++) buf[i] = (float)(rng_n(&r) * 0.25);
        double p = frame_power(buf, FRAME);
        y = (1.0 - alpha) * y + alpha * p;
        if (f >= WARM) {
            double d = to_db(y);
            acc += d; acc2 += d * d; c++;
        }
    }
    double mu = acc / c;
    *sd_db = sqrt((acc2 - c * mu * mu) / (double)(c - 1));
    /* one-pole centroid (1-alpha)/alpha frames, plus the half-frame
     * average sampling age that the boxcar age also includes */
    *age_ms = ((1.0 - alpha) / alpha + 0.5) * TF_MS;
    return 0;
}

typedef struct {
    int N;
    double sd, sd_lo, sd_hi, theory, ratio, neff;
    double t10, t90, tset, gd, age, win;
} Point;

/* Power-domain step: constant-amplitude tone (5 whole cycles per frame),
 * amplitude A0 -> A1 at frame STEPS. */
static int step_run(int N, double* t10, double* t90, double* tset, double* gd_ms) {
    Meter m = {0};
    if (meter_init(&m, FRAME, N) != 0) return -1;
    static double yout[2 * STEPS];
    float buf[FRAME];
    double w = 2.0 * PI * 937.5 / SR;
    for (int f = 0; f < 2 * STEPS; f++) {
        double amp = (f < STEPS ? 0.05 : 0.4) * sqrt(2.0);
        for (int i = 0; i < FRAME; i++) {
            int idx = f * FRAME + i;
            buf[i] = (float)(amp * sin(w * (double)idx));
        }
        yout[f] = meter_push(&m, frame_power(buf, FRAME));
    }
    double finalp = yout[2 * STEPS - 1];
    int k10 = -1, k90 = -1, ks = -1;
    for (int f = STEPS; f < 2 * STEPS; f++) {
        int k = f - STEPS + 1;
        if (k10 < 0 && yout[f] >= 0.1 * finalp) k10 = k;
        if (k90 < 0 && yout[f] >= 0.9 * finalp) k90 = k;
        if (ks < 0 && yout[f] >= finalp * pow(10.0, -0.05)) ks = k;
    }
    if (k10 < 0 || k90 < 0 || ks < 0) return -1;
    *t10 = k10 * TF_MS;
    *t90 = k90 * TF_MS;
    *tset = ks * TF_MS;
    *gd_ms = group_delay_frames(N) * TF_MS;
    return 0;
}

static int point_run(int sig, int N, Point* p) {
    double sds[NSEED];
    for (int s = 0; s < NSEED; s++) {
        if (variance_run(sig, N, SEED0 + 7919u * (uint64_t)s, &sds[s]) != 0) return -1;
    }
    double mean = 0.0;
    for (int s = 0; s < NSEED; s++) mean += sds[s];
    mean /= NSEED;
    double var = 0.0;
    for (int s = 0; s < NSEED; s++) var += (sds[s] - mean) * (sds[s] - mean);
    var /= (NSEED - 1);
    double half = 3.182 * sqrt(var / NSEED);   /* t(0.975, df=3) */
    p->N = N;
    p->sd = mean; p->sd_lo = mean - half; p->sd_hi = mean + half;
    if (step_run(N, &p->t10, &p->t90, &p->tset, &p->gd) != 0) return -1;
    p->age = mean_age_frames(N) * TF_MS;
    p->win = (double)N * TF_MS;
    p->theory = DB_PER_LN * sqrt(2.0 / ((double)FRAME * (double)N));
    p->ratio = p->sd / p->theory;
    return 0;
}

static void table(const char* name, Point* pts, int npts, int white) {
    printf("# signal=%s\n", name);
    printf("N,win_ms,t10_ms,t90_ms,tsettle_ms,gd_ms,age_ms,sigma_db,sigma_lo,sigma_hi,"
           "sigma_theory,ratio,neff,sigma2_age,sigma2_win\n");
    for (int i = 0; i < npts; i++) {
        Point* p = &pts[i];
        /* neff is relative to this signal's own N=1 point */
        double neff = (pts[0].sd / p->sd) * (pts[0].sd / p->sd);
        char tb[32], rb[32];
        if (white) snprintf(tb, sizeof(tb), "%.4f", p->theory), snprintf(rb, sizeof(rb), "%.3f", p->ratio);
        else snprintf(tb, sizeof(tb), "na"), snprintf(rb, sizeof(rb), "na");
        printf("%d,%.3f,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f,%.4f,%.4f,%s,%s,%.2f,%.4f,%.4f\n",
               p->N, p->win, p->t10, p->t90, p->tset, p->gd, p->age,
               p->sd, p->sd_lo, p->sd_hi, tb, rb, neff,
               p->sd * p->sd * p->age, p->sd * p->sd * p->win);
    }
}

static void bench(double budget_ms, double sigma_target) {
    const int ns[3] = { 1, 8, 64 };
    Point w[3], pk[3];
    printf("# frame=%d samples (%.4f ms) sr=%.0f dur=%.0f s seeds=%d budget_settle=%.1f ms sigma_target=%.2f dB\n",
           FRAME, TF_MS, SR, SECS, NSEED, budget_ms, sigma_target);
    printf("# exact boxcar identity: sigma_db^2 * age_ms = (10/ln10)^2/fs = %.5f dB^2*ms\n", K_AGE);
    printf("# information floor for any estimator using a window T: sigma_db >= %.4f/sqrt(T_ms) dB\n",
           sqrt(K_WIN));
    printf("# window needed for sigma_target = %.2f dB: %.2f ms (mean age %.2f ms)\n",
           sigma_target, K_WIN / (sigma_target * sigma_target),
           0.5 * K_WIN / (sigma_target * sigma_target));
    {
        int n_min = 0, n_max = 0;
        for (int N = 1; N <= MAXN; N++) {
            double win = (double)N * TF_MS;
            double sig = sqrt(K_WIN / win);
            if (sig <= sigma_target) { if (!n_min) n_min = N; }
            if (0.8895 * win <= budget_ms) n_max = N;
        }
        printf("# analytic feasible band on white noise: N in [%d, %d] "
               "(sigma target gives N>=%d, %.1f ms settle budget gives N<=%d)\n",
               n_min, n_max, n_min, budget_ms, n_max);
        if (n_min > n_max) printf("#   EMPTY band: spec is infeasible at every frame count.\n");
    }

    for (int sig = 0; sig < 2; sig++) {
        Point* pts = sig ? pk : w;
        for (int i = 0; i < 3; i++) {
            if (point_run(sig, ns[i], &pts[i]) != 0) { printf("# point_run failed\n"); return; }
        }
        table(sig ? "pink" : "white", pts, 3, !sig);
    }

    printf("# verdicts (budget applies to measured 0.5 dB settle time)\n");
    for (int i = 0; i < 3; i++) {
        Point* p = &w[i];
        int lat_ok = p->tset <= budget_ms;
        int var_ok = p->sd_hi <= sigma_target;
        printf("# white N=%d: tsettle %.2f ms %s (budget %.1f), sigma %.4f dB %s (target %.2f) -> %s\n",
               p->N, p->tset, lat_ok ? "OK" : "VIOLATED", budget_ms,
               p->sd, var_ok ? "OK" : "VIOLATED", sigma_target,
               (lat_ok && var_ok) ? "FEASIBLE" : "INFEASIBLE");
        if (!var_ok) {
            double floor_here = sqrt(K_WIN / p->win);
            printf("#   variance floor for a %.2f ms window is %.4f dB (measured %.4f); target %.2f dB needs window >= %.2f ms\n",
                   p->win, floor_here, p->sd, sigma_target,
                   K_WIN / (sigma_target * sigma_target));
            printf("#   no estimator over that window can beat it, so this point cannot reach the target.\n");
        }
        if (!lat_ok) {
            printf("#   mean age %.2f ms and group delay %.2f ms already exceed the %.1f ms budget;\n",
                   p->age, p->gd, budget_ms);
            printf("#   a uniform %d-frame mean cannot report earlier, so this point cannot meet the budget.\n",
                   p->N);
            printf("#   (the measured %.4f dB would require a %.2f ms window by the floor above.)\n",
                   p->sd, K_WIN / (p->sd * p->sd));
        }
    }

    printf("# pink: the same N buys less, because 1/f frames are correlated\n");
    for (int i = 0; i < 3; i++) {
        Point* p = &pk[i];
        printf("#   pink N=%d: sigma %.4f dB, age %.2f ms, sigma^2*age %.2f dB^2*ms (white constant %.2f)\n",
               p->N, p->sd, p->age, p->sd * p->sd * p->age, K_AGE);
    }

    {
        double sd, age;
        if (exp_variance_run(8, &sd, &age) == 0) {
            printf("# cross-check one-pole smoother alpha=1/8 (not a box mean): age %.2f ms, "
                   "sigma %.4f dB, sigma^2*age %.4f dB^2*ms (analytic sigma %.4f dB)\n",
                   age, sd, sd * sd * age, DB_PER_LN * sqrt(2.0 / FRAME * 0.125 / 1.875));
        }
    }
}

/* Sweep a list of frame counts and apply the same verdict, to show where the
 * feasible band starts and ends. */
static void sweep(const int* ns, int npts, double budget_ms, double sigma_target) {
    printf("# sweep budget_settle=%.1f ms sigma_target=%.2f dB\n", budget_ms, sigma_target);
    printf("N,win_ms,tsettle_ms,age_ms,sigma_db,sigma_lo,sigma_hi,lat_ok,var_ok\n");
    for (int i = 0; i < npts; i++) {
        Point p;
        if (point_run(0, ns[i], &p) != 0) { printf("# point_run failed\n"); return; }
        int lat_ok = p.tset <= budget_ms;
        int var_ok = p.sd_hi <= sigma_target;
        printf("%d,%.3f,%.2f,%.2f,%.4f,%.4f,%.4f,%d,%d\n",
               p.N, p.win, p.tset, p.age, p.sd, p.sd_lo, p.sd_hi, lat_ok, var_ok);
    }
}

/* Feasible band as a function of the spec, all analytic. The measured points
 * are checked against these predictions in `bench` and `sweep`. */
static void spec_map(void) {
    const double budgets[] = { 20.0, 30.0, 50.0, 100.0 };
    const double targets[] = { 0.10, 0.15, 0.20, 0.30 };
    printf("# white-noise feasible N band; row = settle budget (ms), col = sigma target (dB)\n");
    printf("budget_ms");
    for (int j = 0; j < 4; j++) printf(",%.2f_dB", targets[j]);
    printf("\n");
    for (int i = 0; i < 4; i++) {
        printf("%.0f", budgets[i]);
        for (int j = 0; j < 4; j++) {
            int n_min = 0, n_max = 0;
            for (int N = 1; N <= MAXN; N++) {
                double win = (double)N * TF_MS;
                if (sqrt(K_WIN / win) <= targets[j] && !n_min) n_min = N;
                if (0.8895 * win <= budgets[i]) n_max = N;
            }
            if (n_min == 0 || n_min > n_max) printf(",none");
            else printf(",%d-%d", n_min, n_max);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    if (argc > 1 && strcmp(argv[1], "map") == 0) { spec_map(); return 0; }
    if (argc > 1 && strcmp(argv[1], "sweep") == 0) {
        int ns[32], npts = 0;
        for (int i = 2; i < argc && npts < 32; i++) ns[npts++] = atoi(argv[i]);
        if (npts == 0) { int d[] = {4, 6, 7, 9, 10, 11, 16}; memcpy(ns, d, sizeof(d)); npts = 7; }
        sweep(ns, npts, 50.0, 0.15);
        return 0;
    }
    if (argc > 1 && strcmp(argv[1], "bench") == 0) {
        double budget = argc > 2 ? atof(argv[2]) : 50.0;
        double sigma = argc > 3 ? atof(argv[3]) : 0.15;
        bench(budget, sigma);
        return 0;
    }
    selftest();
    return failures == 0 ? 0 : 1;
}
