/* doppler_report: drives src/doppler.c to produce the figures quoted in the
 * module docs. Serves as the source of truth for the CSV/JSON under
 * artifacts/doppler/. Not part of the module API.
 *
 *   ./report <outdir>
 *
 * Outputs:
 *   sweep.csv   listener-time trace: analytic f, phase-reference f, errors
 *   alias.csv   per-sample phase step against the Nyquist edge (low rates)
 *   table.csv   per-profile / per-sample-rate error summary
 *   timing.csv  solve cost for the two speed profiles
 *   summary.json scalars used by the write-up
 */
#define _POSIX_C_SOURCE 199309L
#define DOPPLER_NO_MAIN
#include "../src/doppler.c"

#include <stdlib.h>
#include <time.h>

#define T0 (-0.5)
#define PASS_SECONDS 1.0

typedef struct {
    const char* name;
    double x0, v0, a;
} Profile;

static const Profile PROFILES[] = {
    { "slow", 0.0, 30.0, 0.0 },
    { "fast", 0.0, 120.0, 20.0 },
};

static const int rates[] = { 48000, 8000, 1000, 880 };
#define NRATES (int)(sizeof(rates) / sizeof(rates[0]))

typedef struct {
    double fmin, fmax, chirp_max;
    double e_back_max, e_back_rms, e_center_max;
    double drift;
    double vrmax;
    int aliased;
    int nvalid;
    double ratio_span;
} Stats;

/* One pass over n samples. Arrays may be NULL when only stats are wanted.
 * Samples are classified from the third sample on, the first one for which
 * both a backward and a centered phase difference exist. */
static int collect(const Doppler* m, double fs, int n, Stats* s,
                   double* t_out, double* f_out, double* fref_out,
                   double* eback_out, double* ecenter_out, double* chirp_out,
                   int* alias_out) {
    PhaseRef pr;
    double h = 1.0 / fs;
    int eb_n = 0, ec_n = 0;
    double eb_max = 0.0, ec_max = 0.0, eb_sq = 0.0;
    phase_ref_init(&pr, fs);
    s->fmin = INFINITY; s->fmax = -INFINITY;
    s->chirp_max = 0.0; s->vrmax = 0.0;
    s->e_back_max = NAN; s->e_back_rms = NAN; s->e_center_max = NAN;
    s->aliased = 0; s->nvalid = 0;

    for (int i = 0; i < n; i++) {
        double t = T0 + i * h;
        DopplerPoint p;
        double ch = 0.0;
        if (doppler_at(m, t, &p) != 0) return -1;
        if (doppler_chirp(m, t, &ch) != 0) return -1;
        if (phase_ref_push(&pr, m, t) != 0) return -1;

        s->fmin = fmin(s->fmin, p.f);
        s->fmax = fmax(s->fmax, p.f);
        s->chirp_max = fmax(s->chirp_max, fabs(ch));
        s->vrmax = fmax(s->vrmax, fabs(p.vr));
        if (t_out) t_out[i] = t;
        if (f_out) f_out[i] = p.f;
        if (chirp_out) chirp_out[i] = ch;
        if (fref_out) fref_out[i] = NAN;
        if (eback_out) eback_out[i] = NAN;
        if (ecenter_out) ecenter_out[i] = NAN;
        if (alias_out) alias_out[i] = 0;

        /* Per-sample errors over unaliased samples only: a folded step is not
         * an estimator error, it is the method leaving its domain. */
        if (pr.count >= 2) {
            double eb = pr.f_back - p.f;
            if (fref_out) fref_out[i] = pr.f_back;
            if (eback_out) eback_out[i] = eb;
            if (isfinite(eb)) {
                double a = fabs(eb);
                if (a > eb_max) eb_max = a;
                eb_sq += eb * eb; eb_n++;
            }
        }
        if (pr.count >= 3) {
            /* centered estimate sits at t_{i-1}; compare with f there */
            DopplerPoint q;
            double fprev = doppler_at(m, t - h, &q) == 0 ? q.f : p.f;
            double ec = pr.f_center - fprev;
            if (ecenter_out) ecenter_out[i] = ec;
            if (isfinite(ec)) {
                double a = fabs(ec);
                if (a > ec_max) ec_max = a;
                ec_n++;
            }
            if (pr.aliased) s->aliased++;
            else s->nvalid++;
            if (alias_out) alias_out[i] = pr.aliased;
        }
    }

    if (eb_n > 0) { s->e_back_max = eb_max; s->e_back_rms = sqrt(eb_sq / eb_n); }
    if (ec_n > 0) s->e_center_max = ec_max;
    s->ratio_span = s->fmax / s->fmin;
    return 0;
}

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void measure_timing(FILE* fp, const Profile* pf, int fs, int n, int repeats) {
    Doppler m;
    double sink = 0.0;
    doppler_init(&m, 440.0, DOPPLER_C, pf->x0, pf->v0, pf->a, 5.0);
    double t0w = now_sec();
    for (int r = 0; r < repeats; r++) {
        for (int i = 0; i < n; i++) {
            DopplerPoint p;
            if (doppler_at(&m, T0 + i / (double)fs, &p) == 0) sink += p.f;
        }
    }
    double total = now_sec() - t0w;
    double per = total / (double)(n * repeats);
    fprintf(fp, "%s,%d,%d,%d,%.6f,%.6f,%.3f\n",
            pf->name, fs, n, repeats, total, per * 1e6, per * 1e9);
    if (sink == 12345.6789) fprintf(stderr, "unreachable %g\n", sink);
}

int main(int argc, char** argv) {
    const char* dir = argc > 1 ? argv[1] : "artifacts/doppler";
    char path[1024];
    FILE* fp;
    Doppler m;

    int n48 = 48000;
    double* t = malloc((size_t)n48 * sizeof(double));
    double* f = malloc((size_t)n48 * sizeof(double));
    double* fr = malloc((size_t)n48 * sizeof(double));
    double* eb = malloc((size_t)n48 * sizeof(double));
    double* ec = malloc((size_t)n48 * sizeof(double));
    double* ch = malloc((size_t)n48 * sizeof(double));
    int* al = malloc((size_t)n48 * sizeof(int));
    if (!t || !f || !fr || !eb || !ec || !ch || !al) return 1;

    /* ---- sweep: listener-time trace at 48 kHz, one row per millisecond ---- */
    snprintf(path, sizeof(path), "%s/sweep.csv", dir);
    fp = fopen(path, "w");
    if (!fp) { perror(path); return 1; }
    fprintf(fp, "profile,t,f_analytic,f_phase_back,f_err_back,f_err_center,aliased\n");
    for (int pi = 0; pi < 2; pi++) {
        const Profile* pf = &PROFILES[pi];
        Stats s;
        if (doppler_init(&m, 440.0, DOPPLER_C, pf->x0, pf->v0, pf->a, 5.0) != 0) return 1;
        if (collect(&m, 48000.0, n48, &s, t, f, fr, eb, ec, ch, al) != 0) return 1;
        for (int i = 0; i < n48; i += 48)        /* every 1 ms */
            fprintf(fp, "%s,%.9f,%.9f,%.9f,%.9f,%.9f,%d\n",
                    pf->name, t[i], f[i], fr[i], eb[i], ec[i], al[i]);
    }
    fclose(fp);

    /* ---- alias: per-sample phase step vs the Nyquist edge, low rates ---- */
    snprintf(path, sizeof(path), "%s/alias.csv", dir);
    fp = fopen(path, "w");
    if (!fp) { perror(path); return 1; }
    fprintf(fp, "profile,fs,t,f_analytic,phase_step_pi,f_over_nyquist,aliased\n");
    {
        int lowrates[2] = { 1000, 880 };
        for (int pi = 0; pi < 2; pi++) {
            const Profile* pf = &PROFILES[pi];
            for (int k = 0; k < 2; k++) {
                int fs = lowrates[k], n = fs;
                if (doppler_init(&m, 440.0, DOPPLER_C, pf->x0, pf->v0, pf->a, 5.0) != 0) return 1;
                PhaseRef q;
                phase_ref_init(&q, fs);
                for (int i = 0; i < n; i++) {
                    double tt = T0 + i / (double)fs;
                    if (phase_ref_push(&q, &m, tt) != 0) return 1;
                    if (q.count < 3) continue;   /* need one previous step */
                    DopplerPoint p;
                    if (doppler_at(&m, tt, &p) != 0) return 1;
                    fprintf(fp, "%s,%d,%.9f,%.6f,%.6f,%.6f,%d\n",
                            pf->name, fs, tt, p.f, fabs(q.step1) / M_PI,
                            p.f / (fs / 2.0), q.aliased);
                }
            }
        }
    }
    fclose(fp);

    /* ---- table: error summary across profiles and rates ---- */
    snprintf(path, sizeof(path), "%s/table.csv", dir);
    fp = fopen(path, "w");
    if (!fp) { perror(path); return 1; }
    fprintf(fp, "profile,fs,t0,n,f_min,f_max,abs_vr_max,chirp_max_hz_s,"
                "err_back_max,err_back_rms,err_center_max,phase_drift,"
                "h_over_2_chirp,ratio,valid_samples,aliased_samples\n");
    for (int pi = 0; pi < 2; pi++) {
        const Profile* pf = &PROFILES[pi];
        for (int r = 0; r < NRATES; r++) {
            int fs = rates[r], n = fs;
            Stats s;
            double bound;
            if (doppler_init(&m, 440.0, DOPPLER_C, pf->x0, pf->v0, pf->a, 5.0) != 0) return 1;
            if (collect(&m, (double)fs, n, &s, t, f, fr, eb, ec, ch, al) != 0) return 1;
            s.drift = phase_drift(&m, (double)fs, T0, n);
            bound = 0.5 / (double)fs * s.chirp_max;
            fprintf(fp, "%s,%d,%.1f,%d,%.6f,%.6f,%.4f,%.4f,"
                        "%.9g,%.9g,%.9g,%.9g,%.9g,%.9f,%d,%d\n",
                    pf->name, fs, T0, n, s.fmin, s.fmax, s.vrmax, s.chirp_max,
                    s.e_back_max, s.e_back_rms, s.e_center_max, s.drift,
                    bound, s.ratio_span, s.nvalid, s.aliased);
        }
    }
    fclose(fp);

    /* ---- summary: scalars quoted in the write-up ---- */
    Stats slow_s, fast_s;
    if (doppler_init(&m, 440.0, DOPPLER_C, PROFILES[0].x0, PROFILES[0].v0,
                     PROFILES[0].a, 5.0) != 0) return 1;
    if (collect(&m, 48000.0, n48, &slow_s, t, f, fr, eb, ec, ch, al) != 0) return 1;
    slow_s.drift = phase_drift(&m, 48000.0, T0, n48);

    if (doppler_init(&m, 440.0, DOPPLER_C, PROFILES[1].x0, PROFILES[1].v0,
                     PROFILES[1].a, 5.0) != 0) return 1;
    if (collect(&m, 48000.0, n48, &fast_s, t, f, fr, eb, ec, ch, al) != 0) return 1;
    fast_s.drift = phase_drift(&m, 48000.0, T0, n48);

    snprintf(path, sizeof(path), "%s/summary.json", dir);
    fp = fopen(path, "w");
    if (!fp) { perror(path); return 1; }
    fprintf(fp, "{\n  \"f0_hz\": 440.0,\n  \"c_m_s\": 343.0,\n  \"d_m\": 5.0,\n");
    fprintf(fp, "  \"slow\": {\"f_min\": %.6f, \"f_max\": %.6f, \"drift\": %.9g, "
                "\"err_back_max\": %.9g},\n",
            slow_s.fmin, slow_s.fmax, slow_s.drift, slow_s.e_back_max);
    fprintf(fp, "  \"fast\": {\"f_min\": %.6f, \"f_max\": %.6f, \"drift\": %.9g, "
                "\"err_back_max\": %.9g},\n",
            fast_s.fmin, fast_s.fmax, fast_s.drift, fast_s.e_back_max);
    fprintf(fp, "  \"alias_onset_vr_m_s\": {");
    for (int r = 0; r < NRATES; r++) {
        double fs = (double)rates[r];
        double vr = DOPPLER_C * (2.0 * 440.0 / fs - 1.0);   /* f0*c/(c+vr) = fs/2 */
        fprintf(fp, "%s\"%d\": %.6f", r ? ", " : "", rates[r], vr);
    }
    fprintf(fp, "}\n}\n");
    fclose(fp);

    /* ---- timing: per-sample analytic solve for the two profiles ---- */
    snprintf(path, sizeof(path), "%s/timing.csv", dir);
    fp = fopen(path, "w");
    if (!fp) { perror(path); return 1; }
    fprintf(fp, "profile,fs,n,repeats,total_s,us_per_sample,ns_per_sample\n");
    for (int pi = 0; pi < 2; pi++)
        measure_timing(fp, &PROFILES[pi], 48000, 48000, 20);
    fclose(fp);

    free(t); free(f); free(fr); free(eb); free(ec); free(ch); free(al);
    return 0;
}
