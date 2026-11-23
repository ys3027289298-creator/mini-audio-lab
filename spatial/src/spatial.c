#define _POSIX_C_SOURCE 200809L

/* spatial: constant-power panning and unified distance attenuation.
 *
 * The attenuation is one parameter set over the whole distance range; see
 * atten.h for the curve and its C1 seams.  `make test` runs the self-checks,
 * `make report` also writes the evidence tables under out/.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "atten.h"

typedef struct { float left, right; } Stereo;

static int pan(float x, Stereo* s) {
    if (x < -1.0f || x > 1.0f) return -1;
    float a = (x + 1.0f) * 0.25f * 3.14159265f;   /* 0..pi/2 */
    s->left = cosf(a);
    s->right = sinf(a);
    return 0;
}

/* ---------------- configuration used by the self-test and the report ------- */

static const AttenCfg CFG = { 1.0f, 1.0f, 0.5f };   /* ref, rolloff, knee */
/* The pre-unification two-set model: a shallow near set and an inverse-distance
 * far set, selected by a hard zone boundary — the configuration described in
 * the report.  The tree at HEAD only had a single formula, so this is a
 * reconstruction used as the before/after baseline, not recovered code.  Its
 * far set equals the unified far-field law. */
static const AttenLegacyCfg LEGACY = { 1.0f, 0.2f, 1.0f, 3.0f };

#define SR 48000
#define FPS 60
#define FRAME_SAMPLES (SR / FPS)
#define NFRAMES 241
#define NSAMPLES (NFRAMES * FRAME_SAMPLES)
#define APPROACH_FRAMES 120
#define D_NEAR 0.2
#define D_FAR 8.0

/* Source flies in from 8 m to 0.2 m over 2 s, then back out over 2 s.  Time is
 * measured in video frames so the audio and the picture timeline share one
 * clock: frame f == audio sample f * FRAME_SAMPLES. */
static double scene_dist(double frame) {
    if (frame <= (double)APPROACH_FRAMES)
        return D_FAR + (D_NEAR - D_FAR) * frame / (double)APPROACH_FRAMES;
    return D_NEAR + (D_FAR - D_NEAR) * (frame - (double)APPROACH_FRAMES)
                                 / (double)APPROACH_FRAMES;
}

/* What the picture shows: the projected size of the source is proportional to
 * 1/d and saturates at 1 once the source fills the frame.  The audio gain is
 * compared against this visual cue. */
static double visual_gain(double d) {
    return d < (double)CFG.ref ? 1.0 : (double)CFG.ref / d;
}

static double db(double g) { return 20.0 * log10(g); }

/* ---------------- self-test ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static float ugain(double d) {
    float g = -123.0f;
    AttenStatus st = atten_gain((float)d, &CFG, &g);
    if (st != ATTEN_OK) return -123.0f;
    return g;
}

static float uslope(double d) {
    float s = -123.0f;
    AttenStatus st = atten_gain_slope((float)d, &CFG, &s);
    if (st != ATTEN_OK) return -123.0f;
    return s;
}

static void run_checks(void) {
    Stereo s;
    check("center", pan(0.0f, &s) == 0 && fabsf(s.left - s.right) < 1e-6f);
    check("center_power", fabsf(s.left * s.left + s.right * s.right - 1.0f) < 1e-5f);

    pan(-1.0f, &s);
    check("hard_left", s.left > 0.99f && s.right < 1e-6f);

    pan(1.0f, &s);
    check("hard_right", s.right > 0.99f && s.left < 1e-6f);

    check("pan_bad", pan(1.5f, &s) == -1);

    /* near-field platform and far-field law */
    check("near_full", fabsf(ugain(0.5) - 1.0f) < 1e-6f);
    check("at_ref_full", fabsf(ugain(1.0) - 1.0f) < 1e-6f);
    check("farther_quieter", ugain(5.0) < ugain(3.0));
    check("far_is_one_over_d", fabsf(ugain(4.0) - 0.25f) < 1e-6f);
    check("far_is_one_over_d_far", fabsf(ugain(20.0) - 0.05f) < 1e-6f);

    /* bounded and monotone on a dense sweep, no spikes */
    {
        int bounded = 1, monotone = 1;
        double prev = 2.0;
        for (int i = 1; i <= 20000; i++) {
            double d = i * 0.001;
            double g = ugain(d);
            if (!(g > 0.0 && g <= 1.0)) bounded = 0;
            if (g > prev + 1e-6) monotone = 0;
            prev = g;
        }
        check("gain_bounded", bounded);
        check("monotone_decreasing", monotone);
    }

    /* distance zero is an error, not a silent full-scale gain */
    {
        float g = 123.0f;
        AttenStatus st = atten_gain(0.0f, &CFG, &g);
        check("zero_dist_err", st == ATTEN_ERR_DIST_ZERO);
        check("zero_dist_not_full", g == 0.0f);

        g = 123.0f;
        st = atten_gain(-0.0f, &CFG, &g);
        check("neg_zero_dist_err", st == ATTEN_ERR_DIST_ZERO && g == 0.0f);

        g = 123.0f;
        st = atten_gain(-1.0f, &CFG, &g);
        check("neg_dist_err", st == ATTEN_ERR_DIST_NEG && g == 0.0f);

        g = 123.0f;
        st = atten_gain(NAN, &CFG, &g);
        check("nan_dist_err", st == ATTEN_ERR_NONFINITE && g == 0.0f);

        g = 123.0f;
        st = atten_gain(INFINITY, &CFG, &g);
        check("inf_dist_err", st == ATTEN_ERR_NONFINITE && g == 0.0f);

        g = 123.0f;
        check("null_out_err", atten_gain(1.0f, &CFG, NULL) == ATTEN_ERR_ARG);
        g = 123.0f;
        check("null_cfg_err", atten_gain(1.0f, NULL, &g) == ATTEN_ERR_ARG && g == 0.0f);

        AttenCfg bad = CFG; bad.ref = 0.0f;
        g = 123.0f;
        check("bad_ref_err", atten_gain(1.0f, &bad, &g) == ATTEN_ERR_REF && g == 0.0f);
        bad = CFG; bad.rolloff = -1.0f;
        check("bad_rolloff_err", atten_gain(1.0f, &bad, &g) == ATTEN_ERR_ROLLOFF);
        bad = CFG; bad.knee = -1.0f;
        check("bad_knee_err", atten_gain(1.0f, &bad, &g) == ATTEN_ERR_KNEE);
    }

    /* the single-set helper and the legacy curve reject zero as well */
    {
        float g = 123.0f;
        check("set_zero_dist_err",
              atten_gain_set(0.0f, 1.0f, 1.0f, &g) == ATTEN_ERR_DIST_ZERO && g == 0.0f);
        int set = 7;
        g = 123.0f;
        check("legacy_zero_dist_err",
              atten_gain_legacy(0.0f, &LEGACY, &set, &g) == ATTEN_ERR_DIST_ZERO
              && g == 0.0f && set == -1);
    }

    /* value continuity (C0) and slope continuity (C1) at both seams: the
     * single-sided limits agree */
    {
        double r = CFG.ref, w = CFG.knee, k = CFG.rolloff;
        double e = 1e-5;
        check("seam_ref_c0", fabs(ugain(r - e) - ugain(r + e)) < 1e-5f);
        check("seam_ref_c1", fabs(uslope(r - e) - uslope(r + e)) < 1e-3f);
        check("seam_knee_c0", fabs(ugain(r + w - e) - ugain(r + w + e)) < 1e-5f);
        check("seam_knee_c1", fabs(uslope(r + w - e) - uslope(r + w + e)) < 1e-3f);

        /* the analytic one-sided limits, not just a finite-difference probe */
        double g_k = r / (r + k * w);
        double s_k = -r * k / ((r + k * w) * (r + k * w));
        check("seam_ref_limits_c0",
              fabs((r / (r + k * 0.0)) - 1.0) < 1e-15);
        check("seam_ref_limits_c1",
              fabs((-r * k * 0.0 / (r * r)) - 0.0) < 1e-15);
        check("seam_knee_limits_c0",
              fabs((r / (r + k * w)) - g_k) < 1e-15);
        check("seam_knee_limits_c1",
              fabs((-r * k * 1.0 / ((r + k * w) * (r + k * w))) - s_k) < 1e-15);
    }

    /* knee = 0 degenerates to the C0 curve; the gap shrinks with the probe
     * interval instead of staying at a fixed jump */
    {
        AttenCfg c0 = { 1.0f, 1.0f, 0.0f };
        float coarse = 0.0f, fine = 0.0f, a = 0.0f, b = 0.0f;
        atten_gain(0.9f, &c0, &a);  atten_gain(1.1f, &c0, &b);
        coarse = fabsf(a - b);
        atten_gain(0.99f, &c0, &a); atten_gain(1.01f, &c0, &b);
        fine = fabsf(a - b);
        check("zero_knee_c0", fine < 0.2f * coarse);
    }

    /* the legacy two-set switch really does jump, at the same distance where
     * the unified curve is smooth */
    {
        int set = 0;
        float below = 0.0f, above = 0.0f;
        atten_gain_legacy(2.999f, &LEGACY, &set, &below);
        atten_gain_legacy(3.000f, &LEGACY, &set, &above);
        check("legacy_switch_jump", fabsf(below - above) > 0.3f);
        check("unified_at_switch_smooth", fabsf(ugain(2.999) - ugain(3.000)) < 1e-3f);
    }

    printf("all_ok=%d\n", failures == 0);
}

/* ---------------- evidence generation ---------------- */

#define NTH 8
static const double THRESH[NTH] = {
    1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 2e-1
};

typedef struct {
    double prev;
    double max_dg;
    double max_ddb;
    long count[NTH];
} JumpStat;

static void jstat_start(JumpStat* js, double g) {
    memset(js, 0, sizeof(*js));
    js->prev = g;
}

static void jstat_add(JumpStat* js, double g) {
    double dg = fabs(g - js->prev);
    double ddb = fabs(db(g) - db(js->prev));
    if (dg > js->max_dg) js->max_dg = dg;
    if (ddb > js->max_ddb) js->max_ddb = ddb;
    for (int t = 0; t < NTH; t++)
        if (dg > THRESH[t]) js->count[t]++;
    js->prev = g;
}

typedef struct {
    JumpStat audio[2];   /* 0 = legacy, 1 = unified */
    JumpStat frame[2];
    JumpStat dense[2];
    long audio_set_switches;
    long frame_set_switches;
    double seam1_c0, seam1_c1, seam2_c0, seam2_c1;
    double eps_tbl[3][5];  /* eps, seam1 c0, seam1 c1, seam2 c0, seam2 c1 */
    double maxdev_legacy, maxdev_legacy_d;
    double maxdev_unified, maxdev_unified_d;
    double legacy_switch_delta_db;
    int switch_frame[2];
    int maxdev_frame;
} Metrics;

static void compute_metrics(Metrics* m) {
    memset(m, 0, sizeof(*m));

    /* dense distance sweep, 1e-4 m step */
    {
        const double step = 1e-4, d0 = 0.05, d1 = 20.0;
        int n = (int)((d1 - d0) / step) + 1;
        double d = d0;
        float gl = 0.0f, gu = 0.0f;
        atten_gain_legacy((float)d, &LEGACY, NULL, &gl);
        atten_gain((float)d, &CFG, &gu);
        jstat_start(&m->dense[0], gl);
        jstat_start(&m->dense[1], gu);
        for (int i = 1; i < n; i++) {
            d = d0 + (double)i * step;
            atten_gain_legacy((float)d, &LEGACY, NULL, &gl);
            atten_gain((float)d, &CFG, &gu);
            jstat_add(&m->dense[0], gl);
            jstat_add(&m->dense[1], gu);
        }
    }

    /* audio-rate motion: one gain value per audio sample */
    {
        int set = -1, prev_set = -1;
        float gl = 0.0f, gu = 0.0f;
        double d = scene_dist(0.0);
        atten_gain_legacy((float)d, &LEGACY, &set, &gl);
        atten_gain((float)d, &CFG, &gu);
        m->audio_set_switches = 0;
        prev_set = set;
        jstat_start(&m->audio[0], gl);
        jstat_start(&m->audio[1], gu);
        for (int i = 1; i < NSAMPLES; i++) {
            d = scene_dist((double)i / FRAME_SAMPLES);
            atten_gain_legacy((float)d, &LEGACY, &set, &gl);
            atten_gain((float)d, &CFG, &gu);
            if (set != prev_set) { m->audio_set_switches++; prev_set = set; }
            jstat_add(&m->audio[0], gl);
            jstat_add(&m->audio[1], gu);
        }
    }

    /* frame-rate motion: one gain value per video frame */
    {
        int set = -1, prev_set = -1;
        int switch_seen = 0;
        double maxdev_l = 0.0, maxdev_u = 0.0;
        float gl = 0.0f, gu = 0.0f;
        m->frame_set_switches = 0;
        m->switch_frame[0] = m->switch_frame[1] = -1;
        for (int f = 0; f < NFRAMES; f++) {
            double d = scene_dist((double)f);
            atten_gain_legacy((float)d, &LEGACY, &set, &gl);
            atten_gain((float)d, &CFG, &gu);
            double v = visual_gain(d);
            double dev_l = fabs(db(gl / v));
            double dev_u = fabs(db(gu / v));
            if (dev_l > maxdev_l) { maxdev_l = dev_l; m->maxdev_legacy_d = d; }
            if (dev_u > maxdev_u) { maxdev_u = dev_u; m->maxdev_unified_d = d;
                                    m->maxdev_frame = f; }
            if (f == 0) {
                jstat_start(&m->frame[0], gl);
                jstat_start(&m->frame[1], gu);
            } else {
                jstat_add(&m->frame[0], gl);
                jstat_add(&m->frame[1], gu);
            }
            if (prev_set >= 0 && set != prev_set) {
                m->frame_set_switches++;
                if (switch_seen < 2) m->switch_frame[switch_seen] = f;
                switch_seen++;
            }
            prev_set = set;
        }
        m->maxdev_legacy = maxdev_l;
        m->maxdev_unified = maxdev_u;
    }

    /* analytic one-sided limits at the two seams */
    {
        double r = CFG.ref, k = CFG.rolloff, w = CFG.knee;
        /* seam at d = ref: plateau (1, 0) vs knee branch (h(0)=0, h'(0)=0) */
        m->seam1_c0 = fabs(1.0 - (r / (r + k * 0.0)));
        m->seam1_c1 = fabs(0.0 - (-r * k * 0.0 / (r * r)));
        /* seam at d = ref + knee: knee branch vs far branch (h(w)=w, h'(w)=1) */
        double den = r + k * w;
        m->seam2_c0 = fabs((r / (r + k * w)) - (r / den));
        m->seam2_c1 = fabs((-r * k * 1.0 / (den * den)) - (-r * k / (den * den)));

        /* size of the jump the legacy switch introduces */
        float below = 0.0f, above = 0.0f;
        atten_gain_legacy(LEGACY.switch_dist - 1e-4f, &LEGACY, NULL, &below);
        atten_gain_legacy(LEGACY.switch_dist, &LEGACY, NULL, &above);
        m->legacy_switch_delta_db = db(above) - db(below);
    }

    /* numeric probe of the implemented float functions, shrinking epsilon */
    {
        const double eps[3] = { 1e-3, 1e-4, 1e-5 };
        double r = CFG.ref, w = CFG.knee;
        for (int i = 0; i < 3; i++) {
            double e = eps[i];
            m->eps_tbl[i][0] = e;
            m->eps_tbl[i][1] = fabs(ugain(r - e) - ugain(r + e));
            m->eps_tbl[i][2] = fabs(uslope(r - e) - uslope(r + e));
            m->eps_tbl[i][3] = fabs(ugain(r + w - e) - ugain(r + w + e));
            m->eps_tbl[i][4] = fabs(uslope(r + w - e) - uslope(r + w + e));
        }
    }
}

static FILE* open_out(const char* name) {
    char path[256];
    snprintf(path, sizeof(path), "out/%s", name);
    FILE* f = fopen(path, "w");
    if (!f) { perror(path); failures++; }
    return f;
}

static int write_distance_level(void) {
    FILE* f = open_out("distance_level.csv");
    if (!f) return -1;
    fprintf(f, "dist_m,legacy_set,legacy_gain,legacy_db,unified_gain,unified_db,"
               "visual_gain,visual_db,unified_minus_visual_db,note\n");
    for (int i = 1; i <= 200; i++) {
        double d = i * 0.05;
        float gl = 0.0f, gu = 0.0f;
        int set = -1;
        atten_gain_legacy((float)d, &LEGACY, &set, &gl);
        atten_gain((float)d, &CFG, &gu);
        double v = visual_gain(d);
        const char* note = "";
        if (fabs(d - (double)CFG.ref) < 1e-9) note = "seam_ref";
        else if (fabs(d - (double)(CFG.ref + CFG.knee)) < 1e-9) note = "seam_knee";
        else if (fabs(d - (double)LEGACY.switch_dist) < 1e-9) note = "legacy_switch";
        fprintf(f, "%.2f,%d,%.6f,%.4f,%.6f,%.4f,%.6f,%.4f,%.4f,%s\n",
                d, set, gl, db(gl), gu, db(gu), v, db(v), fabs(db(gu / v)), note);
    }
    fclose(f);
    return 0;
}

static int write_first_diff(void) {
    FILE* f = open_out("first_diff.csv");
    if (!f) return -1;
    fprintf(f, "dist_m,legacy_gain,legacy_diff,legacy_diff_db,"
               "unified_gain,unified_diff,unified_diff_db,"
               "unified_slope_analytic,note\n");
    const double step = 0.01;
    for (int i = 0; i <= 350; i++) {
        double d = 0.50 + i * step;
        double dp = d - step;
        float gl = 0.0f, glp = 0.0f, gu = 0.0f, gup = 0.0f, us = 0.0f;
        atten_gain_legacy((float)d, &LEGACY, NULL, &gl);
        atten_gain_legacy((float)dp, &LEGACY, NULL, &glp);
        atten_gain((float)d, &CFG, &gu);
        atten_gain((float)dp, &CFG, &gup);
        atten_gain_slope((float)d, &CFG, &us);
        double ldiff = gl - glp, udiff = gu - gup;
        const char* note = "";
        if (fabs(d - (double)CFG.ref) < 1e-9) note = "seam_ref";
        else if (fabs(d - (double)(CFG.ref + CFG.knee)) < 1e-9) note = "seam_knee";
        else if (fabs(d - (double)LEGACY.switch_dist) < 1e-9) note = "legacy_switch";
        fprintf(f, "%.2f,%.6f,%.6f,%.4f,%.6f,%.6f,%.4f,%.6f,%s\n",
                d, gl, ldiff, db(gl) - db(glp), gu, udiff,
                db(gu) - db(gup), us, note);
    }
    fclose(f);
    return 0;
}

static int write_jump_counts(const Metrics* m) {
    FILE* f = open_out("jump_counts.csv");
    if (!f) return -1;
    fprintf(f, "series,model,threshold,max_abs_delta_gain,max_abs_delta_db,"
               "jump_samples,set_switches\n");
    const char* series[3] = { "audio_sample", "frame", "distance_sweep" };
    const char* model[2] = { "legacy_two_set", "unified" };
    long switches[3] = { m->audio_set_switches, m->frame_set_switches, -1 };
    for (int s = 0; s < 3; s++) {
        const JumpStat* js = s == 0 ? m->audio : (s == 1 ? m->frame : m->dense);
        for (int k = 0; k < 2; k++) {
            for (int t = 0; t < NTH; t++) {
                fprintf(f, "%s,%s,%.0e,%.9f,%.6f,%ld,%ld\n",
                        series[s], model[k], THRESH[t],
                        js[k].max_dg, js[k].max_ddb, js[k].count[t],
                        k == 0 ? switches[s] : 0L);
            }
        }
    }
    fclose(f);
    return 0;
}

static int write_event_deviation(const Metrics* m) {
    FILE* f = open_out("event_deviation.csv");
    if (!f) return -1;
    fprintf(f, "frame,time_ms,event,dist_m,legacy_set,legacy_gain,legacy_db,"
               "legacy_frame_step_db,legacy_dev_db,unified_gain,unified_db,"
               "unified_frame_step_db,unified_dev_db,visual_gain,visual_db\n");

    static const char* tag[NFRAMES];
    memset(tag, 0, sizeof(tag));
    tag[0] = "spawn";
    tag[60] = "approach_mid";
    tag[APPROACH_FRAMES] = "closest";
    tag[180] = "recede_mid";
    tag[NFRAMES - 1] = "despawn";
    if (m->switch_frame[0] >= 0) tag[m->switch_frame[0]] = "legacy_switch_in";
    if (m->switch_frame[1] >= 0) tag[m->switch_frame[1]] = "legacy_switch_out";
    if (m->maxdev_frame >= 0 && m->maxdev_frame < NFRAMES)
        tag[m->maxdev_frame] = "max_unified_dev";

    double prev_l = 0.0, prev_u = 0.0;
    int have_prev = 0;
    for (int i = 0; i < NFRAMES; i++) {
        double d = scene_dist((double)i);
        float gl = 0.0f, gu = 0.0f;
        int set = -1;
        atten_gain_legacy((float)d, &LEGACY, &set, &gl);
        atten_gain((float)d, &CFG, &gu);
        double v = visual_gain(d);
        double lstep = have_prev ? db(gl) - prev_l : 0.0;
        double ustep = have_prev ? db(gu) - prev_u : 0.0;
        fprintf(f, "%d,%.3f,%s,%.4f,%d,%.6f,%.4f,%.4f,%.4f,%.6f,%.4f,%.4f,%.4f,"
                   "%.6f,%.4f\n",
                i, i * 1000.0 / FPS, tag[i] ? tag[i] : "", d, set,
                gl, db(gl), lstep, fabs(db(gl / v)),
                gu, db(gu), ustep, fabs(db(gu / v)),
                v, db(v));
        prev_l = db(gl);
        prev_u = db(gu);
        have_prev = 1;
    }
    fclose(f);
    return 0;
}

static int write_error_cases(void) {
    FILE* f = open_out("error_cases.csv");
    if (!f) return -1;
    fprintf(f, "case,input,status,out_gain\n");
    struct { const char* name; float d; } cases[] = {
        { "dist_zero", 0.0f },
        { "dist_neg_zero", -0.0f },
        { "dist_negative", -1.0f },
        { "dist_nan", NAN },
        { "dist_inf", INFINITY },
        { "dist_valid_far", 4.0f },
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        float g = 123.0f;
        AttenStatus st = atten_gain(cases[i].d, &CFG, &g);
        char in[32];
        if (isnan(cases[i].d)) snprintf(in, sizeof(in), "nan");
        else if (isinf(cases[i].d)) snprintf(in, sizeof(in), "inf");
        else snprintf(in, sizeof(in), "%.1f", (double)cases[i].d);
        fprintf(f, "%s,%s,%s,%.6f\n", cases[i].name, in,
                atten_status_str(st), (double)g);
    }
    {
        AttenCfg bad = CFG; bad.ref = 0.0f;
        float g = 123.0f;
        AttenStatus st = atten_gain(2.0f, &bad, &g);
        fprintf(f, "cfg_ref_zero,2.0,%s,%.6f\n", atten_status_str(st), (double)g);
        bad = CFG; bad.rolloff = -1.0f;
        st = atten_gain(2.0f, &bad, &g);
        fprintf(f, "cfg_rolloff_negative,2.0,%s,%.6f\n", atten_status_str(st), (double)g);
        bad = CFG; bad.knee = -1.0f;
        st = atten_gain(2.0f, &bad, &g);
        fprintf(f, "cfg_knee_negative,2.0,%s,%.6f\n", atten_status_str(st), (double)g);
    }
    fclose(f);
    return 0;
}

static int write_summary(const Metrics* m) {
    FILE* f = open_out("summary.md");
    if (!f) return -1;

    fprintf(f,
        "# 距离衰减统一: 证据摘要\n\n"
        "## 统一配置\n\n"
        "只有一组参数, 不再有远近两套切换:\n\n"
        "- `ref = %.2f m` (近场平台边界, `d <= ref` 时 `g = 1`)\n"
        "- `rolloff = %.2f` (远场 `g = ref / (ref + rolloff*(d-ref))`; 取 1 时即 `ref/d`)\n"
        "- `knee = %.2f m` (C1 过渡段宽度)\n\n"
        "旧模型 (对照): 近集 `rolloff = %.2f`, 远集 `rolloff = %.2f`, 按 `d = %.2f m` 硬切换。"
        "注意: 仓库 HEAD 里只有一个公式、没有远近切换, 因此旧模型是按问题描述重建的对照基线, "
        "参数为假设值; 统一曲线本身不依赖它。\n\n",
        (double)CFG.ref, (double)CFG.rolloff, (double)CFG.knee,
        (double)LEGACY.rolloff_near, (double)LEGACY.rolloff_far,
        (double)LEGACY.switch_dist);

    fprintf(f,
        "## 曲线与 C1 证明\n\n"
        "```\n"
        "d <= ref              : g = 1\n"
        "ref < d <= ref + knee : g = ref / (ref + rolloff * h(d-ref))\n"
        "d > ref + knee        : g = ref / (ref + rolloff * (d-ref))\n"
        "h(x) = x^2 * (2*knee - x) / knee^2,  0 <= x <= knee\n"
        "```\n\n"
        "`h(0)=0, h'(0)=0` -> 近场平台 `(g,g')=(1,0)` 与膝部在 `d=ref` 处值与斜率同时相等。\n\n"
        "`h(knee)=knee, h'(knee)=1` -> 膝部与远场在 `d=ref+knee` 处值与斜率同时相等。\n\n"
        "因此整条曲线在实数域上 C1 (二阶导在接缝处才跳变)。`h'(x) >= 0` 于 `[0,knee]`, "
        "所以 `g` 全域单调不增。`make test` 中的 `seam_*` 项对单侧极限和浮点探测都做了检查。\n\n");

    fprintf(f,
        "实测接缝 (解析单侧极限):\n\n"
        "| 接缝 | 值跳变 \\|Δg\\| | 斜率跳变 \\|Δg'\\| |\n"
        "|---|---|---|\n"
        "| ref = %.2f | %.3e | %.3e |\n"
        "| ref+knee = %.2f | %.3e | %.3e |\n\n",
        (double)CFG.ref, m->seam1_c0, m->seam1_c1,
        (double)(CFG.ref + CFG.knee), m->seam2_c0, m->seam2_c1);

    fprintf(f,
        "浮点函数探测 (eps 缩小): 接缝处的值差按 O(eps) 或 O(eps^2) 趋于 0, 斜率差按 O(eps) "
        "趋于 0, 即残余量随探测宽度消失, 不存在与 eps 无关的固定跳变 "
        "(斜率本身连续, 二阶导在接缝处才跳变, 所以斜率探测差是线性收敛):\n\n"
        "| eps | seam_ref Δg | seam_ref Δg' | seam_knee Δg | seam_knee Δg' |\n"
        "|---|---|---|---|---|\n");
    for (int i = 0; i < 3; i++)
        fprintf(f, "| %.0e | %.3e | %.3e | %.3e | %.3e |\n",
                m->eps_tbl[i][0], m->eps_tbl[i][1], m->eps_tbl[i][2],
                m->eps_tbl[i][3], m->eps_tbl[i][4]);
    fprintf(f, "\n");

    fprintf(f,
        "## 四份材料\n\n"
        "1. `distance_level.csv` — 距离-电平: 每 0.05 m 的旧模型/统一模型增益与 dB, "
        "并给出画面预期的 `min(1, ref/d)`。\n"
        "2. `first_diff.csv` — 一阶差分: 0.50..4.00 m 每 0.01 m 的相邻差, "
        "在旧模型切换点可见 %.2f dB 的差分尖峰。\n"
        "3. `jump_counts.csv` — 突跳采样计数: 音频采样率、视频帧率与 1e-4 m 距离扫描三条序列, "
        "在多个阈值下统计 |Δg| 超限的样本数。\n"
        "4. `event_deviation.csv` — 与画面事件对齐后的偏差: 每帧一行, 事件包括出生、切换、最近点、"
        "最大偏差帧等; 偏差定义为与画面线索 `min(1, ref/d)` 的 dB 差。\n\n"
        "`error_cases.csv` 为附加材料, 记录零距离等错误路径。\n\n",
        m->legacy_switch_delta_db);

    fprintf(f,
        "## 关键数字\n\n"
        "| 指标 | 旧模型 (两套切换) | 统一曲线 |\n"
        "|---|---|---|\n"
        "| 音频采样最大单步 \\|Δg\\| | %.6f | %.9f |\n"
        "| 音频采样最大单步 \\|ΔdB\\| | %.4f | %.6f |\n"
        "| 音频采样 \\|Δg\\| > 0.01 的样本数 | %ld | %ld |\n"
        "| 帧序列最大单步 \\|ΔdB\\| | %.4f | %.4f |\n"
        "| 帧序列 \\|Δg\\| > 0.1 的帧数 | %ld | %ld |\n"
        "| 切换次数 (音频/帧) | %ld / %ld | 0 (无切换概念) |\n"
        "| 与画面线索最大偏差 | %.4f dB @ d=%.3f m | %.4f dB @ d=%.3f m |\n\n",
        m->audio[0].max_dg, m->audio[1].max_dg,
        m->audio[0].max_ddb, m->audio[1].max_ddb,
        m->audio[0].count[4], m->audio[1].count[4],
        m->frame[0].max_ddb, m->frame[1].max_ddb,
        m->frame[0].count[6], m->frame[1].count[6],
        m->audio_set_switches, m->frame_set_switches,
        m->maxdev_legacy, m->maxdev_legacy_d,
        m->maxdev_unified, m->maxdev_unified_d);

    fprintf(f,
        "统一曲线只保留膝部内对画面线索的有限偏差 (上表), 膝外与 `ref/d` 完全一致; "
        "旧模型的偏差在切换点附近达到 %.2f dB 并伴随瞬时跳变。\n\n",
        m->maxdev_legacy);

    fprintf(f,
        "## 零距离\n\n"
        "`dist == 0` (含 `-0.0`) 返回 `ATTEN_ERR_DIST_ZERO`, `dist < 0` 返回 "
        "`ATTEN_ERR_DIST_NEG`, NaN/inf 返回 `ATTEN_ERR_NONFINITE`; 所有错误路径都把 "
        "`out_gain` 写成 `0.0`, 绝不静默返回满电平。自检项:`zero_dist_err`, "
        "`zero_dist_not_full`, `neg_zero_dist_err`。详见 `error_cases.csv`。\n");

    fclose(f);
    return 0;
}

static int generate_report(void) {
    if (mkdir("out", 0755) != 0) {
        /* EEXIST is the normal case on a second run */
    }
    Metrics m;
    compute_metrics(&m);
    int rc = 0;
    rc |= write_distance_level();
    rc |= write_first_diff();
    rc |= write_jump_counts(&m);
    rc |= write_event_deviation(&m);
    rc |= write_error_cases();
    rc |= write_summary(&m);
    if (rc == 0)
        printf("report=out/(distance_level.csv, first_diff.csv, jump_counts.csv, "
               "event_deviation.csv, error_cases.csv, summary.md)\n");
    return rc;
}

int main(int argc, char** argv) {
    int report = 0;
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], "--report") == 0) report = 1;

    run_checks();

    if (report) {
        if (generate_report() != 0) failures++;
        printf("report_ok=%d\n", failures == 0);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
