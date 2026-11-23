/* gain.c: dB-domain bus gain structure. See gain.h for the model. */
#include "gain.h"

#include <float.h>
#include <math.h>
#include <string.h>

#define GAIN_TARGET_MIN_DB (-120.0f)

/* ---- dB helpers ------------------------------------------------------- */

float gain_db_from_amp(float amp) {
    if (!(amp > 0.0f)) return -INFINITY;   /* 0, negative, NaN */
    return 20.0f * log10f(amp);
}

float gain_amp_from_db(float db) {
    if (isinf(db) && db < 0.0f) return 0.0f;
    return powf(10.0f, db / 20.0f);
}

/* Coherent sum: amplitudes add. Silence contributes 0; any NaN poisons. */
float gain_sum_coherent_db(const float* level_db, int n) {
    if (level_db == NULL || n <= 0) return -INFINITY;
    double a = 0.0;
    for (int i = 0; i < n; i++) {
        float db = level_db[i];
        if (isnan(db)) return NAN;
        if (isinf(db)) {
            if (db > 0.0f) return INFINITY;
            continue;                      /* -inf: silent */
        }
        a += pow(10.0, (double)db / 20.0);
    }
    return a > 0.0 ? (float)(20.0 * log10(a)) : -INFINITY;
}

/* Uncorrelated sum: powers add. */
float gain_sum_power_db(const float* level_db, int n) {
    if (level_db == NULL || n <= 0) return -INFINITY;
    double p = 0.0;
    for (int i = 0; i < n; i++) {
        float db = level_db[i];
        if (isnan(db)) return NAN;
        if (isinf(db)) {
            if (db > 0.0f) return INFINITY;
            continue;
        }
        p += pow(10.0, (double)db / 10.0);
    }
    return p > 0.0 ? (float)(10.0 * log10(p)) : -INFINITY;
}

/* ---- plan ------------------------------------------------------------- */

int gain_plan_init(GainPlan* plan, const float* peak_dbfs, int n, float target_dbfs) {
    if (plan == NULL || peak_dbfs == NULL) return -1;
    if (n <= 0 || n > GAIN_MAX_CH) return -1;
    if (isnan(target_dbfs) || isinf(target_dbfs)) return -1;
    if (target_dbfs > 0.0f || target_dbfs < GAIN_TARGET_MIN_DB) return -1;

    for (int i = 0; i < n; i++) {
        float db = peak_dbfs[i];
        if (isnan(db)) return -1;
        if (isinf(db) && db > 0.0f) return -1;              /* +inf level */
        if (!isinf(db) && (db > 0.0f || db < GAIN_TARGET_MIN_DB)) return -1;
    }

    memset(plan, 0, sizeof(*plan));
    plan->channels = n;
    plan->target_dbfs = target_dbfs;
    plan->coherent_sum_dbfs = gain_sum_coherent_db(peak_dbfs, n);
    plan->power_sum_dbfs = gain_sum_power_db(peak_dbfs, n);

    if (isinf(plan->coherent_sum_dbfs) && plan->coherent_sum_dbfs < 0.0f) {
        /* every channel silent: nothing to trim, nothing to allocate */
        plan->bus_headroom_db = INFINITY;
        plan->rms_margin_db = INFINITY;
        plan->post_sum_dbfs = -INFINITY;
        for (int i = 0; i < n; i++) {
            plan->silent[i] = 1;
            plan->trim_db[i] = 0.0f;
            plan->share[i] = 0.0f;
            plan->alloc_dbfs[i] = -INFINITY;
            plan->margin_db[i] = 0.0f;
        }
        return 0;
    }

    float head = target_dbfs - plan->coherent_sum_dbfs;
    plan->bus_headroom_db = head;
    plan->rms_margin_db = target_dbfs - plan->power_sum_dbfs;

    /* Shared trim: move the coherent sum to the target. Under target this is
     * make-up gain (so a lone quiet source becomes audible), capped so a
     * near-silent input is not amplified without limit. */
    float trim = head;
    if (trim > GAIN_MAX_MAKEUP_DB) trim = GAIN_MAX_MAKEUP_DB;
    double s = 0.0;
    for (int i = 0; i < n; i++)
        if (!isinf(peak_dbfs[i])) s += pow(10.0, (double)peak_dbfs[i] / 20.0);

    for (int i = 0; i < n; i++) {
        float db = peak_dbfs[i];
        if (isinf(db) && db < 0.0f) {
            plan->silent[i] = 1;
            plan->trim_db[i] = 0.0f;
            plan->share[i] = 0.0f;
            plan->alloc_dbfs[i] = -INFINITY;
            plan->margin_db[i] = 0.0f;
            continue;
        }
        double amp = pow(10.0, (double)db / 20.0);
        plan->trim_db[i] = trim;
        plan->share[i] = s > 0.0 ? (float)(amp / s) : 0.0f;
        plan->alloc_dbfs[i] = db + trim;
        plan->margin_db[i] = -plan->alloc_dbfs[i];   /* dB below 0 dBFS */
    }
    plan->post_sum_dbfs = plan->coherent_sum_dbfs + trim;
    return 0;
}

/* ---- meters ----------------------------------------------------------- */

void gain_meter_reset(GainMeter* m) {
    if (m == NULL) return;
    m->samples = 0;
    m->peak = 0.0f;
    m->sum_sq = 0.0;
    m->clipped = 0;
}

void gain_meter_feed(GainMeter* m, const float* x, int n) {
    if (m == NULL || x == NULL || n <= 0) return;
    for (int i = 0; i < n; i++) {
        float a = fabsf(x[i]);
        if (a > m->peak) m->peak = a;
        if (a > 1.0f) m->clipped++;
        m->sum_sq += (double)x[i] * (double)x[i];
        m->samples++;
    }
}

float gain_meter_peak_db(const GainMeter* m) {
    if (m == NULL || m->samples <= 0 || !(m->peak > 0.0f)) return -INFINITY;
    return gain_db_from_amp(m->peak);
}

float gain_meter_rms_db(const GainMeter* m) {
    if (m == NULL || m->samples <= 0 || !(m->sum_sq > 0.0)) return -INFINITY;
    return gain_db_from_amp((float)sqrt(m->sum_sq / (double)m->samples));
}

/* ---- bus summing ------------------------------------------------------ */

static void span_open(ClipSpan* s, int i, float v) {
    s->start = i;
    s->peak = fabsf(v);
    s->lo = v;
    s->hi = v;
    s->pos = v > 0.0f ? 1 : 0;
    s->neg = v < 0.0f ? 1 : 0;
}

static void span_extend(ClipSpan* s, float v) {
    float a = fabsf(v);
    if (a > s->peak) s->peak = a;
    if (v < s->lo) s->lo = v;
    if (v > s->hi) s->hi = v;
    if (v > 0.0f) s->pos++;
    if (v < 0.0f) s->neg++;
}

static void span_close(ClipList* list, ClipSpan* s, int end,
                       ClipSpan* spans, int span_cap) {
    s->end = end;
    list->count++;
    if (spans != NULL && list->written < span_cap) spans[list->written++] = *s;
    else if (span_cap > 0) list->overflow = 1;
}

int gain_bus_sum(const float* const* ch, int nch, int frames,
                 const GainPlan* plan, int ramp_frames, float* out,
                 GainMeter* meter, ClipSpan* spans, int span_cap, ClipList* list) {
    if (ch == NULL || out == NULL) return -1;
    if (nch <= 0 || nch > GAIN_MAX_CH || frames <= 0) return -1;
    if (ramp_frames < 0 || ramp_frames > frames) return -1;
    if (plan != NULL && plan->channels != nch) return -1;
    if (span_cap < 0) return -1;
    if (span_cap > 0 && spans == NULL) return -1;
    if (spans != NULL && list == NULL) return -1;
    for (int c = 0; c < nch; c++) if (ch[c] == NULL) return -1;

    if (list != NULL) { list->count = 0; list->written = 0; list->overflow = 0; }

    float gain[GAIN_MAX_CH];
    for (int c = 0; c < nch; c++)
        gain[c] = plan != NULL ? gain_amp_from_db(plan->trim_db[c]) : 1.0f;

    GainMeter local;
    GainMeter* m = meter != NULL ? meter : &local;
    if (meter == NULL) gain_meter_reset(m);

    ClipSpan cur;
    int have_span = 0;

    for (int i = 0; i < frames; i++) {
        float t = ramp_frames > 0 ? (float)i / (float)ramp_frames : 1.0f;
        if (t > 1.0f) t = 1.0f;
        float v = 0.0f;
        for (int c = 0; c < nch; c++) v += ch[c][i] * gain[c] * t;

        if (list != NULL || spans != NULL) {
            int over = v > 1.0f || v < -1.0f;
            if (over) {
                if (!have_span) { span_open(&cur, i, v); have_span = 1; }
                else span_extend(&cur, v);
            } else if (have_span) {
                span_close(list, &cur, i, spans, span_cap);
                have_span = 0;
            }
        }
        if (v > 1.0f || v < -1.0f) m->clipped++;
        if (v > 1.0f) v = 1.0f;
        else if (v < -1.0f) v = -1.0f;
        out[i] = v;
        m->peak = fabsf(v) > m->peak ? fabsf(v) : m->peak;
        m->sum_sq += (double)v * (double)v;
        m->samples++;
    }
    if (have_span) span_close(list, &cur, frames, spans, span_cap);

    return frames;
}

int mix_bus_auto(const float* const* ch, int nch, int frames, float target_dbfs,
                 float* out, GainPlan* used_plan, GainMeter* meter,
                 ClipSpan* spans, int span_cap, ClipList* list) {
#if MIX_GAIN_GROUP
    if (ch == NULL || nch <= 0 || nch > GAIN_MAX_CH || frames <= 0) return -1;
    float level[GAIN_MAX_CH];
    for (int c = 0; c < nch; c++) {
        if (ch[c] == NULL) return -1;
        float peak = 0.0f;
        for (int i = 0; i < frames; i++) {
            float a = fabsf(ch[c][i]);
            if (a > peak) peak = a;
        }
        level[c] = gain_db_from_amp(peak);
    }
    GainPlan plan;
    if (gain_plan_init(&plan, level, nch, target_dbfs) != 0) return -1;
    if (used_plan != NULL) *used_plan = plan;
    return gain_bus_sum(ch, nch, frames, &plan, 0, out, meter, spans, span_cap, list);
#else
    (void)target_dbfs; (void)used_plan;
    return gain_bus_sum(ch, nch, frames, NULL, 0, out, meter, spans, span_cap, list);
#endif
}
