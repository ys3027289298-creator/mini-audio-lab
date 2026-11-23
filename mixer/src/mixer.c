
/* mixer: float mixing with ramp gain and clipping detection.
 *
 * Clipping is reported as maximal contiguous sample intervals (spans) with a
 * half-open [start, end) range, plus a waveform summary. Build with
 * -DMIX_CLIP_LOCATE=0 to restore the original occurrence-only detector; that
 * build's output is unchanged.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define N 128

#ifndef MIX_CLIP_LOCATE
#define MIX_CLIP_LOCATE 1
#endif

/* One maximal run of clipped samples and its pre-clamp extent. */
typedef struct {
    int start, end;      /* half-open [start, end) */
    float peak;          /* max |v| before clamping */
    float lo, hi;        /* min/max v before clamping */
    int pos, neg;        /* clipped samples by polarity */
} ClipSpan;

typedef struct {
    float peak;
    int clipped;
    int samples;
#if MIX_CLIP_LOCATE
    int span_count;      /* total clipped intervals found */
    int span_written;    /* intervals copied to the caller's buffer */
    int span_overflow;   /* caller's buffer was too small */
    int first_clip;      /* sample index, -1 when none clipped */
    int last_clip;
    int clipped_pos;
    int clipped_neg;
    float rms;
    float dc;
    float crest;
#endif
} MixStat;

#if MIX_CLIP_LOCATE
/* Mix and locate every clipped sample interval. spans may be NULL; span_cap
 * bounds how many intervals are copied out while span_count still reports the
 * true total. */
static int mix_ramp_locate(const float* a, const float* b, float gain_a, float gain_b,
                           float* out, int n, int ramp_frames, MixStat* st,
                           ClipSpan* spans, int span_cap) {
    if (n <= 0) return -1;
    if (ramp_frames < 0 || ramp_frames > n) return -1;
    if (gain_a < 0.0f || gain_b < 0.0f) return -1;
    if (span_cap < 0 || (span_cap > 0 && spans == NULL)) return -1;

    float peak = 0.0f;
    double sq = 0.0, sum = 0.0;
    int clipped = 0, cpos = 0, cneg = 0;
    int first = -1, last = -1;
    int span_count = 0, span_written = 0, span_overflow = 0;
    int in_span = 0;
    ClipSpan cur;

    for (int i = 0; i < n; i++) {
        float t = ramp_frames > 0 ? (float)i / (float)ramp_frames : 1.0f;
        if (t > 1.0f) t = 1.0f;
        float ga = gain_a * t, gb = gain_b * t;
        float v = a[i] * ga + b[i] * gb;

        if (v > 1.0f || v < -1.0f) {
            clipped++;
            if (v > 0.0f) cpos++; else cneg++;
            if (first < 0) first = i;
            last = i;
            if (!in_span) {
                in_span = 1;
                cur.start = i;
                cur.peak = fabsf(v);
                cur.lo = v;
                cur.hi = v;
                cur.pos = 0;
                cur.neg = 0;
            } else {
                float m = fabsf(v);
                if (m > cur.peak) cur.peak = m;
                if (v < cur.lo) cur.lo = v;
                if (v > cur.hi) cur.hi = v;
            }
            if (v > 0.0f) cur.pos++; else cur.neg++;
        } else if (in_span) {
            cur.end = i;
            span_count++;
            if (spans != NULL && span_written < span_cap) {
                spans[span_written++] = cur;
            } else if (span_cap > 0) {
                span_overflow = 1;
            }
            in_span = 0;
        }

        if (v > 1.0f) v = 1.0f;
        else if (v < -1.0f) v = -1.0f;
        out[i] = v;
        if (fabsf(v) > peak) peak = fabsf(v);
        double d = (double)v;
        sq += d * d;
        sum += d;
    }
    if (in_span) {
        cur.end = n;
        span_count++;
        if (spans != NULL && span_written < span_cap) {
            spans[span_written++] = cur;
        } else if (span_cap > 0) {
            span_overflow = 1;
        }
    }

    if (st) {
        st->peak = peak;
        st->clipped = clipped;
        st->samples = n;
        st->span_count = span_count;
        st->span_written = span_written;
        st->span_overflow = span_overflow;
        st->first_clip = first;
        st->last_clip = last;
        st->clipped_pos = cpos;
        st->clipped_neg = cneg;
        st->rms = (float)sqrt(sq / (double)n);
        st->dc = (float)(sum / (double)n);
        st->crest = st->rms > 0.0f ? peak / st->rms : 0.0f;
    }
    return n;
}
#endif

static int mix_ramp(const float* a, const float* b, float gain_a, float gain_b,
                    float* out, int n, int ramp_frames, MixStat* st) {
#if MIX_CLIP_LOCATE
    return mix_ramp_locate(a, b, gain_a, gain_b, out, n, ramp_frames, st, NULL, 0);
#else
    if (n <= 0) return -1;
    if (ramp_frames < 0 || ramp_frames > n) return -1;
    if (gain_a < 0.0f || gain_b < 0.0f) return -1;
    float peak = 0.0f; int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = ramp_frames > 0 ? (float)i / (float)ramp_frames : 1.0f;
        if (t > 1.0f) t = 1.0f;
        float ga = gain_a * t, gb = gain_b * t;
        float v = a[i] * ga + b[i] * gb;
        if (v > 1.0f || v < -1.0f) { clipped++; if (v > 1.0f) v = 1.0f; else v = -1.0f; }
        out[i] = v;
        if (fabsf(v) > peak) peak = fabsf(v);
    }
    if (st) { st->peak = peak; st->clipped = clipped; st->samples = n; }
    return n;
#endif
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#if MIX_CLIP_LOCATE
/* ---- fault-injection evidence harness ---------------------------------- */

#define EV_N 512
#define EV_MAXF 8
#define EV_MAXS 8

typedef struct { int start[EV_MAXF], end[EV_MAXF], n; } Fact;

typedef struct {
    const char* name;
    Fact fact;
    MixStat st;
    ClipSpan span[EV_MAXS];
    int kept;
    int exact;
    float recall, precision;
} Case;

enum { C_CLEAN, C_EARLY, C_LATE, C_MIXED, C_ADJ, C_THRESH, C_GAP, C_RAMP, C_SPEC, C_OVF, C_NCASE };

static void ev_base(float* a, float* b, int n) {
    for (int i = 0; i < n; i++) {
        a[i] = 0.30f * sinf(2.0f * 3.14159265f * (float)i / 64.0f);
        b[i] = 0.20f * cosf(2.0f * 3.14159265f * (float)i / 96.0f);
    }
}

static void ev_inject(float* a, float* b, int s, int e, float va, float vb) {
    for (int i = s; i < e; i++) { a[i] = va; b[i] = vb; }
}

static void ev_fact(Fact* f, int s, int e) {
    f->start[f->n] = s;
    f->end[f->n] = e;
    f->n++;
}

/* Plant known clipping defects; fact spans are the expected detector output. */
static void ev_prepare(int id, float* a, float* b, int* ramp, int* cap, Fact* f) {
    ev_base(a, b, EV_N);
    memset(f, 0, sizeof(*f));
    *ramp = 0;
    *cap = EV_MAXS;
    switch (id) {
    case C_CLEAN:
        break;
    case C_EARLY:                                /* 0.95 + 0.35 = 1.30 */
        ev_inject(a, b, 40, 56, 0.95f, 0.35f);
        ev_fact(f, 40, 56);
        break;
    case C_LATE:                                 /* same count, other place */
        ev_inject(a, b, 400, 416, 0.95f, 0.35f);
        ev_fact(f, 400, 416);
        break;
    case C_MIXED:
        ev_inject(a, b, 40, 56, 0.95f, 0.35f);
        ev_fact(f, 40, 56);
        ev_inject(a, b, 100, 103, -0.95f, -0.35f);
        ev_fact(f, 100, 103);
        ev_inject(a, b, 300, 302, 0.80f, 0.40f);
        ev_fact(f, 300, 302);
        ev_inject(a, b, 400, 460, -0.90f, -0.30f);
        ev_fact(f, 400, 460);
        break;
    case C_ADJ:                                  /* touching runs merge */
        ev_inject(a, b, 200, 204, 0.90f, 0.40f);
        ev_inject(a, b, 204, 208, 0.90f, 0.40f);
        ev_fact(f, 200, 208);
        break;
    case C_THRESH:
        ev_inject(a, b, 60, 64, 0.60f, 0.40f);
        ev_inject(a, b, 80, 84, 0.70f, 0.30f);
        ev_inject(a, b, 100, 104, 0.70f, 0.31f);
        ev_fact(f, 100, 104);
        break;
    case C_GAP:                                  /* one clean sample splits */
        ev_inject(a, b, 220, 224, 0.90f, 0.40f);
        ev_fact(f, 220, 224);
        ev_inject(a, b, 225, 229, 0.90f, 0.40f);
        ev_fact(f, 225, 229);
        break;
    case C_RAMP:                                 /* full ramp: clip when 1.8*t > 1 */
        for (int i = 0; i < EV_N; i++) { a[i] = 0.9f; b[i] = 0.9f; }
        *ramp = EV_N;
        ev_fact(f, (int)ceilf((float)EV_N / 1.8f), EV_N);
        break;
    case C_SPEC:                                 /* hot decoys at 0.99 must stay clean */
        ev_inject(a, b, 60, 80, 0.60f, 0.39f);
        ev_inject(a, b, 150, 170, 0.60f, 0.39f);
        ev_inject(a, b, 450, 470, 0.60f, 0.39f);
        ev_inject(a, b, 350, 353, 0.95f, 0.35f);
        ev_fact(f, 350, 353);
        break;
    case C_OVF:                                  /* 4 spans, buffer holds 2 */
        ev_inject(a, b, 40, 56, 0.95f, 0.35f);
        ev_fact(f, 40, 56);
        ev_inject(a, b, 100, 103, -0.95f, -0.35f);
        ev_fact(f, 100, 103);
        ev_inject(a, b, 300, 302, 0.80f, 0.40f);
        ev_fact(f, 300, 302);
        ev_inject(a, b, 400, 460, -0.90f, -0.30f);
        ev_fact(f, 400, 460);
        *cap = 2;
        break;
    }
}

static void mask_from_fact(const Fact* f, unsigned char* m, int n) {
    memset(m, 0, (size_t)n);
    for (int k = 0; k < f->n; k++)
        for (int i = f->start[k]; i < f->end[k]; i++) m[i] = 1;
}

static void mask_from_spans(const Case* c, unsigned char* m, int n) {
    memset(m, 0, (size_t)n);
    for (int k = 0; k < c->kept; k++)
        for (int i = c->span[k].start; i < c->span[k].end; i++) m[i] = 1;
}

static void eval_capture(Case* c) {
    unsigned char fm[EV_N], dm[EV_N];
    mask_from_fact(&c->fact, fm, EV_N);
    mask_from_spans(c, dm, EV_N);
    long fs = 0, ds = 0, hit = 0;
    for (int i = 0; i < EV_N; i++) { fs += fm[i]; ds += dm[i]; hit += (fm[i] && dm[i]); }
    c->recall = fs > 0 ? (float)hit / (float)fs : 0.0f;
    c->precision = ds > 0 ? (float)hit / (float)ds : 0.0f;
    c->exact = (c->st.span_count == c->fact.n);
    for (int k = 0; k < c->fact.n && c->exact; k++)
        if (c->span[k].start != c->fact.start[k] || c->span[k].end != c->fact.end[k])
            c->exact = 0;
}

static int spans_equal(const ClipSpan* a, int na, const ClipSpan* b, int nb) {
    if (na != nb) return 0;
    for (int k = 0; k < na; k++)
        if (a[k].start != b[k].start || a[k].end != b[k].end) return 0;
    return 1;
}

static int run_evidence(void) {
    static const char* names[C_NCASE] = {
        "clean", "pos_early", "pos_late", "mixed", "adjacent", "threshold",
        "gap", "ramp", "specificity", "overflow"
    };
    Case cs[C_NCASE];
    float a[EV_N], b[EV_N], out[EV_N];

    for (int id = 0; id < C_NCASE; id++) {
        int ramp, cap;
        memset(&cs[id], 0, sizeof(cs[id]));
        ev_prepare(id, a, b, &ramp, &cap, &cs[id].fact);
        cs[id].name = names[id];
        mix_ramp_locate(a, b, 1.0f, 1.0f, out, EV_N, ramp, &cs[id].st, cs[id].span, cap);
        cs[id].kept = cs[id].st.span_written;
    }

    /* --- material 1: clipped interval list --- */
    printf("[clip_spans]\n");
    for (int id = 0; id < C_NCASE; id++) {
        const Case* c = &cs[id];
        if (id == C_OVF) {
            printf("case=%s spans=%d kept=%d overflow=%d\n",
                   c->name, c->st.span_count, c->kept, c->st.span_overflow);
            for (int k = 0; k < c->kept; k++)
                printf("case=%s span=%d start=%d end=%d len=%d\n", c->name, k,
                       c->span[k].start, c->span[k].end, c->span[k].end - c->span[k].start);
            continue;
        }
        printf("case=%s spans=%d\n", c->name, c->st.span_count);
        for (int k = 0; k < c->kept; k++) {
            const ClipSpan* s = &c->span[k];
            printf("case=%s span=%d start=%d end=%d len=%d peak=%.4f lo=%.4f hi=%.4f pos=%d neg=%d\n",
                   c->name, k, s->start, s->end, s->end - s->start,
                   s->peak, s->lo, s->hi, s->pos, s->neg);
        }
    }

    /* --- material 2: waveform summary --- */
    printf("[wave_summary]\n");
    for (int id = 0; id < C_NCASE; id++) {
        const MixStat* s = &cs[id].st;
        printf("case=%s samples=%d peak=%.4f rms=%.4f dc=%.4f crest=%.4f "
               "clipped=%d pos=%d neg=%d first=%d last=%d spans=%d\n",
               cs[id].name, s->samples, s->peak, s->rms, s->dc, s->crest,
               s->clipped, s->clipped_pos, s->clipped_neg, s->first_clip, s->last_clip,
               s->span_count);
    }

    /* --- material 3: capture rate after injection --- */
    long fact_spans = 0, hit_spans = 0, fact_samp = 0, det_samp = 0, hit_samp = 0;
    static const int cap_ids[] = { C_EARLY, C_LATE, C_MIXED, C_ADJ, C_THRESH, C_GAP, C_RAMP, C_SPEC };
    for (size_t k = 0; k < sizeof(cap_ids) / sizeof(cap_ids[0]); k++) {
        Case* c = &cs[cap_ids[k]];
        eval_capture(c);
        fact_spans += c->fact.n;
        hit_spans += c->exact ? c->fact.n : 0;
        unsigned char fm[EV_N], dm[EV_N];
        mask_from_fact(&c->fact, fm, EV_N);
        mask_from_spans(c, dm, EV_N);
        for (int i = 0; i < EV_N; i++) {
            fact_samp += fm[i];
            det_samp += dm[i];
            hit_samp += (fm[i] && dm[i]);
        }
        printf("[injection_capture]\n");
        printf("case=%s injected=%d detected=%d exact=%d sample_recall=%.4f sample_precision=%.4f\n",
               c->name, c->fact.n, c->st.span_count, c->exact, c->recall, c->precision);
    }
    printf("capture_rate=%.4f fact_spans=%ld located_exact=%ld\n",
           fact_spans ? (double)hit_spans / (double)fact_spans : 0.0, fact_spans, hit_spans);
    printf("sample_recall=%.4f sample_precision=%.4f\n",
           fact_samp ? (double)hit_samp / (double)fact_samp : 0.0,
           det_samp ? (double)hit_samp / (double)det_samp : 0.0);

    /* Control: the occurrence-only entry point reports the same clipped counts
     * but writes no intervals, so its located-span capture rate is zero. */
    {
        long ctrl_fact = 0, ctrl_located = 0;
        int parity = 1;
        for (size_t k = 0; k < sizeof(cap_ids) / sizeof(cap_ids[0]); k++) {
            Case* c = &cs[cap_ids[k]];
            int ramp, cap; Fact f;
            MixStat legacy;
            memset(&legacy, 0, sizeof(legacy));
            ev_prepare(cap_ids[k], a, b, &ramp, &cap, &f);
            mix_ramp(a, b, 1.0f, 1.0f, out, EV_N, ramp, &legacy);
            ctrl_fact += c->fact.n;
            ctrl_located += legacy.span_written;   /* stays 0 by construction */
            if (legacy.clipped != c->st.clipped || legacy.peak != c->st.peak) parity = 0;
        }
        printf("control_no_location_rate=%.4f fact_spans=%ld count_parity=%d\n",
               ctrl_fact ? (double)ctrl_located / (double)ctrl_fact : 0.0,
               ctrl_fact, parity);
    }
    /* Control: the validator must reject a located list shifted by one sample,
     * otherwise exact=1 above would not be evidence of correct placement. */
    {
        Case tmp = cs[C_MIXED];
        for (int k = 0; k < tmp.kept; k++) { tmp.span[k].start++; tmp.span[k].end++; }
        printf("control_shift_rejected=%d\n",
               spans_equal(tmp.span, tmp.kept, cs[C_MIXED].span, cs[C_MIXED].kept) == 0);
    }

    /* --- material 4: coverage checklist --- */
    printf("[coverage]\n");
    check("cov_clean_no_spans", cs[C_CLEAN].st.span_count == 0 && cs[C_CLEAN].st.clipped == 0);
    check("cov_single_span_located",
          cs[C_EARLY].exact && cs[C_EARLY].span[0].start == 40 && cs[C_EARLY].span[0].end == 56);
    check("cov_multi_span_located", cs[C_MIXED].exact && cs[C_MIXED].st.span_count == 4);
    check("cov_adjacent_merges", cs[C_ADJ].exact &&
          cs[C_ADJ].span[0].start == 200 && cs[C_ADJ].span[0].end == 208);
    check("cov_gap_splits", cs[C_GAP].exact && cs[C_GAP].st.span_count == 2);
    check("cov_ramp_onset_located", cs[C_RAMP].exact &&
          cs[C_RAMP].span[0].start == (int)ceilf((float)EV_N / 1.8f));
    check("cov_span_cap_truncates",
          cs[C_OVF].st.span_count == 4 && cs[C_OVF].kept == 2 && cs[C_OVF].st.span_overflow == 1);
    check("cov_threshold_boundary",
          cs[C_THRESH].st.span_count == 1 && cs[C_THRESH].st.clipped == 4 &&
          cs[C_THRESH].span[0].start == 100 && cs[C_THRESH].span[0].end == 104);
    check("cov_hot_decoy_clean",
          cs[C_SPEC].st.span_count == 1 && cs[C_SPEC].st.clipped == 3 &&
          cs[C_SPEC].span[0].start == 350 && cs[C_SPEC].span[0].end == 353 &&
          cs[C_SPEC].st.peak <= 1.0f + 1e-6f);
    check("cov_output_clamped", cs[C_MIXED].st.peak <= 1.0f + 1e-6f);
    check("cov_polarity_counts",
          cs[C_MIXED].st.clipped_pos == 18 && cs[C_MIXED].st.clipped_neg == 63);
    check("cov_first_last_clip",
          cs[C_MIXED].st.first_clip == 40 && cs[C_MIXED].st.last_clip == 459);
    check("cov_wave_summary", cs[C_MIXED].st.samples == EV_N &&
          cs[C_MIXED].st.rms > 0.0f && cs[C_MIXED].st.crest >= 1.0f - 1e-5f);
    /* occurrence-only report cannot tell these two apart; spans can. */
    check("cov_equal_count_needs_location",
          cs[C_EARLY].st.clipped == cs[C_LATE].st.clipped &&
          cs[C_EARLY].span[0].start != cs[C_LATE].span[0].start);
    {
        /* same numbers through the legacy entry point */
        int ramp, cap; Fact f;
        MixStat via_legacy;
        ev_prepare(C_MIXED, a, b, &ramp, &cap, &f);
        memset(&via_legacy, 0, sizeof(via_legacy));
        int r = mix_ramp(a, b, 1.0f, 1.0f, out, EV_N, ramp, &via_legacy);
        check("cov_legacy_entry_parity",
              r == EV_N && via_legacy.clipped == cs[C_MIXED].st.clipped &&
              via_legacy.peak == cs[C_MIXED].st.peak &&
              via_legacy.span_count == 4 && via_legacy.span_written == 0 &&
              via_legacy.span_overflow == 0);
    }
    {
        ClipSpan sp[EV_MAXS];
        int ramp, cap; Fact f;
        ev_prepare(C_MIXED, a, b, &ramp, &cap, &f);
        check("cov_invalid_span_args",
              mix_ramp_locate(a, b, 1.0f, 1.0f, out, EV_N, ramp, NULL, sp, -1) == -1 &&
              mix_ramp_locate(a, b, 1.0f, 1.0f, out, EV_N, ramp, NULL, NULL, 4) == -1);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif /* MIX_CLIP_LOCATE */

int main(int argc, char** argv) {
#if MIX_CLIP_LOCATE
    if (argc > 1 && strcmp(argv[1], "--evidence") == 0) return run_evidence();
#else
    (void)argc; (void)argv;
#endif
    float a[N], b[N], out[N];
    for (int i = 0; i < N; i++) { a[i] = 0.5f; b[i] = 0.25f; }

    MixStat st;
    check("mix_ok", mix_ramp(a, b, 1.0f, 1.0f, out, N, 0, &st) == N);
    check("mix_sum", fabsf(out[N - 1] - 0.75f) < 1e-5f);
    check("no_clip", st.clipped == 0);

    check("ramp_bad", mix_ramp(a, b, 1.0f, 1.0f, out, N, N + 1, &st) == -1);
    check("neg_gain", mix_ramp(a, b, -1.0f, 1.0f, out, N, 0, &st) == -1);

    float loud[N];
    for (int i = 0; i < N; i++) loud[i] = 0.9f;
    mix_ramp(loud, loud, 1.0f, 1.0f, out, N, 0, &st);
    check("clipped_detected", st.clipped == N && st.peak <= 1.0f + 1e-6f);

    /* ramp must start near zero */
    mix_ramp(a, b, 1.0f, 1.0f, out, N, N, &st);
    check("ramp_start_silent", fabsf(out[0]) < 1e-6f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
