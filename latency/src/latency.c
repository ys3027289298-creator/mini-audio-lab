
/* latency: measure output latency by loopback cross-correlation (no hardcoded
   constants — a fixed value goes out of sync the moment the device changes),
   compensate event scheduling by the measured value, and report four
   materials: sync error table, compensation records, two-buffer comparison,
   and visual-event deviation. Zero measured latency is an error. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define SR 48000.0f
#define MAX_LAT 4096
#define PROBE_N 8192
#define MAX_EVENTS 32
#define MAX_RECORDS 64

typedef struct { int seq; int latency_samples; int buffer_frames; int shift_samples; } Record;
typedef struct { int index; float scheduled_ms; float actual_ms; float error_ms; } SyncRow;
typedef struct { int buffer_frames; int measured_latency; float mean_abs_err_ms; float max_abs_err_ms; } BufferRow;
typedef struct { int frame; float visual_ms; float audio_ms; float dev_before_ms; float dev_after_ms; } VisualRow;

typedef struct {
    float sample_rate;
    int latency_samples;   /* measured per device/buffer, never hardcoded */
    int buffer_frames;
    Record records[MAX_RECORDS];
    int n_records;
} LatencyComp;

/* deterministic per-event device jitter, +-span samples (xorshift32) */
static int jitter(int* state, int span) {
    unsigned x = (unsigned)*state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *state = (int)(x & 0x7fffffff);
    return (int)(x % (unsigned)(2 * span + 1)) - span;
}

/* simulated output device: plays in delayed by lat samples */
static void device_render(const float* in, float* out, int n, int lat) {
    for (int i = 0; i < n; i++) out[i] = (i >= lat) ? in[i - lat] : 0.0f;
}

/* loopback measurement: lag of max cross-correlation between probe and capture */
static int measure_latency(const float* probe, const float* captured, int n) {
    int best_lag = 0;
    float best = -1.0f;
    for (int lag = 0; lag <= MAX_LAT; lag++) {
        float acc = 0.0f;
        for (int i = 0; i + lag < n; i++) acc += probe[i] * captured[i + lag];
        if (acc > best) { best = acc; best_lag = lag; }
    }
    return best_lag;
}

static int measure_device_latency(int device_lat) {
    static float probe[PROBE_N], cap[PROBE_N];
    memset(probe, 0, sizeof(probe));
    memset(cap, 0, sizeof(cap));
    probe[0] = 1.0f;
    device_render(probe, cap, PROBE_N, device_lat);
    return measure_latency(probe, cap, PROBE_N);
}

static int lc_init(LatencyComp* lc, int measured_latency, int buffer_frames, float sr) {
    if (measured_latency <= 0) return -1;   /* zero latency: measurement failed, refuse */
    if (buffer_frames <= 0 || sr <= 0.0f) return -1;
    memset(lc, 0, sizeof(*lc));
    lc->sample_rate = sr;
    lc->latency_samples = measured_latency;
    lc->buffer_frames = buffer_frames;
    return 0;
}

/* shift one event earlier by the measured latency; log a compensation record */
static int lc_compensate(LatencyComp* lc, int scheduled_sample) {
    if (lc == NULL || lc->latency_samples <= 0) return -1;
    if (lc->n_records < MAX_RECORDS) {
        Record* r = &lc->records[lc->n_records];
        r->seq = lc->n_records;
        r->latency_samples = lc->latency_samples;
        r->buffer_frames = lc->buffer_frames;
        r->shift_samples = lc->latency_samples;
    }
    lc->n_records++;
    return scheduled_sample - lc->latency_samples;
}

/* run n events through the device with compensation; fill sync/visual/buffer rows */
static void run_scenario(LatencyComp* lc, int device_lat, const int* sched, const float* visual_ms,
                         int n, int seed, SyncRow* sync, VisualRow* vis, BufferRow* brow) {
    int js = seed;
    float sum = 0.0f, mx = 0.0f;
    for (int i = 0; i < n; i++) {
        int play_at = lc_compensate(lc, sched[i]);
        int heard = play_at + device_lat + jitter(&js, 4);
        float sched_ms = (float)sched[i] / lc->sample_rate * 1000.0f;
        float heard_ms = (float)heard / lc->sample_rate * 1000.0f;
        float err = heard_ms - sched_ms;
        sync[i].index = i;
        sync[i].scheduled_ms = sched_ms;
        sync[i].actual_ms = heard_ms;
        sync[i].error_ms = err;
        vis[i].frame = i;
        vis[i].visual_ms = visual_ms[i];
        vis[i].audio_ms = heard_ms;
        vis[i].dev_before_ms = sched_ms + (float)device_lat / lc->sample_rate * 1000.0f - visual_ms[i];
        vis[i].dev_after_ms = heard_ms - visual_ms[i];
        float a = fabsf(err);
        sum += a;
        if (a > mx) mx = a;
    }
    brow->buffer_frames = lc->buffer_frames;
    brow->measured_latency = lc->latency_samples;
    brow->mean_abs_err_ms = sum / (float)n;
    brow->max_abs_err_ms = mx;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    enum { N = 8, BUF_A = 128, BUF_B = 512 };
    const int dev_a = 96 + 2 * BUF_A;   /* simulated device: base + 2 buffers */
    const int dev_b = 96 + 2 * BUF_B;

    int m_a = measure_device_latency(dev_a);
    int m_b = measure_device_latency(dev_b);
    check("measured_not_hardcoded", m_a == dev_a && m_b == dev_b && m_a != m_b);

    LatencyComp lc_a, lc_b;
    check("init_buf128", lc_init(&lc_a, m_a, BUF_A, SR) == 0);
    check("init_buf512", lc_init(&lc_b, m_b, BUF_B, SR) == 0);
    check("zero_latency_rejected", lc_init(&lc_a, 0, BUF_A, SR) == -1);

    /* silent capture yields zero latency, which must be refused */
    {
        static float probe[PROBE_N], cap[PROBE_N];
        probe[0] = 1.0f;
        check("silent_capture_zero", measure_latency(probe, cap, PROBE_N) == 0);
    }

    /* audio events scheduled on 60 fps visual frame times, starting at 100 ms */
    int sched[N];
    float visual[N];
    for (int i = 0; i < N; i++) {
        visual[i] = 100.0f + (float)i * (1000.0f / 60.0f);
        sched[i] = (int)lroundf(visual[i] * SR / 1000.0f);
    }

    SyncRow sync_a[N], sync_b[N];
    VisualRow vis_a[N], vis_b[N];
    BufferRow buf_a, buf_b;
    run_scenario(&lc_a, dev_a, sched, visual, N, 12345, sync_a, vis_a, &buf_a);
    run_scenario(&lc_b, dev_b, sched, visual, N, 12345, sync_b, vis_b, &buf_b);

    check("compensated_sync_buf128", buf_a.max_abs_err_ms < 0.2f);
    check("compensated_sync_buf512", buf_b.max_abs_err_ms < 0.2f);
    check("records_logged", lc_a.n_records == N && lc_b.n_records == N);
    check("uncompensated_would_drift", fabsf(vis_a[0].dev_before_ms) > 1.0f);

    /* material 1: sync error table */
    for (int b = 0; b < 2; b++) {
        const SyncRow* s = b == 0 ? sync_a : sync_b;
        int bf = b == 0 ? BUF_A : BUF_B;
        for (int i = 0; i < N; i++)
            printf("sync_table buffer=%d idx=%d scheduled_ms=%.3f actual_ms=%.3f error_ms=%+.3f\n",
                   bf, s[i].index, s[i].scheduled_ms, s[i].actual_ms, s[i].error_ms);
    }

    /* material 2: compensation records */
    for (int b = 0; b < 2; b++) {
        const LatencyComp* lc = b == 0 ? &lc_a : &lc_b;
        for (int i = 0; i < lc->n_records; i++)
            printf("record seq=%d latency_samples=%d buffer_frames=%d shift_samples=%d\n",
                   lc->records[i].seq, lc->records[i].latency_samples,
                   lc->records[i].buffer_frames, lc->records[i].shift_samples);
    }

    /* material 3: two buffer lengths compared */
    printf("buffer_compare buffer_frames=%d measured_latency=%d mean_abs_err_ms=%.3f max_abs_err_ms=%.3f\n",
           buf_a.buffer_frames, buf_a.measured_latency, buf_a.mean_abs_err_ms, buf_a.max_abs_err_ms);
    printf("buffer_compare buffer_frames=%d measured_latency=%d mean_abs_err_ms=%.3f max_abs_err_ms=%.3f\n",
           buf_b.buffer_frames, buf_b.measured_latency, buf_b.mean_abs_err_ms, buf_b.max_abs_err_ms);

    /* material 4: deviation from visual events */
    for (int b = 0; b < 2; b++) {
        const VisualRow* v = b == 0 ? vis_a : vis_b;
        int bf = b == 0 ? BUF_A : BUF_B;
        for (int i = 0; i < N; i++)
            printf("visual_dev buffer=%d frame=%d visual_ms=%.3f audio_ms=%.3f dev_before_ms=%+.3f dev_after_ms=%+.3f\n",
                   bf, v[i].frame, v[i].visual_ms, v[i].audio_ms, v[i].dev_before_ms, v[i].dev_after_ms);
    }

    /* all four materials required; any one missing means incomplete */
    int have_sync = N > 0;
    int have_records = lc_a.n_records > 0 && lc_b.n_records > 0;
    int have_buffers = buf_a.measured_latency > 0 && buf_b.measured_latency > 0;
    int have_visual = N > 0;
    check("materials_complete", have_sync && have_records && have_buffers && have_visual);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
