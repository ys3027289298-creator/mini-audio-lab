/* gain.h: channel-count-aware bus gain structure, level metering and clip
 * interval location.
 *
 * The bus is summed in the decibel domain:
 *   coherent (worst case, all channels peak together):
 *       S_peak = 20 * log10( sum_i 10^(L_i / 20) )
 *   uncorrelated (power sum, program material):
 *       S_rms  = 10 * log10( sum_i 10^(L_i / 10) )
 * where L_i is the measured peak level of channel i in dBFS.
 *
 * Bus headroom is derived from the coherent sum, so a plan can guarantee that
 * no sample of the summed bus exceeds target_dbfs even when every channel
 * peaks on the same frame. The same trim is applied to every channel, which
 * keeps the mix balance intact. A single source therefore gets the whole bus
 * budget (no more "one channel is too quiet"), while N sources share it
 * (no more "stacking N channels clips").
 */
#ifndef MIX_GAIN_H
#define MIX_GAIN_H

#include <math.h>

/* Build-time switch: 0 restores the pre-change unity-gain summing path.
 * The binary built with MIX_GAIN_GROUP=0 must be byte-identical to the
 * original mixer output. */
#ifndef MIX_GAIN_GROUP
#define MIX_GAIN_GROUP 1
#endif

#define GAIN_MAX_CH 64

/* Upper bound on make-up gain used to lift a quiet bus to its target. Stops a
 * near-silent input from being amplified without limit. */
#define GAIN_MAX_MAKEUP_DB 24.0f

/* ---- level meter ------------------------------------------------------ */

typedef struct {
    long samples;
    float peak;       /* linear |x| peak */
    double sum_sq;    /* double accumulator */
    int clipped;      /* |x| > 1.0 samples */
} GainMeter;

/* ---- clip interval location ------------------------------------------- */

/* One maximal run of samples with |v| > 1.0 before clamping. */
typedef struct {
    int start, end;    /* half-open [start, end) */
    float peak;        /* max |v| pre-clamp */
    float lo, hi;      /* min/max v pre-clamp */
    long pos, neg;     /* clipped samples by polarity */
} ClipSpan;

typedef struct {
    int count;         /* true number of intervals */
    int written;       /* intervals copied to the caller's buffer */
    int overflow;      /* buffer was too small */
} ClipList;

/* ---- dB helpers ------------------------------------------------------- */

/* linear amplitude -> dBFS; 0 (and NaN) -> -inf */
float gain_db_from_amp(float amp);
/* dBFS -> linear amplitude; -inf -> 0 */
float gain_amp_from_db(float db);

/* dB-domain sums; -inf when every input is silent, NaN for a NaN input. */
float gain_sum_coherent_db(const float* level_db, int n);
float gain_sum_power_db(const float* level_db, int n);

/* ---- bus gain plan ---------------------------------------------------- */

typedef struct {
    int channels;
    float target_dbfs;          /* bus ceiling, dBFS, <= 0 */
    float coherent_sum_dbfs;    /* S_peak of the raw inputs */
    float power_sum_dbfs;       /* S_rms of the raw inputs */
    float bus_headroom_db;      /* target - S_peak; < 0 means over ceiling */
    float rms_margin_db;        /* target - S_rms; headroom left for uncorrelated program */
    float trim_db[GAIN_MAX_CH];       /* applied per-channel gain */
    float share[GAIN_MAX_CH];         /* channel's fraction of the coherent sum */
    float alloc_dbfs[GAIN_MAX_CH];    /* allocated bus slice = post-trim peak */
    float margin_db[GAIN_MAX_CH];     /* 0 dBFS - alloc (dB below full scale) */
    float post_sum_dbfs;        /* coherent sum after trim */
    int silent[GAIN_MAX_CH];    /* channel had no signal (-inf level) */
} GainPlan;

/* Build a plan from per-channel measured peak levels.
 * Returns 0 on success, -1 on invalid arguments:
 *   n outside 1..GAIN_MAX_CH, NULL pointers, NaN/inf(+)/out-of-range level,
 *   target not finite, target > 0 dBFS or < -120 dBFS.
 * Silent channels (-inf level) are left alone (0 dB trim) and contribute 0. */
int gain_plan_init(GainPlan* plan, const float* peak_dbfs, int n, float target_dbfs);

/* ---- meters ----------------------------------------------------------- */

void gain_meter_reset(GainMeter* m);
void gain_meter_feed(GainMeter* m, const float* x, int n);
float gain_meter_peak_db(const GainMeter* m);   /* -inf when empty/silent */
float gain_meter_rms_db(const GainMeter* m);

/* ---- bus summing ------------------------------------------------------ */

/* Sum nch channels into out.
 *   plan != NULL : apply plan->trim_db (grouping enabled)
 *   plan == NULL : unity gains (grouping disabled; this is the legacy path)
 * ramp_frames fades all gains in linearly over the first ramp_frames samples,
 * matching mix_ramp(). Output is hard-clamped to [-1, 1]. meter (optional)
 * measures post-clamp output. spans/list (optional) locate maximal contiguous
 * intervals where the pre-clamp sum exceeded +/-1.0; spans may be NULL only
 * when span_cap == 0, list must be non-NULL if spans is.
 * Returns frames on success, -1 on invalid arguments. */
int gain_bus_sum(const float* const* ch, int nch, int frames,
                 const GainPlan* plan, int ramp_frames, float* out,
                 GainMeter* meter, ClipSpan* spans, int span_cap, ClipList* list);

/* One-shot convenience: measure each channel's peak over the block, derive a
 * plan for target_dbfs, then sum. used_plan (optional) receives the plan.
 * With MIX_GAIN_GROUP=0 this skips planning and sums at unity gain. */
int mix_bus_auto(const float* const* ch, int nch, int frames, float target_dbfs,
                 float* out, GainPlan* used_plan, GainMeter* meter,
                 ClipSpan* spans, int span_cap, ClipList* list);

#endif /* MIX_GAIN_H */
