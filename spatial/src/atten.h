/* atten.h: unified distance attenuation.
 *
 * One parameter set covers the whole distance range.  There is no near/far
 * set switch; the curve is C1 at both seams (d = ref and d = ref + knee):
 *
 *   d <= ref              : g = 1
 *   ref < d <= ref + knee : g = ref / (ref + rolloff * h(d - ref))
 *   d > ref + knee        : g = ref / (ref + rolloff * (d - ref))
 *
 *   h(x) = x^2 * (2*knee - x) / knee^2      (0 <= x <= knee)
 *
 * h(0) = 0 and h'(0) = 0, so the near-field platform (g = 1, g' = 0) joins
 * the knee with matching value and slope.  h(knee) = knee and h'(knee) = 1,
 * so the knee joins the inverse-distance far field with matching value and
 * slope.  With knee = 0 the curve degenerates to C0 (value-continuous, slope
 * kink at d = ref); the default uses a non-zero knee for C1.
 *
 * rolloff = 1 makes the far field exactly ref/d, i.e. the same 1/d law as the
 * linear size of an object on screen.
 *
 * dist == 0 is an error: a source coincident with the listener has no defined
 * attenuation; callers get ATTEN_ERR_DIST_ZERO instead of a silent full-scale
 * gain.  On any error *out_gain is written as 0.0f, never 1.0f.
 */
#ifndef ATTEN_H
#define ATTEN_H

typedef enum {
    ATTEN_OK = 0,
    ATTEN_ERR_DIST_ZERO = -1,   /* dist == 0 (incl. -0.0) */
    ATTEN_ERR_DIST_NEG = -2,    /* dist < 0 */
    ATTEN_ERR_NONFINITE = -3,   /* NaN / inf distance */
    ATTEN_ERR_REF = -4,         /* ref <= 0 or non-finite */
    ATTEN_ERR_ROLLOFF = -5,     /* rolloff < 0 or non-finite */
    ATTEN_ERR_KNEE = -6,        /* knee < 0 or non-finite */
    ATTEN_ERR_ARG = -7,         /* NULL cfg / out pointer missing */
    ATTEN_ERR_SWITCH = -8       /* legacy switch_dist <= ref */
} AttenStatus;

typedef struct {
    float ref;      /* reference distance, > 0 */
    float rolloff;  /* dimensionless falloff, >= 0 */
    float knee;     /* C1 transition width after ref, >= 0 (0 = C0) */
} AttenCfg;

/* The pre-unification two-set model: two rolloff values selected by a hard
 * zone boundary.  Kept only so the evidence tables can show the jump that the
 * unified curve removes.  set_used: 0 = near set, 1 = far set, -1 = error. */
typedef struct {
    float ref;
    float rolloff_near;
    float rolloff_far;
    float switch_dist;  /* near set below switch_dist, far set at/above it */
} AttenLegacyCfg;

/* Single-set inverse-distance law (no switch, no knee).  This is the baseline
 * the legacy model switches between, and the far-field limit of atten_gain. */
AttenStatus atten_gain_set(float dist, float ref, float rolloff, float* out_gain);

/* Unified curve described above. */
AttenStatus atten_gain(float dist, const AttenCfg* cfg, float* out_gain);

/* d(gain)/d(dist) of the unified curve, same domain and status rules. */
AttenStatus atten_gain_slope(float dist, const AttenCfg* cfg, float* out_slope);

/* Legacy two-set curve, for before/after evidence only. */
AttenStatus atten_gain_legacy(float dist, const AttenLegacyCfg* cfg,
                              int* set_used, float* out_gain);

/* 20*log10(gain); -INFINITY for gain <= 0. */
float atten_db(float gain);

const char* atten_status_str(AttenStatus st);

#endif
