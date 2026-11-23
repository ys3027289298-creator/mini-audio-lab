/* atten.c: unified distance attenuation, C1 across the whole range. */
#include "atten.h"
#include <math.h>
#include <stddef.h>

/* knee easing and its first derivative on [0, w] */
static double knee_h(double x, double w) {
    return x * x * (2.0 * w - x) / (w * w);
}

static double knee_hp(double x, double w) {
    return 4.0 * x / w - 3.0 * x * x / (w * w);
}

static AttenStatus check_cfg(const AttenCfg* c) {
    if (!isfinite(c->ref) || c->ref <= 0.0f) return ATTEN_ERR_REF;
    if (!isfinite(c->rolloff) || c->rolloff < 0.0f) return ATTEN_ERR_ROLLOFF;
    if (!isfinite(c->knee) || c->knee < 0.0f) return ATTEN_ERR_KNEE;
    return ATTEN_OK;
}

/* NaN first: NaN is neither == 0 nor < 0. -0.0f compares equal to 0.0f. */
static AttenStatus check_dist(float d) {
    if (!isfinite(d)) return ATTEN_ERR_NONFINITE;
    if (d == 0.0f) return ATTEN_ERR_DIST_ZERO;
    if (d < 0.0f) return ATTEN_ERR_DIST_NEG;
    return ATTEN_OK;
}

AttenStatus atten_gain_set(float dist, float ref, float rolloff, float* out_gain) {
    if (out_gain) *out_gain = 0.0f;
    if (!out_gain) return ATTEN_ERR_ARG;
    AttenCfg c = { ref, rolloff, 0.0f };
    AttenStatus st = check_cfg(&c);
    if (st != ATTEN_OK) return st;
    st = check_dist(dist);
    if (st != ATTEN_OK) return st;
    double r = ref, k = rolloff, d = dist;
    *out_gain = d <= r ? 1.0f : (float)(r / (r + k * (d - r)));
    return ATTEN_OK;
}

AttenStatus atten_gain(float dist, const AttenCfg* cfg, float* out_gain) {
    if (out_gain) *out_gain = 0.0f;
    if (!out_gain || !cfg) return ATTEN_ERR_ARG;
    AttenStatus st = check_cfg(cfg);
    if (st != ATTEN_OK) return st;
    st = check_dist(dist);
    if (st != ATTEN_OK) return st;

    double r = cfg->ref, k = cfg->rolloff, w = cfg->knee, d = dist;
    double g;
    if (d <= r) {
        g = 1.0;
    } else if (w > 0.0 && d < r + w) {
        g = r / (r + k * knee_h(d - r, w));
    } else {
        g = r / (r + k * (d - r));
    }
    *out_gain = (float)g;
    return ATTEN_OK;
}

AttenStatus atten_gain_slope(float dist, const AttenCfg* cfg, float* out_slope) {
    if (out_slope) *out_slope = 0.0f;
    if (!out_slope || !cfg) return ATTEN_ERR_ARG;
    AttenStatus st = check_cfg(cfg);
    if (st != ATTEN_OK) return st;
    st = check_dist(dist);
    if (st != ATTEN_OK) return st;

    double r = cfg->ref, k = cfg->rolloff, w = cfg->knee, d = dist;
    double s;
    if (d <= r) {
        s = 0.0;
    } else if (w > 0.0 && d < r + w) {
        double den = r + k * knee_h(d - r, w);
        s = -r * k * knee_hp(d - r, w) / (den * den);
    } else {
        double den = r + k * (d - r);
        s = -r * k / (den * den);
    }
    *out_slope = (float)s;
    return ATTEN_OK;
}

AttenStatus atten_gain_legacy(float dist, const AttenLegacyCfg* cfg,
                              int* set_used, float* out_gain) {
    if (out_gain) *out_gain = 0.0f;
    if (set_used) *set_used = -1;
    if (!out_gain || !cfg) return ATTEN_ERR_ARG;
    if (!isfinite(cfg->ref) || cfg->ref <= 0.0f) return ATTEN_ERR_REF;
    if (!isfinite(cfg->rolloff_near) || cfg->rolloff_near < 0.0f)
        return ATTEN_ERR_ROLLOFF;
    if (!isfinite(cfg->rolloff_far) || cfg->rolloff_far < 0.0f)
        return ATTEN_ERR_ROLLOFF;
    if (!isfinite(cfg->switch_dist) || cfg->switch_dist <= cfg->ref)
        return ATTEN_ERR_SWITCH;
    AttenStatus st = check_dist(dist);
    if (st != ATTEN_OK) return st;

    double r = cfg->ref, sw = cfg->switch_dist, d = dist;
    double g;
    int set;
    if (d <= r) {
        g = 1.0;
        set = 0;
    } else if (d < sw) {
        double k = cfg->rolloff_near;
        g = r / (r + k * (d - r));
        set = 0;
    } else {
        double k = cfg->rolloff_far;
        g = r / (r + k * (d - r));
        set = 1;
    }
    *out_gain = (float)g;
    if (set_used) *set_used = set;
    return ATTEN_OK;
}

float atten_db(float gain) {
    return gain > 0.0f ? 20.0f * log10f(gain) : -INFINITY;
}

const char* atten_status_str(AttenStatus st) {
    switch (st) {
    case ATTEN_OK: return "ATTEN_OK";
    case ATTEN_ERR_DIST_ZERO: return "ATTEN_ERR_DIST_ZERO";
    case ATTEN_ERR_DIST_NEG: return "ATTEN_ERR_DIST_NEG";
    case ATTEN_ERR_NONFINITE: return "ATTEN_ERR_NONFINITE";
    case ATTEN_ERR_REF: return "ATTEN_ERR_REF";
    case ATTEN_ERR_ROLLOFF: return "ATTEN_ERR_ROLLOFF";
    case ATTEN_ERR_KNEE: return "ATTEN_ERR_KNEE";
    case ATTEN_ERR_ARG: return "ATTEN_ERR_ARG";
    case ATTEN_ERR_SWITCH: return "ATTEN_ERR_SWITCH";
    default: return "ATTEN_ERR_UNKNOWN";
    }
}
