
/* limiter: peak limiter with attack/release envelope and gain reduction report. */
#ifndef MINI_LIMITER_H
#define MINI_LIMITER_H

#include <math.h>
#include <string.h>

typedef struct {
    float attack_ms;
    float release_ms;
    float sample_rate;
    float env;
    int engaged;
} Limiter;

typedef struct { float peak_in; float peak_out; int limited_samples; float max_reduction; } Report;

static void lim_init(Limiter* l, float attack_ms, float release_ms, float sr) {
    memset(l, 0, sizeof(*l));
    l->attack_ms = attack_ms;
    l->release_ms = release_ms;
    l->sample_rate = sr;
    l->env = 1.0f;
}

/* one-pole coefficients from time constants */
static float lim_coef(float ms, float sr) {
    if (ms <= 0.0f || sr <= 0.0f) return 0.0f;
    return expf(-1.0f / (ms * 0.001f * sr));
}

static int lim_process(Limiter* l, const float* in, float* out, int n, float ceiling, Report* rep) {
    if (l == NULL || n <= 0) return -1;
    if (ceiling <= 0.0f || ceiling > 1.0f) return -1;
    if (l->attack_ms <= 0.0f || l->release_ms <= 0.0f) return -1;

    float ca = lim_coef(l->attack_ms, l->sample_rate);
    float cr = lim_coef(l->release_ms, l->sample_rate);
    float pin = 0.0f, pout = 0.0f, red = 0.0f;
    int lim = 0;

    for (int i = 0; i < n; i++) {
        float x = in[i];
        if (fabsf(x) > pin) pin = fabsf(x);
        float target = fabsf(x) > ceiling ? ceiling / fabsf(x) : 1.0f;
        /* attack when target below env, release otherwise */
        float c = target < l->env ? ca : cr;
        l->env = c * l->env + (1.0f - c) * target;
        if (l->env < 1.0f) { lim++; float r = 1.0f - l->env; if (r > red) red = r; }
        float y = x * l->env;
        if (y > ceiling) y = ceiling;
        if (y < -ceiling) y = -ceiling;
        out[i] = y;
        if (fabsf(y) > pout) pout = fabsf(y);
    }
    if (rep) { rep->peak_in = pin; rep->peak_out = pout; rep->limited_samples = lim; rep->max_reduction = red; }
    l->engaged = lim > 0;
    return n;
}

#endif
