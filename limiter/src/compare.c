
/* compare: old single-rate release vs new envelope-split release on the same
 * material, plus a pure-impulse overshoot check. Writes CSVs into out/ for
 * tools/plot_compare.py. The old algorithm is frozen here on purpose so the
 * "before" curve does not drift with limiter.c. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SR 48000.0f
#define CEIL 0.8f

#define ATTACK_MS 5.0f
#define OLD_RELEASE_MS 300.0f   /* old: one release rate for everything */
#define NEW_FAST_MS 60.0f      /* new: fast while the level envelope is hot */
#define NEW_SLOW_MS 400.0f     /* new: slow while it is quiet */
#define NEW_SPLIT 0.25f
#define LEVEL_ATTACK_MS 1.0f
#define LEVEL_RELEASE_MS 200.0f
#define BLEND_MS 10.0f

static float coef(float ms) {
    return expf(-1.0f / (ms * 0.001f * (float)SR));
}

/* old limiter: single release rate, gain envelope recorded per sample */
static void old_process(const float* in, float* gain, int n) {
    float ca = coef(ATTACK_MS), cr = coef(OLD_RELEASE_MS);
    float env = 1.0f;
    for (int i = 0; i < n; i++) {
        float ax = fabsf(in[i]);
        float target = ax > CEIL ? CEIL / ax : 1.0f;
        float c = target < env ? ca : cr;
        env = c * env + (1.0f - c) * target;
        gain[i] = env;
    }
}

/* new limiter: envelope-split release, mirrors src/limiter.c */
static void new_process(const float* in, float* gain, int n) {
    float ca = coef(ATTACK_MS), cf = coef(NEW_FAST_MS), cs = coef(NEW_SLOW_MS);
    float cw = coef(BLEND_MS), cla = coef(LEVEL_ATTACK_MS), clr = coef(LEVEL_RELEASE_MS);
    float env = 1.0f, level = 0.0f, w = 0.0f;
    for (int i = 0; i < n; i++) {
        float ax = fabsf(in[i]);
        float cl = ax > level ? cla : clr;
        level = cl * level + (1.0f - cl) * ax;
        float target = ax > CEIL ? CEIL / ax : 1.0f;
        if (target < env) {
            env = ca * env + (1.0f - ca) * target;
        } else {
            float wt = level >= NEW_SPLIT ? 1.0f : 0.0f;
            w = cw * w + (1.0f - cw) * wt;
            float cr = w * cf + (1.0f - w) * cs;
            env = cr * env + (1.0f - cr) * target;
        }
        gain[i] = env;
    }
}

/* material: 0.4 s loud burst (over ceiling) then 1.6 s quiet passage -
 * the case where the old limiter drags accumulated reduction into the
 * quiet part and sounds muddy */
static float material(int i, unsigned* seed) {
    float t = (float)i / (float)SR;
    if (t < 0.4f) return 1.2f * sinf(2.0f * (float)M_PI * 110.0f * t);
    *seed = *seed * 1664525u + 1013904223u;
    float noise = (float)(*seed >> 8) / 8388608.0f - 1.0f;
    return 0.06f * sinf(2.0f * (float)M_PI * 440.0f * t) + 0.01f * noise;
}

int main(void) {
    const int n = (int)(2.0 * SR);
    float* in = malloc(sizeof(float) * n);
    float* g_old = malloc(sizeof(float) * n);
    float* g_new = malloc(sizeof(float) * n);
    if (!in || !g_old || !g_new) return 1;

    unsigned seed = 22222u;
    for (int i = 0; i < n; i++) in[i] = material(i, &seed);
    old_process(in, g_old, n);
    new_process(in, g_new, n);

    FILE* f = fopen("out/material_gain.csv", "w");
    if (!f) return 1;
    fprintf(f, "t,in,gain_old,gain_new\n");
    float red_old = 0.0f, red_new = 0.0f;
    for (int i = 0; i < n; i++) {
        fprintf(f, "%.6f,%.6f,%.6f,%.6f\n", (float)i / (float)SR, in[i], g_old[i], g_new[i]);
        if (1.0f - g_old[i] > red_old) red_old = 1.0f - g_old[i];
        if (1.0f - g_new[i] > red_new) red_new = 1.0f - g_new[i];
    }
    fclose(f);
    printf("material_max_reduction_old=%.4f\n", red_old);
    printf("material_max_reduction_new=%.4f\n", red_new);

    /* pure impulse: one full-scale sample at 10 ms, silence elsewhere */
    const int ni = (int)(0.2 * SR);
    for (int i = 0; i < ni; i++) in[i] = 0.0f;
    in[(int)(0.01 * SR)] = 1.0f;
    old_process(in, g_old, ni);
    new_process(in, g_new, ni);

    f = fopen("out/impulse_gain.csv", "w");
    if (!f) return 1;
    fprintf(f, "t,gain_old,gain_new\n");
    float max_old = 0.0f, max_new = 0.0f, min_old = 1.0f, min_new = 1.0f;
    for (int i = 0; i < ni; i++) {
        fprintf(f, "%.6f,%.6f,%.6f\n", (float)i / (float)SR, g_old[i], g_new[i]);
        if (g_old[i] > max_old) max_old = g_old[i];
        if (g_new[i] > max_new) max_new = g_new[i];
        if (g_old[i] < min_old) min_old = g_old[i];
        if (g_new[i] < min_new) min_new = g_new[i];
    }
    fclose(f);
    printf("impulse_gain_range_old=[%.6f,%.6f]\n", min_old, max_old);
    printf("impulse_gain_range_new=[%.6f,%.6f]\n", min_new, max_new);
    printf("impulse_overshoot_old=%d\n", max_old > 1.0f + 1e-6f);
    printf("impulse_overshoot_new=%d\n", max_new > 1.0f + 1e-6f);

    free(in); free(g_old); free(g_new);
    return 0;
}
