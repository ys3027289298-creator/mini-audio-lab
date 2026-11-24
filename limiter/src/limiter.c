
/* limiter: self-test for limiter.h. */
#include <stdio.h>
#include <math.h>
#include "limiter.h"

#define BLK 128

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Limiter l;
    float in[BLK], out[BLK];
    Report rep;

    for (int i = 0; i < BLK; i++) in[i] = 0.5f;
    lim_init(&l, 5.0f, 50.0f, 48000.0f);
    check("quiet_passthrough", lim_process(&l, in, out, BLK, 0.9f, &rep) == BLK && rep.limited_samples == 0);
    check("quiet_peak", fabsf(rep.peak_out - 0.5f) < 1e-4f);

    for (int i = 0; i < BLK; i++) in[i] = 2.0f;
    lim_init(&l, 5.0f, 50.0f, 48000.0f);
    lim_process(&l, in, out, BLK, 0.8f, &rep);
    check("loud_limited", rep.limited_samples > 0);
    check("below_ceiling", rep.peak_out <= 0.8f + 1e-4f);
    check("reduction_positive", rep.max_reduction > 0.0f);

    check("bad_ceiling", lim_process(&l, in, out, BLK, 0.0f, &rep) == -1);
    lim_init(&l, 0.0f, 50.0f, 48000.0f);
    check("bad_attack", lim_process(&l, in, out, BLK, 0.8f, &rep) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
