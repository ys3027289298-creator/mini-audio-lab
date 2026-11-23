
/* meter: unified level metering shared by export and playback paths. */
#include <stdio.h>
#include <math.h>
#include "meter.h"

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* near-full-scale sine: peak ~0 dBFS, rms ~-3.01 dBFS */
    float x[4800];
    for (int i = 0; i < 4800; i++)
        x[i] = 0.9999f * sinf(2.0f * 3.14159265358979f * 440.0f * (float)i / 48000.0f);

    Meter a, b;
    meter_reset(&a);
    meter_feed(&a, x, 4800);
    check("sine_peak_0db", fabsf(meter_peak_dbfs(&a)) < 0.01f);
    check("sine_rms_minus3db", fabsf(meter_rms_dbfs(&a) - -3.0103f) < 0.01f);
    check("sine_not_clipped", a.clipped == 0);

    /* block-by-block feed must equal one-shot feed (统一窗口与分块无关, 含尾块) */
    meter_reset(&b);
    for (int off = 0; off < 4800; off += 128) {
        int m = 4800 - off < 128 ? 4800 - off : 128;
        meter_feed(&b, x + off, m);
    }
    check("block_eq_oneshot_peak", meter_peak_dbfs(&a) == meter_peak_dbfs(&b));
    check("block_eq_oneshot_rms", fabsf(meter_rms_dbfs(&a) - meter_rms_dbfs(&b)) < 1e-5f);

    /* dc 0.5: peak == rms == -6.02 dB */
    float dc[256];
    for (int i = 0; i < 256; i++) dc[i] = 0.5f;
    Meter c;
    meter_reset(&c);
    meter_feed(&c, dc, 256);
    check("dc_peak_eq_rms", fabsf(meter_peak_dbfs(&c) - meter_rms_dbfs(&c)) < 1e-5f);
    check("dc_level", fabsf(meter_rms_dbfs(&c) - -6.0206f) < 0.01f);

    /* empty meter: -inf, no crash */
    Meter e;
    meter_reset(&e);
    check("empty_peak_ninf", isinf(meter_peak_dbfs(&e)) && meter_peak_dbfs(&e) < 0.0f);
    check("empty_rms_ninf", isinf(meter_rms_dbfs(&e)) && meter_rms_dbfs(&e) < 0.0f);

    /* clipping counter */
    float hot[4] = { 1.0f, -1.0f, 0.5f, 1.2f };
    Meter d;
    meter_reset(&d);
    meter_feed(&d, hot, 4);
    check("clip_count", d.clipped == 3);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
