
/* lookahead_study: transient peak retention vs limiter lookahead length.
 *
 * Model: audio path delayed by L samples (lookahead), envelope follower
 * reads the undelayed sidechain, gain = envelope applied to delayed audio.
 * No hard clip on output: the clip would pin every peak to the ceiling and
 * hide the effect being measured. Latency introduced = L samples.
 *
 * usage: lookahead_study <out_dir>
 * writes: sweep.csv, peaks.csv, wave_<L>ms.csv
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SR 48000.0f
#define PI_F 3.14159265358979f
#define CEIL 0.8f
#define PEAK_TARGET 1.2f
#define ATTACK_MS 3.0f
#define RELEASE_MS 80.0f
#define MAX_L 961                     /* 20 ms + 1 @ 48 kHz */
#define N 24000                       /* 0.5 s of drums */
#define HIT0_MS 50.0f
#define MEAS_MS 10.0f                 /* transient-peak search window */
#define RETENTION_MIN 0.60f           /* 验收线: 瞬态被压掉不得超过 40% */

static float coef(float ms) { return expf(-1.0f / (ms * 0.001f * SR)); }

/* synthesize two kick hits: pitch-dropping sine + onset click */
static void synth_drums(float* x, int n) {
    const float hits[] = { HIT0_MS, 300.0f };
    for (int h = 0; h < 2; h++) {
        int t0 = (int)(hits[h] * 0.001f * SR);
        /* onset click: damped 3 kHz burst, this is the transient peak */
        for (int i = 0; i < (int)(0.002f * SR) && t0 + i < n; i++) {
            float t = i / SR;
            x[t0 + i] += 1.15f * cosf(2.0f * PI_F * 3000.0f * t) * expf(-t / 0.0006f);
        }
        /* body: pitch-dropping sine, 2 ms ramp in, sustained above ceiling */
        float ph = 0.0f;
        for (int i = 0; i < (int)(0.25f * SR) && t0 + i < n; i++) {
            float t = i / SR;
            float f = 45.0f + 75.0f * expf(-t / 0.03f);
            ph += 2.0f * PI_F * f / SR;
            float ramp = t < 0.002f ? 0.5f - 0.5f * cosf(PI_F * t / 0.002f) : 1.0f;
            x[t0 + i] += 1.05f * ramp * sinf(ph) * expf(-t / 0.4f);
        }
    }
    float pk = 0.0f;
    for (int i = 0; i < n; i++) if (fabsf(x[i]) > pk) pk = fabsf(x[i]);
    for (int i = 0; i < n; i++) x[i] *= PEAK_TARGET / pk;
}

/* lookahead limiter: delay line on audio, envelope on sidechain.
 * Sidechain has a peak detector (instant attack, 5 ms decay) so the
 * envelope target is not duty-cycled by the waveform's own oscillation. */
#define PKDECAY_MS 5.0f
typedef struct { float dly[MAX_L]; int L, idx; float env, pk; } LookLim;

static void ll_init(LookLim* l, int L) {
    memset(l, 0, sizeof(*l));
    l->L = L;
    l->env = 1.0f;
}

static float ll_tick(LookLim* l, float x, float ca, float cr, float cp) {
    float a = fabsf(x);
    l->pk = a > l->pk * cp ? a : l->pk * cp;
    float target = l->pk > CEIL ? CEIL / l->pk : 1.0f;
    float c = target < l->env ? ca : cr;
    l->env = c * l->env + (1.0f - c) * target;
    float d;
    if (l->L == 0) d = x;
    else { d = l->dly[l->idx]; l->dly[l->idx] = x; l->idx = (l->idx + 1) % l->L; }
    return d * l->env;
}

static float peak_in_window(const float* x, int s, int e, int* argmax) {
    float p = 0.0f;
    for (int i = s; i < e; i++) if (fabsf(x[i]) > p) { p = fabsf(x[i]); if (argmax) *argmax = i; }
    return p;
}

int main(int argc, char** argv) {
    const char* dir = argc > 1 ? argv[1] : ".";
    static float in[N], out[N];
    float ca = coef(ATTACK_MS), cr = coef(RELEASE_MS), cp = coef(PKDECAY_MS);
    char path[512];

    synth_drums(in, N);

    int t0 = (int)(HIT0_MS * 0.001f * SR);
    int mw = (int)(MEAS_MS * 0.001f * SR);
    int tp = 0;                       /* index of the transient peak sample */
    float peak_in = peak_in_window(in, t0, t0 + mw, &tp);

    /* sweep lookahead lengths */
    const float sweep[] = { 0.0f, 0.25f, 0.5f, 1.0f, 2.0f, 3.0f, 5.0f, 8.0f, 12.0f, 16.0f };
    const int nsweep = sizeof(sweep) / sizeof(sweep[0]);
    const float rec[] = { 0.5f, 3.0f, 12.0f };   /* 三档录波形 */
    const int nrec = sizeof(rec) / sizeof(rec[0]);

    snprintf(path, sizeof path, "%s/sweep.csv", dir);
    FILE* fs = fopen(path, "w");
    fprintf(fs, "lookahead_ms,latency_ms,peak_in,peak_out,retention\n");

    snprintf(path, sizeof path, "%s/peaks.csv", dir);
    FILE* fp = fopen(path, "w");
    fprintf(fp, "lookahead_ms,latency_ms,peak_in,peak_out,retention,crushed_pct,pass\n");

    int all_pass = 1;
    for (int s = 0; s < nsweep; s++) {
        int L = (int)(sweep[s] * 0.001f * SR + 0.5f);
        LookLim l; ll_init(&l, L);
        for (int i = 0; i < N; i++) out[i] = ll_tick(&l, in[i], ca, cr, cp);
        /* track the transient peak sample itself through the delay line:
         * in[tp] exits at tp+L with gain env(tp+L) */
        float po = fabsf(out[tp + L]);
        float r = po / peak_in;
        float lat_ms = 1000.0f * L / SR;
        fprintf(fs, "%.2f,%.3f,%.4f,%.4f,%.4f\n", sweep[s], lat_ms, peak_in, po, r);

        int is_rec = 0;
        for (int k = 0; k < nrec; k++)
            if (fabsf(sweep[s] - rec[k]) < 1e-6f) is_rec = 1;
        if (is_rec) {
            int pass = r >= RETENTION_MIN;
            if (!pass) all_pass = 0;
            fprintf(fp, "%.2f,%.3f,%.4f,%.4f,%.4f,%.1f,%s\n",
                    sweep[s], lat_ms, peak_in, po, r, 100.0f * (1.0f - r),
                    pass ? "PASS" : "FAIL");
            snprintf(path, sizeof path, "%s/wave_%.1fms.csv", dir, sweep[s]);
            FILE* fw = fopen(path, "w");
            fprintf(fw, "t_ms,input,output\n");
            int a = t0 - (int)(5.0f * 0.001f * SR);
            int b = t0 + (int)(50.0f * 0.001f * SR);
            for (int i = a; i < b; i++)
                fprintf(fw, "%.3f,%.4f,%.4f\n", 1000.0f * i / SR, in[i], out[i]);
            fclose(fw);
        }
    }
    fclose(fs);
    fclose(fp);

    printf("peak_in=%.4f ceiling=%.2f retention_min=%.2f\n", peak_in, CEIL, RETENTION_MIN);
    printf("all_recorded_lengths_pass=%d\n", all_pass);
    return all_pass ? 0 : 1;
}
