
/* sidechain: key-input envelope trigger with attack/release follower and hold state machine.
 *
 * Two build modes in one binary:
 *   no args  -> built-in self test (used by `make test`)
 *   args     -> stream mode:
 *       sidechain in.f32 sr thresh_db attack_ms release_ms hold_ms env_out.f32 gate_out.u8
 *     in.f32 is raw mono float32. env_out gets the per-sample envelope (float32),
 *     gate_out gets the per-sample gate state (uint8). Stats are printed to stdout.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct {
    float sample_rate;
    float attack_ms;
    float release_ms;
    float hold_ms;
    float threshold;      /* linear amplitude */
    float env;
    int engaged;
    int hold_samples;
    int hold_left;
    long events;          /* number of upward threshold crossings */
} Sidechain;

static int sc_init(Sidechain* s, float sr, float attack_ms, float release_ms,
                   float hold_ms, float threshold) {
    if (s == NULL || sr <= 0.0f || attack_ms <= 0.0f || release_ms <= 0.0f ||
        hold_ms < 0.0f || threshold <= 0.0f)
        return -1;
    memset(s, 0, sizeof(*s));
    s->sample_rate = sr;
    s->attack_ms = attack_ms;
    s->release_ms = release_ms;
    s->hold_ms = hold_ms;
    s->threshold = threshold;
    s->hold_samples = (int)(hold_ms * 0.001f * sr + 0.5f);
    s->env = 0.0f;
    s->engaged = 0;
    return 0;
}

/* one-pole coefficient from a time constant */
static float sc_coef(float ms, float sr) {
    return expf(-1.0f / (ms * 0.001f * sr));
}

/* process one sample; returns 1 when a new trigger (upward crossing) occurs */
static int sc_step(Sidechain* s, float x) {
    float ca = sc_coef(s->attack_ms, s->sample_rate);
    float cr = sc_coef(s->release_ms, s->sample_rate);
    float ax = fabsf(x);
    /* fast attack, slow release envelope follower */
    float c = ax > s->env ? ca : cr;
    s->env = c * s->env + (1.0f - c) * ax;

    int new_trig = 0;
    if (!s->engaged) {
        if (s->env >= s->threshold) {
            s->engaged = 1;
            s->hold_left = s->hold_samples;
            s->events++;
            new_trig = 1;
        }
    } else {
        if (s->env >= s->threshold) {
            s->hold_left = s->hold_samples;   /* re-arm hold while above */
        } else if (s->hold_left > 0) {
            s->hold_left--;                   /* stay engaged during hold */
        } else {
            s->engaged = 0;
        }
    }
    return new_trig;
}

static int sc_process(Sidechain* s, const float* in, int n, float* env_out,
                      unsigned char* gate_out) {
    if (s == NULL || in == NULL || n <= 0) return -1;
    long engaged = 0;
    for (int i = 0; i < n; i++) {
        sc_step(s, in[i]);
        if (env_out) env_out[i] = s->env;
        if (gate_out) gate_out[i] = (unsigned char)s->engaged;
        if (s->engaged) engaged++;
    }
    return (int)engaged;
}

static float* read_f32(const char* path, long* n_out) {
    FILE* f = fopen(path, "rb");
    if (!f) { perror(path); return NULL; }
    fseek(f, 0, SEEK_END);
    long bytes = ftell(f);
    fseek(f, 0, SEEK_SET);
    long n = bytes / (long)sizeof(float);
    float* buf = (float*)malloc((size_t)bytes + 1);
    if (!buf) { fclose(f); return NULL; }
    if (fread(buf, 1, (size_t)bytes, f) != (size_t)bytes) {
        free(buf); fclose(f); return NULL;
    }
    fclose(f);
    *n_out = n;
    return buf;
}

static int write_bytes(const char* path, const void* buf, size_t bytes) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return -1; }
    int ok = fwrite(buf, 1, bytes, f) == bytes ? 0 : -1;
    fclose(f);
    return ok;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void selftest(void) {
    enum { N = 4800 };           /* 100 ms at 48 kHz */
    float in[N], env[N];
    unsigned char gate[N];
    Sidechain s;

    /* silence must never trigger */
    for (int i = 0; i < N; i++) in[i] = 0.0f;
    check("init_ok", sc_init(&s, 48000.0f, 2.0f, 20.0f, 30.0f, 0.1f) == 0);
    sc_process(&s, in, N, env, gate);
    check("silence_no_trig", s.events == 0 && s.engaged == 0);

    /* loud step must trigger exactly once and stay engaged */
    for (int i = 0; i < N; i++) in[i] = 0.5f;
    sc_init(&s, 48000.0f, 2.0f, 20.0f, 30.0f, 0.1f);
    sc_process(&s, in, N, env, gate);
    check("loud_triggers", s.events == 1);
    check("loud_engaged", s.engaged == 1);
    check("envelope_tracks", fabsf(env[N - 1] - 0.5f) < 1e-2f);

    /* dip shorter than hold keeps the gate engaged; long dip releases */
    sc_init(&s, 48000.0f, 2.0f, 20.0f, 30.0f, 0.1f);
    for (int i = 0; i < N; i++) in[i] = 0.5f;
    sc_process(&s, in, 2400, env, gate);                 /* 50 ms on */
    int was_engaged = s.engaged;
    for (int i = 0; i < 480; i++) in[i] = 0.0f;         /* 10 ms dip < 30 ms hold */
    sc_process(&s, in, 480, env, gate);
    check("hold_keeps_gate", was_engaged == 1 && s.engaged == 1);
    for (int i = 0; i < N; i++) in[i] = 0.0f;
    sc_process(&s, in, N, env, gate);                    /* long quiet stretch */
    check("long_dip_releases", s.engaged == 0);

    /* re-trigger after release counts a second event */
    long before = s.events;
    for (int i = 0; i < N; i++) in[i] = 0.5f;
    sc_process(&s, in, N, env, gate);
    check("retrig_counts", s.events == before + 1);

    /* parameter validation */
    check("bad_sr", sc_init(&s, 0.0f, 2.0f, 20.0f, 30.0f, 0.1f) == -1);
    check("bad_attack", sc_init(&s, 48000.0f, 0.0f, 20.0f, 30.0f, 0.1f) == -1);
    check("bad_hold", sc_init(&s, 48000.0f, 2.0f, 20.0f, -1.0f, 0.1f) == -1);
    check("bad_threshold", sc_init(&s, 48000.0f, 2.0f, 20.0f, 30.0f, 0.0f) == -1);
    check("bad_process", sc_process(&s, NULL, N, env, gate) == -1);

    printf("all_ok=%d\n", failures == 0);
}

int main(int argc, char** argv) {
    if (argc == 1) {
        selftest();
        return failures == 0 ? 0 : 1;
    }
    if (argc != 9) {
        fprintf(stderr,
                "usage: %s in.f32 sr thresh_db attack_ms release_ms hold_ms env_out.f32 gate_out.u8\n",
                argv[0]);
        return 2;
    }
    float sr = (float)atof(argv[2]);
    float thresh_db = (float)atof(argv[3]);
    float attack_ms = (float)atof(argv[4]);
    float release_ms = (float)atof(argv[5]);
    float hold_ms = (float)atof(argv[6]);
    float threshold = powf(10.0f, thresh_db / 20.0f);

    long n = 0;
    float* in = read_f32(argv[1], &n);
    if (!in) return 2;
    float* env = (float*)malloc((size_t)n * sizeof(float));
    unsigned char* gate = (unsigned char*)malloc((size_t)n);
    if (!env || !gate) { free(in); free(env); free(gate); return 2; }

    Sidechain s;
    if (sc_init(&s, sr, attack_ms, release_ms, hold_ms, threshold) != 0) {
        fprintf(stderr, "bad parameters\n");
        return 2;
    }
    long engaged = sc_process(&s, in, (int)n, env, gate);
    if (write_bytes(argv[7], env, (size_t)n * sizeof(float)) != 0 ||
        write_bytes(argv[8], gate, (size_t)n) != 0)
        return 2;

    printf("samples=%ld\n", n);
    printf("threshold_db=%.3f\n", thresh_db);
    printf("trigger_events=%ld\n", s.events);
    printf("engaged_samples=%ld\n", engaged);
    printf("all_ok=1\n");

    free(in); free(env); free(gate);
    return 0;
}
