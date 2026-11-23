/* dither_export: render test signals and export them as 16-bit PCM WAV,
   with or without TPDF dither. Alongside each WAV it writes a .f32 sidecar
   with the exact float input, so an independent verifier can check the
   dither-off path byte for byte. Output is reproducible across runs. */
#include "dither.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FS 48000

/* ---- WAV writer: 16-bit PCM, mono ---- */
static void put_u32(FILE* f, uint32_t v) {
    fputc(v & 0xff, f); fputc((v >> 8) & 0xff, f);
    fputc((v >> 16) & 0xff, f); fputc((v >> 24) & 0xff, f);
}
static void put_u16(FILE* f, uint16_t v) { fputc(v & 0xff, f); fputc((v >> 8) & 0xff, f); }

static int write_wav16(const char* path, const int16_t* codes, int n) {
    FILE* f = fopen(path, "wb");
    if (!f) return -1;
    uint32_t data_bytes = (uint32_t)n * 2u;
    fwrite("RIFF", 1, 4, f); put_u32(f, 36u + data_bytes); fwrite("WAVE", 1, 4, f);
    fwrite("fmt ", 1, 4, f); put_u32(f, 16u); put_u16(f, 1); put_u16(f, 1);
    put_u32(f, FS); put_u32(f, FS * 2u); put_u16(f, 2); put_u16(f, 16);
    fwrite("data", 1, 4, f); put_u32(f, data_bytes);
    for (int i = 0; i < n; i++) put_u16(f, (uint16_t)codes[i]);
    fclose(f);
    return 0;
}

/* ---- Sidecar: the exact float input, raw little-endian float32, so the
   independent verifier can re-derive the dither-free codes byte for byte. ---- */
static int write_f32(const char* wav_path, const float* x, int n) {
    char side[1024];
    if (snprintf(side, sizeof(side), "%s.f32", wav_path) >= (int)sizeof(side))
        return -1;
    FILE* f = fopen(side, "wb");
    if (!f) return -1;
    for (int i = 0; i < n; i++) {
        uint32_t u;
        memcpy(&u, &x[i], 4);
        put_u32(f, u);
    }
    fclose(f);
    return 0;
}

/* ---- Test signals ---- */

/* Sine at an exact frequency, amplitude given in dBFS. 3 s. */
static void gen_tone(float* x, int n, float freq, float dbfs) {
    float a = powf(10.0f, dbfs / 20.0f);
    for (int i = 0; i < n; i++)
        x[i] = a * sinf(2.0f * 3.14159265f * freq * (float)i / FS);
}

/* Sine decaying linearly -6 -> -90 dBFS over the buffer: the classic
   "dither is level-independent" demonstration in one file. */
static void gen_sweep(float* x, int n, float freq, float db_hi, float db_lo) {
    for (int i = 0; i < n; i++) {
        float frac = (float)i / (float)(n - 1);
        float db = db_hi + (db_lo - db_hi) * frac;
        x[i] = powf(10.0f, db / 20.0f) * sinf(2.0f * 3.14159265f * freq * (float)i / FS);
    }
}

/* Tone gated on/off (100 ms each) at -70 dBFS: shows the "silence" that
   plain rounding turns into hard digital zero, and the noise that dither
   leaves behind so the gaps do not read as dropouts. */
static void gen_gated(float* x, int n, float freq, float dbfs) {
    float a = powf(10.0f, dbfs / 20.0f);
    for (int i = 0; i < n; i++) {
        int phase = (i / (FS / 10)) % 2;   /* 100 ms blocks */
        x[i] = phase ? a * sinf(2.0f * 3.14159265f * freq * (float)i / FS) : 0.0f;
    }
}

/* ---- CLI ---- */
static int usage(const char* argv0) {
    fprintf(stderr,
        "usage: %s <mode> <out.wav> [freq_hz]\n"
        "  mode: tone60   -60 dBFS, 3 s\n"
        "        sweep    -6 -> -90 dBFS, 3 s\n"
        "        gated    -70 dBFS gated 100/100 ms, 3 s\n"
        "  freq_hz defaults to 750, which is exactly 64 samples per cycle at\n"
        "  48 kHz, so it is coherent for any power-of-two FFT and harmonic and\n"
        "  signal metrics land on exact bins instead of smearing with leakage.\n"
        "  DITHER=on|off is read from the environment (default on),\n"
        "  DITHER_SEED=<uint64> sets the stream seed (default 1234), and\n"
        "  TONE_DB=<level> overrides the tone level for tone60 (default -60).\n",
        argv0);
    return 2;
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) return usage(argv[0]);
    const char* mode = argv[1];
    const char* path = argv[2];
    float freq = argc == 4 ? (float)atof(argv[3]) : 750.0f;
    if (!(freq > 0.0f && freq < (float)FS / 2.0f)) return usage(argv[0]);

    int n = FS * 3;
    float* x = (float*)malloc(sizeof(float) * (size_t)n);
    int16_t* codes = (int16_t*)malloc(sizeof(int16_t) * (size_t)n);
    if (!x || !codes) { free(x); free(codes); return 1; }

    float tone_db = -60.0f;
    const char* lv = getenv("TONE_DB");
    if (lv) tone_db = (float)atof(lv);
    if (strcmp(mode, "tone60") == 0)       gen_tone(x, n, freq, tone_db);
    else if (strcmp(mode, "sweep") == 0)   gen_sweep(x, n, freq, -6.0f, -90.0f);
    else if (strcmp(mode, "gated") == 0)   gen_gated(x, n, freq, -70.0f);
    else { usage(argv[0]); free(x); free(codes); return 2; }

    const char* d = getenv("DITHER");
    int dither = !(d && strcmp(d, "off") == 0);
    uint64_t seed = 1234;
    const char* s = getenv("DITHER_SEED");
    if (s) seed = strtoull(s, NULL, 10);

    DitherRng rng;
    dither_rng_seed(&rng, seed, 0x5eed);
    DitherReport rep;
    if (dither_quantize16_code(x, n, codes, dither, &rng, &rep) != n) {
        fprintf(stderr, "quantize failed\n"); free(x); free(codes); return 1;
    }
    if (write_wav16(path, codes, n) != 0) {
        fprintf(stderr, "cannot write %s\n", path); free(x); free(codes); return 1;
    }
    if (write_f32(path, x, n) != 0) {
        fprintf(stderr, "cannot write %s.f32\n", path); free(x); free(codes); return 1;
    }
    printf("mode=%s freq=%.1f dither=%d seed=%llu peak=%.6f err_rms=%.4f "
           "err_dc=%+.5f clipped=%lld changed=%lld\n",
           mode, (double)freq, dither, (unsigned long long)seed, rep.peak,
           rep.err_rms, rep.err_dc, (long long)rep.clipped,
           (long long)rep.changed);
    free(x); free(codes);
    return 0;
}
