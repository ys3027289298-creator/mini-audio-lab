
/* bench_blockfft: block-length vs. error/latency trade-off measurement.
 *
 * For each sample rate and candidate block length it streams a fixed
 * multi-tone + noise signal through the block analyzer and compares the
 * estimated one-sided PSD against the full-segment (blocking disabled)
 * reference on a common frequency grid. Emits one CSV row per run. */
#include "blockfft.c"
#include <time.h>

/* deterministic LCG noise so runs are reproducible */
static float lcg_noise(unsigned* s) {
    *s = *s * 1664525u + 1013904223u;
    return ((float)(*s >> 8) / (float)(1u << 24)) - 0.5f;
}

/* fixed test signal: three tones (including a close pair for resolution),
 * plus low-level white noise; n samples at rate sr */
static void gen_signal(float* x, int n, int sr) {
    unsigned seed = 20261124u;
    for (int i = 0; i < n; i++) {
        double t = (double)i / sr;
        double v = 0.6 * cos(2.0 * M_PI * 1000.0 * t)
                 + 0.4 * cos(2.0 * M_PI * 2000.0 * t)
                 + 0.4 * cos(2.0 * M_PI * 2350.0 * t);   /* 350 Hz apart */
        v += 0.02 * (double)lcg_noise(&seed);            /* ~ -46 dBFS noise floor */
        x[i] = (float)v;
    }
}

#if BLOCKING
/* fraction of bin k's cell [k-0.5, k+0.5]*df overlapping [flo, fhi] */
static double cell_overlap(int k, double df, double flo, double fhi) {
    double lo = (k - 0.5) * df, hi = (k + 0.5) * df;
    double a = lo > flo ? lo : flo;
    double b = hi < fhi ? hi : fhi;
    if (b <= a) return 0.0;
    return (b - a) / df;
}

/* integrated band power (PSD * df) in [flo, fhi], grid agnostic */
static double band_power(const double* psd, int n, double df, double flo, double fhi) {
    double s = 0.0;
    int k0 = (int)(flo / df - 0.5), k1 = (int)(fhi / df + 0.5);
    if (k0 < 0) k0 = 0;
    if (k1 > n) k1 = n;
    for (int k = k0; k <= k1; k++)
        s += psd[k] * df * cell_overlap(k, df, flo, fhi);
    return s;
}

/* RMS band-power error (dB) vs. the full-segment reference over 250 Hz
 * bands up to 8 kHz; converges to 0 as B approaches the segment length */
static double band_error_db(const double* psd, int B, int sr,
                            const double* ref, int N) {
    const double bw = 250.0, fmax = 8000.0;
    int nbands = (int)(fmax / bw);
    double err2 = 0.0;
    double df_b = (double)sr / B, df_r = (double)sr / N;
    for (int b = 0; b < nbands; b++) {
        double pb = band_power(psd, B / 2, df_b, b * bw, (b + 1) * bw);
        double pr = band_power(ref, N / 2, df_r, b * bw, (b + 1) * bw);
        double d = 10.0 * log10((pb + 1e-12) / (pr + 1e-12));
        err2 += d * d;
    }
    return sqrt(err2 / nbands);
}
#endif /* BLOCKING — only the block-only band error metric; peak/resolve
         metrics below take a raw PSD and are reused by the full path */

/* peak-frequency error for the 1 kHz tone: strongest bin within one bin
 * width of 1 kHz (the interval always contains a bin center) */
static double peak_error_hz(const double* psd, int B, int sr) {
    double df = (double)sr / B;
    int best = -1; double bestp = -1.0;
    for (int k = 0; k <= B / 2; k++) {
        double f = k * df;
        if (fabs(f - 1000.0) > df) continue;
        if (psd[k] > bestp) { bestp = psd[k]; best = k; }
    }
    return fabs((best >= 0 ? best * df : 0.0) - 1000.0);
}

/* can the 2.0/2.35 kHz tone pair be resolved? two local maxima with a
 * 3 dB+ valley between them in the 1.7-2.7 kHz window */
static int pair_resolved(const double* psd, int B, int sr) {
    double df = (double)sr / B;
    int k0 = (int)(1700.0 / df), k1 = (int)(2700.0 / df);
    if (k1 > B / 2) k1 = B / 2;
    int peaks = 0; double peakval[8]; int npk = 0;
    for (int k = k0 + 1; k < k1; k++) {
        if (psd[k] > psd[k - 1] && psd[k] >= psd[k + 1]) {
            peaks++;
            if (npk < 8) peakval[npk++] = psd[k];
        }
    }
    if (peaks < 2) return 0;
    double lo = psd[k0];
    for (int k = k0; k <= k1; k++) if (psd[k] < lo) lo = psd[k];
    double hi = peakval[0];
    for (int i = 1; i < npk; i++) if (peakval[i] > hi) hi = peakval[i];
    return 10.0 * log10(hi / (lo + 1e-20)) >= 3.0;
}

/* time a single radix-2 FFT of size n, repeated enough to be measurable;
 * returns microseconds per transform */
static double time_one_fft(int n, int reps) {
    C* z = (C*)malloc((size_t)n * sizeof(C));
    for (int i = 0; i < n; i++) { z[i].re = (float)(i % 13) * 0.1f; z[i].im = 0.0f; }
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int r = 0; r < reps; r++) fft(z, n, 0);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    free(z);
    double us = (t1.tv_sec - t0.tv_sec) * 1e6 + (t1.tv_nsec - t0.tv_nsec) / 1e3;
    return us / reps;
}

int main(void) {
    const int rates[2] = { 48000, 96000 };
    const int N = 4096;                       /* fixed 0.5 s signal at 8 kHz-ish... */
    const int blocks[] = {64, 128, 256, 512, 1024, 2048};
    const int nblk = (int)(sizeof(blocks) / sizeof(blocks[0]));

    float* x = (float*)malloc((size_t)N * sizeof(float));
    double* ref = (double*)malloc((size_t)(N / 2 + 1) * sizeof(double));

    printf("rate_hz,block,mode,navg,latency_ms,band_err_db,peak_err_hz,resolved,fft_us,total_ms\n");

    for (int ri = 0; ri < 2; ri++) {
        int sr = rates[ri];
        gen_signal(x, N, sr);
        full_spectrum(x, N, (float)sr, ref);  /* blocking-disabled reference */

        /* blocking disabled: one full-segment transform, whole-signal latency */
        double full_us = time_one_fft(N, 200);
        double full_peak = peak_error_hz(ref, N, sr);
        printf("%d,%d,full,1,%.3f,0.000,%.3f,%d,%.3f,%.3f\n",
               sr, N, 1000.0 * N / sr, full_peak, pair_resolved(ref, N, sr),
               full_us, full_us / 1000.0);

#if BLOCKING
        for (int bi = 0; bi < nblk; bi++) {
            int B = blocks[bi];
            BlkFFT st;
            if (blk_init(&st, B, (float)sr) != 0) { fprintf(stderr, "bad block %d\n", B); return 1; }

            struct timespec t0, t1;
            clock_gettime(CLOCK_MONOTONIC, &t0);
            blk_push(&st, x, N);
            blk_flush(&st);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double total_ms = (t1.tv_sec - t0.tv_sec) * 1e3 + (t1.tv_nsec - t0.tv_nsec) / 1e6;

            double lsd = band_error_db(st.psd, B, sr, ref, N);
            double peak = peak_error_hz(st.psd, B, sr);
            int resolved = pair_resolved(st.psd, B, sr);
            double lat_ms = 1000.0 * B / sr;                 /* one-block algorithmic latency */
            double fft_us = time_one_fft(B, 2000);
            printf("%d,%d,block,%d,%.3f,%.3f,%.3f,%d,%.3f,%.3f\n",
                   sr, B, st.navg, lat_ms, lsd, peak, resolved, fft_us, total_ms);
            blk_free(&st);
        }
#else
        /* blocking disabled: every nominal block length routes through the
         * single full-segment FFT, i.e. the pre-blocking per-sample path */
        (void)nblk;
        for (int bi = 0; bi < nblk; bi++) {
            int B = blocks[bi];
            struct timespec t0, t1;
            clock_gettime(CLOCK_MONOTONIC, &t0);
            full_spectrum(x, N, (float)sr, ref);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double total_ms = (t1.tv_sec - t0.tv_sec) * 1e3 + (t1.tv_nsec - t0.tv_nsec) / 1e6;
            printf("%d,%d,block,0,%.3f,0.000,0.000,0,%.3f,%.3f\n",
                   sr, B, 1000.0 * B / sr, time_one_fft(B, 2000), total_ms);
        }
#endif
    }
    free(x); free(ref);
    return 0;
}
