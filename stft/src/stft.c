
/* stft: STFT framing timestamp alignment — analytic deviation budget vs measured.
 *
 * Frame m covers samples [mH, mH+N) with a length-N Hann window.
 * Reported timestamps: frame start mH/fs, frame center (mH+(N-1)/2)/fs.
 *
 * Timestamp deviation budget (two independent worst-case terms):
 *   grid quantization  : H/(2*fs)   event phase against the hop grid
 *   window smearing    : N/(2*fs)   event anywhere inside the window half-width
 *   delta_budget       = (H + N) / (2*fs)
 *
 * Oversampling factor rho = N/H gives delta_budget = N*(1+1/rho)/(2*fs).
 *
 * Two assignment conventions are measured on a dense sweep of fractional event
 * positions (real onsets need not land on samples):
 *   A nearest-center : frame whose center is closest to the event
 *   B onset-frame    : first frame whose window covers the event, center stamp
 * An actual Hann-windowed energy framing on an interpolated impulse checks that
 * the max-energy frame always equals the nearest-center frame.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define FS 48000.0
#define DELTA_MAX_MS 5.0          /* pass/fail timestamp budget */
#define SWEEP_CELLS 64           /* hop cells swept per config */
#define PTS_PER_CELL 128         /* event positions per hop cell */

typedef struct { int N, H; int expect_pass; const char* tag; } Cfg;

static const Cfg CFGS[] = {
    {256, 128, 1, "chosen  N=256 H=128 rho=2"},
    {256,  64, 1, "dense   N=256 H=64  rho=4"},
    {128,  64, 1, "small   N=128 H=64  rho=2"},
    {256, 256, 0, "critical N=256 H=256 rho=1"},
    {512, 256, 0, "wide    N=512 H=256 rho=2"},
    {512, 512, 0, "critbig N=512 H=512 rho=1"},
};
#define NCFG (int)(sizeof(CFGS) / sizeof(CFGS[0]))

static double hann(int n, int N) {
    return 0.5 * (1.0 - cos(2.0 * PI * n / (N - 1)));
}

/* frame index whose center is nearest to fractional event position p */
static int nearest_frame(double p, int N, int H) {
    double c0 = (N - 1) / 2.0;
    int m = (int)lround((p - c0) / H);
    return m;
}

/* first frame whose window [mH, mH+N) contains p */
static int onset_frame(double p, int N, int H) {
    return (int)floor((p - N) / H) + 1;
}

static double center_ms(int m, int N, int H) {
    return (m * H + (N - 1) / 2.0) / FS * 1000.0;
}

int main(int argc, char** argv) {
    int failures = 0;
    FILE* f_sum = NULL;
    FILE* f_phase = NULL;
    FILE* f_ts_pass = NULL;
    FILE* f_ts_fail = NULL;
    if (argc > 1) {
        char path[512];
        snprintf(path, sizeof(path), "%s/summary.csv", argv[1]);
        f_sum = fopen(path, "w");
        snprintf(path, sizeof(path), "%s/phase_sweep.csv", argv[1]);
        f_phase = fopen(path, "w");
        snprintf(path, sizeof(path), "%s/timestamps_pass.csv", argv[1]);
        f_ts_pass = fopen(path, "w");
        snprintf(path, sizeof(path), "%s/timestamps_fail.csv", argv[1]);
        f_ts_fail = fopen(path, "w");
        if (!f_sum || !f_phase || !f_ts_pass || !f_ts_fail) {
            printf("csv_open=0\n");
            return 1;
        }
        fprintf(f_sum, "tag,N,H,rho,grid_ms,smear_ms,budget_ms,meas_nearest_ms,"
                       "meas_onset_ms,delta_max_ms,verdict\n");
        fprintf(f_phase, "phase_frac,nearest_ms,onset_ms,event_ms\n");
        fprintf(f_ts_pass, "event,true_ms,nearest_frame,nearest_center_ms,"
                           "nearest_dev_ms,onset_frame,onset_center_ms,"
                           "onset_dev_ms,within_budget\n");
        fprintf(f_ts_fail, "event,true_ms,nearest_frame,nearest_center_ms,"
                           "nearest_dev_ms,onset_frame,onset_center_ms,"
                           "onset_dev_ms,within_budget\n");
    }

    for (int ci = 0; ci < NCFG; ci++) {
        int N = CFGS[ci].N, H = CFGS[ci].H;
        double rho = (double)N / H;
        int npts = SWEEP_CELLS * PTS_PER_CELL;
        double p0 = 2.0 * N;                 /* keep all used frames interior */
        double step = (double)H / PTS_PER_CELL;
        int nframes = SWEEP_CELLS + 2 * (N / H) + 4;
        int buf_len = nframes * H + N + 4;

        double* x = (double*)calloc((size_t)buf_len, sizeof(double));
        double max_near = 0.0, max_onset = 0.0;
        int winner_mismatch = 0;

        for (int k = 0; k < npts; k++) {
            double p = p0 + k * step;      /* fractional event position */
            memset(x, 0, (size_t)buf_len * sizeof(double));
            int i0 = (int)floor(p);
            double frac = p - i0;
            x[i0] += 1.0 - frac;           /* unit-mass linear-interp impulse */
            x[i0 + 1] += frac;

            int mA = nearest_frame(p, N, H);
            int mB = onset_frame(p, N, H);
            double dA = p / FS * 1000.0 - center_ms(mA, N, H);
            double dB = p / FS * 1000.0 - center_ms(mB, N, H);
            if (fabs(dA) > max_near) max_near = fabs(dA);
            if (fabs(dB) > max_onset) max_onset = fabs(dB);

            /* real Hann-windowed frame energies: winner must be nearest frame.
             * At exact half-hop phases two frames tie analytically (identical
             * Hann energy); a tolerance absorbs FP summation-order rounding so
             * the tie is resolved in favor of mA. */
            double eA = 0.0;
            for (int n = 0; n < N; n++) {
                double s = x[mA * H + n];
                double w = hann(n, N);
                eA += w * w * s * s;
            }
            double e_tol = 1e-12 * (eA + 1.0);
            for (int m = mA - 2; m <= mA + 2; m++) {
                if (m == mA) continue;
                double E = 0.0;
                for (int n = 0; n < N; n++) {
                    double s = x[m * H + n];
                    double w = hann(n, N);
                    E += w * w * s * s;
                }
                if (E > eA + e_tol) winner_mismatch++;
            }

            /* phase sweep of the chosen config, one hop cell */
            if (f_phase && N == 256 && H == 128 && k < PTS_PER_CELL) {
                fprintf(f_phase, "%.6f,%.6f,%.6f,%.6f\n",
                        (p - p0) / H, dA, dB, p / FS * 1000.0);
            }
        }

        double grid_ms = H / (2.0 * FS) * 1000.0;
        double smear_ms = N / (2.0 * FS) * 1000.0;
        double budget_ms = grid_ms + smear_ms;
        int pass = budget_ms <= DELTA_MAX_MS;
        double tol_ms = 0.05 / FS * 1000.0;   /* sweep resolution margin */

        /* self-tests: measured envelopes sit inside the analytic terms and
         * below the budget; sweep is fractional so maxima approach but do not
         * overshoot the theoretical H/2 and N/2 boundaries */
        int ok_cover = (H <= N);
        int ok_winner = (winner_mismatch == 0);
        int ok_near_env = max_near <= grid_ms + tol_ms;
        int ok_onset_env = max_onset <= smear_ms + tol_ms;
        int ok_under_budget = max_near <= budget_ms + tol_ms &&
                              max_onset <= budget_ms + tol_ms;
        int ok_verdict = (pass == CFGS[ci].expect_pass);
        if (!ok_cover || !ok_winner || !ok_near_env || !ok_onset_env ||
            !ok_under_budget || !ok_verdict)
            failures++;

        printf("%s: cover=%d winner=%d near_env=%d onset_env=%d under_budget=%d "
               "verdict=%d (meas near=%.4fms onset=%.4fms budget=%.4fms)\n",
               CFGS[ci].tag, ok_cover, ok_winner, ok_near_env, ok_onset_env,
               ok_under_budget, ok_verdict, max_near, max_onset, budget_ms);

        if (f_sum) {
            fprintf(f_sum, "%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%s\n",
                    CFGS[ci].tag, N, H, rho, grid_ms, smear_ms, budget_ms,
                    max_near, max_onset, (double)DELTA_MAX_MS,
                    pass ? "PASS" : "FAIL");
        }
        free(x);
    }

    /* named timestamp table: adversarial phases for chosen and critical config */
    if (f_ts_pass && f_ts_fail) {
        const int pairs[][2] = {{256, 128}, {256, 256}};
        for (int pi = 0; pi < 2; pi++) {
            int N = pairs[pi][0], H = pairs[pi][1];
            FILE* f = pi == 0 ? f_ts_pass : f_ts_fail;
            double budget_ms = (H + N) / (2.0 * FS) * 1000.0;
            /* fractional offsets spanning one hop cell, incl. both boundaries */
            double off[9] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.999};
            double base = 4.0 * N + 0.37 * H;
            for (int k = 0; k < 9; k++) {
                double p = base + off[k] * H;
                int mA = nearest_frame(p, N, H);
                int mB = onset_frame(p, N, H);
                double dA = p / FS * 1000.0 - center_ms(mA, N, H);
                double dB = p / FS * 1000.0 - center_ms(mB, N, H);
                fprintf(f, "e%d,%.6f,%d,%.6f,%.6f,%d,%.6f,%.6f,%d\n", k,
                        p / FS * 1000.0, mA, center_ms(mA, N, H), dA, mB,
                        center_ms(mB, N, H), dB,
                        fabs(dA) <= budget_ms && fabs(dB) <= budget_ms);
            }
        }
    }

    if (f_sum) fclose(f_sum);
    if (f_phase) fclose(f_phase);
    if (f_ts_pass) fclose(f_ts_pass);
    if (f_ts_fail) fclose(f_ts_fail);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
