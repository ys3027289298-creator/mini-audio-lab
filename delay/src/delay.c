
/* delay: cross-correlation delay estimation with a Cramer-Rao lower bound.

   Model.  Reference frame x[n] (unit mean square) and delayed frame y[n]:

       y[n] = a*x(n - tau) + w[n],   n = 0..N-1,  tau in samples (fractional)

   a is a gain and w is white noise with variance sigma^2.  tau is found by
   locating the peak of the circular cross-correlation

       r[k] = sum_n x[n]*y[n+k],

   then refining it to a fraction of a sample by interpolating r at arbitrary
   lag from the cross-spectrum.

   Accuracy floor.  For a periodic, band-limited x the Fisher information on
   tau is

       J(tau) = a^2 * sum_n x'(n)^2 / sigma^2 = SNR * sum_n x'(n)^2,

   and by Parseval sum_n x'(n)^2 = sum_k (2*pi*k/N)^2 |X[k]|^2.  Hence

       CRLB(tau) = 1 / (SNR * N * <W^2>),          samples^2
       <W^2>     = sum_k (2*pi*k/N)^2 |X[k]|^2 / sum_k |X[k]|^2,
       SNR       = a^2 * mean(x^2) / sigma^2.

   This is a lower bound on any unbiased estimator's variance; efficiency is
   CRLB / measured variance (<= 1 in expectation).

   Detectability.  With no delayed copy the normalised peak

       z = sqrt(N) * max_tau |r(tau)| / sqrt(E_x * E_y)

   is distributed like the maximum of ~2*B*N unit-normal lags: it stays near
   sqrt(2*ln N) (about 3.5 here) and falls off smoothly.  By Cauchy-Schwarz
   z <= sqrt(N) always.  z below Z_MIN means the peak is indistinguishable
   from noise, so the delay is reported unmeasurable instead of as a number. */

#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

#define N 512          /* frame length (power of two) */
#define BAND 0.35      /* occupied band edge, cycles/sample */
#define FINE 32        /* fine-search steps per sample around the peak */

#define Z_MIN 6.0      /* detection gate on the normalised peak */

typedef struct {
    double tau;        /* fractional delay of y relative to x, samples */
    double z;          /* normalised correlation peak (detection statistic) */
    double crlb;       /* CRLB variance from this frame, samples^2 (NaN if not measured) */
    int    measured;   /* 1 when z >= Z_MIN */
} DelayEstimate;

/* ---- minimal radix-2 FFT on doubles ------------------------------------- */

typedef struct { double re, im; } C;

static void fft(C* a, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { C t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * PI / len * (inverse ? 1.0 : -1.0);
        C wl = { cos(ang), sin(ang) };
        for (int i = 0; i < n; i += len) {
            C w = { 1.0, 0.0 };
            for (int k = 0; k < len / 2; k++) {
                C u = a[i + k];
                C v = { a[i + k + len / 2].re * w.re - a[i + k + len / 2].im * w.im,
                        a[i + k + len / 2].re * w.im + a[i + k + len / 2].im * w.re };
                a[i + k].re = u.re + v.re; a[i + k].im = u.im + v.im;
                a[i + k + len / 2].re = u.re - v.re; a[i + k + len / 2].im = u.im - v.im;
                double nr = w.re * wl.re - w.im * wl.im;
                double ni = w.re * wl.im + w.im * wl.re;
                w.re = nr; w.im = ni;
            }
        }
    }
    if (inverse) { for (int i = 0; i < n; i++) { a[i].re /= n; a[i].im /= n; } }
}

/* ---- correlation and its interpolation ---------------------------------- */

/* Cross-spectrum conj(X)*Y, and the real energies of both frames. */
static void cross_spectrum(const float* x, const float* y, C* X, C* Y, C* Cxy,
                           double* ex, double* ey) {
    for (int n = 0; n < N; n++) {
        X[n] = (C){ x[n], 0.0 };
        Y[n] = (C){ y[n], 0.0 };
    }
    fft(X, N, 0);
    fft(Y, N, 0);
    for (int k = 0; k < N; k++)
        Cxy[k] = (C){ X[k].re * Y[k].re + X[k].im * Y[k].im,
                      X[k].re * Y[k].im - X[k].im * Y[k].re };
    double sx = 0.0, sy = 0.0;
    for (int n = 0; n < N; n++) { sx += (double)x[n] * x[n]; sy += (double)y[n] * y[n]; }
    *ex = sx; *ey = sy;
}

/* Signed DFT bin index: k in [0,N/2] is k, k in (N/2,N) is k-N.  Band-limited
   interpolation and phase ramps must use this so that negative frequencies
   stay conjugate-symmetric. */
static int k_signed(int k) { return (k <= N / 2) ? k : k - N; }

/* Interpolated correlation at fractional lag tau:
   r(tau) = (1/N) sum_k Cxy[k] exp(+2*pi*i*k_signed(k)*tau/N), real part. */
static double r_at(const C* Cxy, double tau) {
    double s = 0.0;
    for (int k = 0; k < N; k++) {
        double ang = 2.0 * PI * k_signed(k) * tau / (double)N;
        s += Cxy[k].re * cos(ang) - Cxy[k].im * sin(ang);
    }
    return s / (double)N;
}

/* Mean-square angular frequency <W^2> of a real frame's spectrum. */
static double mean_sq_freq(const C* X) {
    double num = 0.0, den = 0.0;
    for (int k = 0; k < N; k++) {
        double wk = 2.0 * PI * k_signed(k) / (double)N;
        double p = X[k].re * X[k].re + X[k].im * X[k].im;
        num += wk * wk * p;
        den += p;
    }
    return den > 0.0 ? num / den : 0.0;
}

/* ---- estimator ---------------------------------------------------------- */

/* Estimate tau of y relative to x. Returns 0 on success, -1 on bad input.
   On success fills est; *coarse_k receives the integer-lag peak.
   Note: the scratch buffers below are static to keep this off the stack, so
   this function is not reentrant; call it from one thread at a time. */
static int estimate_delay(const float* x, const float* y, DelayEstimate* est, int* coarse_k) {
    if (x == NULL || y == NULL || est == NULL) return -1;
    memset(est, 0, sizeof(*est));

    static C X[N], Y[N], Cxy[N], R[N];
    double ex, ey;
    cross_spectrum(x, y, X, Y, Cxy, &ex, &ey);
    if (ex <= 0.0 || ey <= 0.0) return -1;

    /* coarse integer peak via inverse transform of the cross-spectrum */
    memcpy(R, Cxy, sizeof(R));
    fft(R, N, 1);
    int k0 = 0;
    double best = -HUGE_VAL;
    for (int k = 0; k < N; k++) if (R[k].re > best) { best = R[k].re; k0 = k; }
    if (coarse_k) *coarse_k = k0;

    /* fine search +-1 sample around k0, then a parabolic fit on three points */
    double step = 1.0 / FINE;
    double tbest = k0, vbest = -HUGE_VAL;
    for (int i = -FINE; i <= FINE; i++) {
        double t = k0 + i * step;
        double v = r_at(Cxy, t);
        if (v > vbest) { vbest = v; tbest = t; }
    }
    double vm = r_at(Cxy, tbest - step);
    double vp = r_at(Cxy, tbest + step);
    double den = vm - 2.0 * vbest + vp;
    double delta = (den != 0.0) ? 0.5 * (vm - vp) / den * step : 0.0;
    if (delta < -step) delta = -step;
    if (delta >  step) delta =  step;
    double tau = tbest + delta;

    est->tau = tau;

    /* detection statistic from the interpolated peak (z <= sqrt(N)) */
    double rpeak = r_at(Cxy, tau);
    double z = sqrt((double)N) * rpeak / sqrt(ex * ey);
    if (z < 0.0) z = 0.0;
    est->z = z;
    est->measured = (z >= Z_MIN);

    /* data-driven CRLB: SNR_est = z^2/(N - z^2), so
       CRLB = (N - z^2) / (N * z^2 * <W^2>) */
    double w2 = mean_sq_freq(X);
    double z2 = z * z;
    if (w2 > 0.0 && z2 > 0.0 && z2 < (double)N)
        est->crlb = ((double)N - z2) / ((double)N * z2 * w2);
    else
        est->crlb = HUGE_VAL;

    if (!est->measured) est->crlb = NAN;
    return 0;
}

/* Analytic CRLB variance for a known SNR and reference spectrum. */
static double crlb_for_snr(const C* X, double snr) {
    double w2 = mean_sq_freq(X);
    return (snr > 0.0 && w2 > 0.0) ? 1.0 / (snr * (double)N * w2) : HUGE_VAL;
}

/* ---- deterministic test harness ----------------------------------------- */

static unsigned long long rng_state = 0x9E3779B97F4A7C15ULL;
static unsigned rng_u32(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (unsigned)(rng_state >> 32);
}
static double rng_uniform(void) { return (rng_u32() + 0.5) / 4294967296.0; }
static double rng_normal(void) {
    double u1 = rng_uniform(), u2 = rng_uniform();
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* Build a band-limited reference x (unit mean square) and its exact delayed
   copy xd[n] = x(n - tau0), both from the fixed seed. */
static void build_signal(float* x, float* xd, C* X, double tau0) {
    static C Xd[N];
    int K = (int)lround(BAND * N);
    memset(X, 0, sizeof(C) * N);
    for (int k = 1; k <= K; k++) {
        double re = rng_normal(), im = rng_normal();
        X[k] = (C){ re, im };
        X[N - k] = (C){ re, -im };
    }
    memcpy(Xd, X, sizeof(Xd));
    for (int k = 0; k < N; k++) {
        double ang = -2.0 * PI * k_signed(k) * tau0 / (double)N;
        double c = cos(ang), s = sin(ang);
        Xd[k] = (C){ X[k].re * c - X[k].im * s, X[k].re * s + X[k].im * c };
    }
    C t[N];
    memcpy(t, X, sizeof(t));  fft(t, N, 1);
    for (int n = 0; n < N; n++) x[n] = (float)t[n].re;
    memcpy(t, Xd, sizeof(t)); fft(t, N, 1);
    for (int n = 0; n < N; n++) xd[n] = (float)t[n].re;

    /* scale x and xd together to unit mean square, so SNR is in sample power */
    double ms = 0.0;
    for (int n = 0; n < N; n++) ms += (double)x[n] * x[n];
    ms /= N;
    float g = (float)(1.0 / sqrt(ms));
    for (int n = 0; n < N; n++) { x[n] *= g; xd[n] *= g; }
}

int main(void) {
    const double tau0 = 7.37;               /* delay used to build xd */
    float x[N], xd[N], y[N];
    C X[N];
    build_signal(x, xd, X, tau0);

    /* Verify the construction: the noiseless cross-correlation of x with xd
       must peak at tau0.  Take this measured peak as ground truth so a broken
       signal model cannot silently pass the accuracy tests below. */
    double truth;
    {
        static C Y0[N], C0[N];
        double ex, ey;
        cross_spectrum(x, xd, X, Y0, C0, &ex, &ey);
        if (fabs(ex - ey) > 1e-3 * ex)
            printf("# warning: |x|^2 and |xd|^2 differ (%.3f vs %.3f)\n", ex, ey);
        DelayEstimate e0;
        estimate_delay(x, xd, &e0, NULL);
        truth = e0.tau;
        check("construction_peak", fabs(truth - tau0) < 0.01);
        check("construction_z", e0.z > 0.99 * sqrt((double)N));
    }
    const double tau_true = truth;          /* measured delay of xd vs x */

    /* ---- 1. cross-correlation curve at a clean SNR ---- */
    const double snr_curve_db = 20.0;
    double a_curve = pow(10.0, snr_curve_db / 20.0);
    for (int n = 0; n < N; n++)
        y[n] = (float)(a_curve * xd[n] + rng_normal());

    {
        static C Y[N], Cxy[N];
        double ex, ey;
        cross_spectrum(x, y, X, Y, Cxy, &ex, &ey);
        double norm = sqrt(ex * ey);
        printf("# cross-correlation curve at SNR=%.0f dB (lag, normalised r)\n", snr_curve_db);
        for (int i = -30; i <= 30; i++) {
            double tt = tau_true + i * 0.1;
            printf("xcorr,%.2f,%.6f\n", tt, r_at(Cxy, tt) / norm);
        }
        printf("# integer-lag samples (what a coarse argmax sees)\n");
        for (int k = 4; k <= 11; k++)
            printf("coarse,%d,%.6f\n", k, r_at(Cxy, k) / norm);
        printf("# end curve\n");
    }

    /* ---- 2. fractional vs integer estimate at a good SNR ---- */
    {
        DelayEstimate e;
        int k0;
        check("estimate_ok", estimate_delay(x, y, &e, &k0) == 0);
        printf("estimate,%d,%.5f,%.5f\n", k0, e.tau, tau_true);
        printf("z_peak=%.2f\n", e.z);
        check("fractional_near_true", fabs(e.tau - tau_true) < 0.05);
        check("beats_integer_rounding", fabs(e.tau - tau_true) < fabs((double)k0 - tau_true));
    }

    /* ---- 3. Monte-Carlo efficiency at two SNRs ---- */
    const double snr_lo = 0.0, snr_hi = 10.0;     /* two evaluation points */
    int T = 400;
    double crlb_hi = crlb_for_snr(X, pow(10.0, snr_hi / 10.0));
    double crlb_lo = crlb_for_snr(X, pow(10.0, snr_lo / 10.0));
    double mse_hi = 0.0, mse_lo = 0.0, mse_int = 0.0;

    rng_state = 0x243F6A8885A308D3ULL;             /* fresh noise stream */
    for (int t = 0; t < T; t++) {
        double a = pow(10.0, snr_hi / 20.0);
        for (int n = 0; n < N; n++) y[n] = (float)(a * xd[n] + rng_normal());
        DelayEstimate e; int k0;
        estimate_delay(x, y, &e, &k0);
        if (e.measured) {
            double d = e.tau - tau_true;                  /* unwrap into [-N/2, N/2) */
            d = fmod(d + N / 2.0, (double)N);
            if (d < 0) d += (double)N;
            d -= N / 2.0;
            mse_hi += d * d;
        }
        double eint = (double)k0 - tau_true;              /* integer-only error */
        mse_int += eint * eint;
    }
    for (int t = 0; t < T; t++) {
        double a = pow(10.0, snr_lo / 20.0);
        for (int n = 0; n < N; n++) y[n] = (float)(a * xd[n] + rng_normal());
        DelayEstimate e; int k0;
        estimate_delay(x, y, &e, &k0);
        if (e.measured) {
            double d = e.tau - tau_true;
            d = fmod(d + N / 2.0, (double)N);
            if (d < 0) d += (double)N;
            d -= N / 2.0;
            mse_lo += d * d;
        }
    }
    mse_hi /= T; mse_lo /= T; mse_int /= T;

    double rms_hi = sqrt(mse_hi), rms_lo = sqrt(mse_lo);
    double eff_hi = crlb_hi / mse_hi, eff_lo = crlb_lo / mse_lo;

    printf("crlb_sd_10db=%.5f\n", sqrt(crlb_hi));
    printf("rmse_10db=%.5f\n", rms_hi);
    printf("efficiency_10db=%.3f\n", eff_hi);
    printf("crlb_sd_0db=%.5f\n", sqrt(crlb_lo));
    printf("rmse_0db=%.5f\n", rms_lo);
    printf("efficiency_0db=%.3f\n", eff_lo);
    printf("integer_only_rmse=%.4f\n", sqrt(mse_int));

    /* ---- 3b. efficiency across a range of SNRs (for the figure) ---- */
    {
        static const double sweep[] = { 20, 15, 10, 5, 0, -5, -10, -12, -15, -20, -25 };
        int ns = (int)(sizeof(sweep) / sizeof(sweep[0])), Ts = 1500;
        printf("# snr sweep: snr_db,crlb_sd,rmse,eff,measured_frac,gross_frac\n");
        for (int s = 0; s < ns; s++) {
            double snr_lin = pow(10.0, sweep[s] / 10.0);
            double a = pow(10.0, sweep[s] / 20.0);
            double cl = crlb_for_snr(X, snr_lin), mse = 0.0; int det = 0, gross = 0;
            for (int t = 0; t < Ts; t++) {
                for (int n = 0; n < N; n++) y[n] = (float)(a * xd[n] + rng_normal());
                DelayEstimate e; estimate_delay(x, y, &e, NULL);
                if (!e.measured) continue;
                det++;
                double d = e.tau - tau_true;
                d = fmod(d + N / 2.0, (double)N);
                if (d < 0) d += (double)N;
                d -= N / 2.0;
                mse += d * d;
                if (fabs(d) >= 0.5) gross++;   /* wrong correlation lobe -> > half a sample off */
            }
            double eff = det > 0 ? (cl / (mse / det)) : NAN;
            printf("sweep,%.0f,%.5f,%.5f,%.3f,%.3f,%.4f\n", sweep[s], sqrt(cl),
                   det > 0 ? sqrt(mse / det) : NAN, eff, (double)det / Ts,
                   (double)gross / Ts);
        }
        printf("# end sweep\n");
    }

    check("efficiency_10db_near_one", eff_hi > 0.5 && eff_hi < 1.5);
    check("efficiency_0db_near_one", eff_lo > 0.5 && eff_lo < 1.5);
    check("bound_holds_10db", mse_hi >= 0.8 * crlb_hi);
    check("bound_holds_0db", mse_lo >= 0.8 * crlb_lo);
    check("beats_integer_by_10x", rms_lo < 0.1 * sqrt(mse_int));

    /* CRLB must fall as 1/SNR: 10 dB more SNR -> variance down by 10x */
    double ratio = crlb_lo / crlb_hi;
    printf("crlb_ratio=%.2f\n", ratio);
    check("crlb_scales_as_1_over_snr", fabs(ratio - 10.0) < 0.1);

    /* ---- 4. too noisy: must report unmeasurable ---- */
    {
        const double snr_dead_db = -30.0;
        double a = pow(10.0, snr_dead_db / 20.0);
        int not_measured = 0;
        double zsum = 0.0;
        DelayEstimate last = {0};
        for (int t = 0; t < 50; t++) {
            for (int n = 0; n < N; n++) y[n] = (float)(a * xd[n] + rng_normal());
            DelayEstimate e;
            estimate_delay(x, y, &e, NULL);
            if (!e.measured) not_measured++;
            zsum += e.z;
            last = e;
        }
        printf("z_mean_low_snr=%.2f\n", zsum / 50.0);
        printf("unmeasurable_frac=%.2f\n", not_measured / 50.0);
        check("low_snr_unmeasurable", not_measured == 50);
        check("low_snr_below_gate", zsum / 50.0 < Z_MIN);
        check("low_snr_crlb_nan", isnan(last.crlb));
    }

    /* ---- 4b. the gate separates the two regimes: at +10 dB almost every
       trial is measurable, at -20 dB almost none is ---- */
    {
        int ok_hi = 0, ok_lo = 0, T2 = 40;
        for (int t = 0; t < T2; t++) {
            double a = pow(10.0, 10.0 / 20.0);
            for (int n = 0; n < N; n++) y[n] = (float)(a * xd[n] + rng_normal());
            DelayEstimate e; estimate_delay(x, y, &e, NULL);
            ok_hi += e.measured;
        }
        for (int t = 0; t < T2; t++) {
            double a = pow(10.0, -20.0 / 20.0);
            for (int n = 0; n < N; n++) y[n] = (float)(a * xd[n] + rng_normal());
            DelayEstimate e; estimate_delay(x, y, &e, NULL);
            ok_lo += e.measured;
        }
        printf("measured_frac_10db=%.2f\n", (double)ok_hi / T2);
        printf("measured_frac_m20db=%.2f\n", (double)ok_lo / T2);
        check("gate_high_snr_measurable", ok_hi == T2);
        check("gate_low_snr_rejected", ok_lo == 0);
    }

    /* ---- 5. input validation ---- */
    {
        DelayEstimate e;
        float zeros[N];
        memset(zeros, 0, sizeof(zeros));
        check("null_input", estimate_delay(NULL, y, &e, NULL) == -1);
        check("silent_frame", estimate_delay(zeros, y, &e, NULL) == -1);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
