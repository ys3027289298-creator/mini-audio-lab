
/* aliasing_bound: aliasing residual of decimation vs FIR order.
 *
 * Model: input rate fs, output rate r*fs with r = 1/D (integer decimation).
 * Higher r -> more output samples (slower) but wider filter transition band
 * (less aliasing). Anti-alias filter: Kaiser-windowed sinc, order ORDER,
 * passband edge wp = PB*pi/D (preserve PB of the new Nyquist band),
 * stopband edge ws = pi/D (fold frequency), transition dw = (1-PB)*pi/D.
 *
 * Bound: content aliasing into the preserved band [0, D*wp] comes only from
 * input bands [2*pi*k/D - wp, 2*pi*k/D + wp] (k >= 1), all inside the
 * stopband since 2*pi/D - wp > pi/D. With |H| <= ds in the stopband and
 * Kaiser's order formula N = (A-8)/(2.285*dw), A = -20*log10(ds):
 *
 *     A(N,D)  = 2.285*N*dw + 8 = 2.285*N*(1-PB)*pi/D + 8
 *     eps(D) <= 10^(-A/20)                          (analytic upper bound)
 *
 * Inverted for target residual eps_t (A_req = -20*log10(eps_t)):
 *
 *     r_min = (A_req - 8) / (2.285*N*(1-PB)*pi)
 *
 * Realizable gears are r = 1/D, so the minimum-ratio gear is D0 = floor(1/r_min).
 * This program measures aliasing at D0-1, D0, D0+1 and checks the results.
 */
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

#define ORDER 128          /* FIR order N (length N+1) */
#define PB 0.8             /* preserved fraction of the output Nyquist band */
#define EPS_TARGET 1e-3    /* target aliasing residual (amplitude ratio) */
#define LIN (1 << 18)      /* input length for the tone measurements */
#define NGRID 20001        /* frequency-scan grid points */

static double bessel_i0(double x) {
    double term = 1.0, sum = 1.0, y = (x * 0.5) * (x * 0.5);
    for (int k = 1; k < 1000; k++) {
        term *= y / ((double)k * k);
        sum += term;
        if (term < 1e-13 * sum) break;
    }
    return sum;
}

static double kaiser_beta(double A) {
    if (A > 50.0) return 0.1102 * (A - 8.7);
    if (A >= 21.0) return 0.5842 * pow(A - 21.0, 0.4) + 0.07886 * (A - 21.0);
    return 0.0;
}

/* Design Kaiser-windowed sinc lowpass for decimation by D; returns the
 * attenuation predicted by the inverted Kaiser order formula. */
static double design(int D, int order, double* h) {
    double wp = PB * PI / D, ws = PI / D;
    double wc = 0.5 * (wp + ws);
    double A = 2.285 * order * (ws - wp) + 8.0;
    double beta = kaiser_beta(A);
    double i0b = bessel_i0(beta);
    double sum = 0.0;
    for (int n = 0; n <= order; n++) {
        double m = n - order * 0.5;
        double s = (m == 0.0) ? 1.0 : sin(wc * m) / (wc * m);
        double u = 2.0 * n / order - 1.0;
        double w = bessel_i0(beta * sqrt(1.0 - u * u)) / i0b;
        h[n] = s * (wc / PI) * w;
        sum += h[n];
    }
    for (int n = 0; n <= order; n++) h[n] /= sum;
    return A;
}

static double h_mag(const double* h, int len, double w) {
    double re = 0.0, im = 0.0;
    for (int n = 0; n < len; n++) { re += h[n] * cos(w * n); im -= h[n] * sin(w * n); }
    return hypot(re, im);
}

static double fold_freq(double w, int D) {
    double Om = fmod(w * D, 2.0 * PI);
    if (Om < 0.0) Om += 2.0 * PI;
    if (Om > PI) Om = 2.0 * PI - Om;
    return Om;
}

/* Scan the stopband: peak |H| anywhere in [ws, pi] (conservative worst case)
 * and peak |H| over the subset that folds into the preserved band. */
static void scan_stopband(const double* h, int len, int D,
                          double* full_peak, double* w_full,
                          double* alias_peak, double* w_alias) {
    double ws = PI / D;
    *full_peak = *alias_peak = 0.0;
    *w_full = *w_alias = ws;
    for (int i = 0; i < NGRID; i++) {
        double w = ws + (PI - ws) * i / (NGRID - 1);
        double m = h_mag(h, len, w);
        if (m > *full_peak) { *full_peak = m; *w_full = w; }
        if (fold_freq(w, D) <= PB * PI + 1e-9 && m > *alias_peak) {
            *alias_peak = m; *w_alias = w;
        }
    }
}

/* Measure aliasing: filter a tone at w_tone, decimate by D, measure the
 * amplitude at the folded frequency over the steady-state section.
 * cplx=1 uses a complex exponential x[n]=e^(jwn): a real tone folding exactly
 * onto output frequency 0 or pi is sampling-degenerate (cos(pi*m+phi) =
 * (-1)^m*cos(phi) loses the sine quadrature), and its conjugate image leaks
 * into nearby DFT bins; the complex tone measures |H(w)| directly in all
 * cases. cplx=0 (real cosine) is kept as a sanity check at non-degenerate
 * frequencies, where both must agree. */
static double x_re[LIN], x_im[LIN];
static double y_re[LIN / 2 + 2 * ORDER + 64];
static double y_im[LIN / 2 + 2 * ORDER + 64];

static double tone_residual(const double* h, int order, int D, double w_tone, int cplx) {
    int len = order + 1;
    for (int n = 0; n < LIN; n++) {
        x_re[n] = cos(w_tone * n);
        x_im[n] = cplx ? sin(w_tone * n) : 0.0;
    }
    int m_out = (LIN - order) / D;
    for (int m = 0; m < m_out; m++) {
        double are = 0.0, aim = 0.0;
        int base = m * D;
        for (int k = 0; k < len; k++) { are += h[k] * x_re[base + k]; aim += h[k] * x_im[base + k]; }
        y_re[m] = are; y_im[m] = aim;
    }
    int g = 3 * order / D + 8;
    /* complex tone keeps its signed frequency; real tone folds to [0, pi] */
    double Om = cplx ? fmod(w_tone * D, 2.0 * PI) : fold_freq(w_tone, D);
    if (Om < 0.0) Om += 2.0 * PI;
    double re = 0.0, im = 0.0;
    int cnt = 0;
    for (int m = g; m < m_out - g; m++) {
        double c = cos(Om * m), s = -sin(Om * m);
        re += y_re[m] * c - y_im[m] * s;
        im += y_re[m] * s + y_im[m] * c;
        cnt++;
    }
    return (cplx ? 1.0 : 2.0) * hypot(re, im) / cnt;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    static double h[ORDER + 1], hm[ORDER * 2 + 2];

    double A_req = -20.0 * log10(EPS_TARGET);
    double r_min = (A_req - 8.0) / (2.285 * ORDER * (1.0 - PB) * PI);
    int D0 = (int)floor(1.0 / r_min);
    printf("A_req=%.1fdB r_min=%.4f D_gear=%d (r=1/%d=%.4f)\n",
           A_req, r_min, D0, D0, 1.0 / D0);
    printf("D      r     A_pred_dB      bound   meas_alias  meas_sbpeak  A_ach_dB  dev_dB\n");

    double bound[3] = {0}, dev[3] = {0};
    int meas_ok = 1, chain_ok = 1, dev_ok = 1, margin_ok = 1;
    for (int gi = 0; gi < 3; gi++) {
        int D = D0 - 1 + gi;
        double A_pred = design(D, ORDER, h);
        bound[gi] = pow(10.0, -A_pred / 20.0);
        double fp, wf, ap, wa;
        scan_stopband(h, ORDER + 1, D, &fp, &wf, &ap, &wa);
        double meas_alias = tone_residual(h, ORDER, D, wa, 1);
        double meas_full = tone_residual(h, ORDER, D, wf, 1);
        double meas_alias_real = tone_residual(h, ORDER, D, wa, 0);
        double A_ach = -20.0 * log10(fp);
        dev[gi] = A_pred - A_ach;
        printf("%d  %.4f   %7.1f  %.3e   %.3e    %.3e   %7.1f  %+.2f\n",
               D, 1.0 / D, A_pred, bound[gi], meas_alias, meas_full, A_ach, dev[gi]);
        if (meas_alias > bound[gi]) meas_ok = 0;
        if (fabs(meas_alias - ap) > 1e-3 * ap) chain_ok = 0;   /* cplx tone == |H| */
        if (fabs(meas_full - fp) > 1e-3 * fp) chain_ok = 0;
        if (fabs(meas_alias_real - ap) > 1e-3 * ap) chain_ok = 0; /* real == cplx */
        if (fabs(dev[gi]) > 1.2) dev_ok = 0;                   /* Kaiser fit accuracy */
        /* +12% order margin must make achieved attenuation >= prediction */
        int order2 = (int)ceil(ORDER * 1.12);
        design(D, order2, hm);
        double f2, w2, a2, w3;
        scan_stopband(hm, order2 + 1, D, &f2, &w2, &a2, &w3);
        if (-20.0 * log10(f2) < A_pred - 0.05) margin_ok = 0;
    }

    check("bound_shrinks_with_ratio", bound[0] < bound[1] && bound[1] < bound[2]);
    check("inversion_minimal", bound[1] <= EPS_TARGET && bound[2] > EPS_TARGET);
    check("meas_alias_within_bound", meas_ok);
    check("meas_matches_response", chain_ok);
    check("deviation_within_kaiser_fit", dev_ok);
    check("order_margin_corrects", margin_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
