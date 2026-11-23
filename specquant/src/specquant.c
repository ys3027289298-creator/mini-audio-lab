
/* specquant: pick a spectral-quantization bit-width from the
 * spectral-entropy information lower bound.
 *
 * Method:
 *   1. one-sided power spectrum -> normalized p_k -> spectral entropy
 *      H = -sum p_k log2 p_k (bit);  K_eff = 2^H effective bins.
 *   2. lower bound L = H bit/coefficient: the entropy of the spectral
 *      distribution floors the average code length per coefficient
 *      (Shannon source coding). A Gaussian entropy-rate estimate is
 *      printed alongside as a stricter reference.
 *   3. sweep bit-width B: quantize the one-sided FFT coefficients
 *      (one shared peak scale), measure MSE/SNR and package bytes.
 *   4. pick the largest B with B <= 2*L.
 *   5. compare against time-domain uniform quantization at the same
 *      package size.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define MAXN 8192
#define PI 3.14159265358979323846
#define E 2.71828182845904523536

static double log2d(double x) { return log(x) / 0.6931471805599453; }

/* iterative radix-2 FFT on separate re/im arrays */
static void fft(double* re, double* im, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double t;
            t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * PI / len * (inverse ? 1.0 : -1.0);
        double wr = cos(ang), wi = sin(ang);
        for (int i = 0; i < n; i += len) {
            double cr = 1.0, ci = 0.0;
            for (int k = 0; k < len / 2; k++) {
                double ur = re[i + k], ui = im[i + k];
                double vr = re[i + k + len / 2] * cr - im[i + k + len / 2] * ci;
                double vi = re[i + k + len / 2] * ci + im[i + k + len / 2] * cr;
                re[i + k] = ur + vr; im[i + k] = ui + vi;
                re[i + k + len / 2] = ur - vr; im[i + k + len / 2] = ui - vi;
                double nr = cr * wr - ci * wi, ni = cr * wi + ci * wr;
                cr = nr; ci = ni;
            }
        }
    }
    if (inverse) {
        for (int i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
    }
}

static uint32_t rng_state;
static void rng_seed(uint32_t s) { rng_state = s ? s : 1u; }
static double frand(void) { /* uniform in [-1,1) */
    rng_state = rng_state * 1664525u + 1013904223u;
    return (double)(rng_state >> 8) / 8388608.0 - 1.0;
}

/* game-SFX-like material: three decaying partials over a noise floor */
static void gen_material(double* x, int n, uint32_t seed) {
    static const double f[3] = { 0.0203, 0.0617, 0.1371 }; /* cycles/sample */
    static const double a[3] = { 1.0, 0.5, 0.25 };
    double tau = n / 4.0;
    rng_seed(seed);
    for (int i = 0; i < n; i++) {
        double env = exp(-i / tau);
        double s = 0.0;
        for (int k = 0; k < 3; k++) s += a[k] * sin(2.0 * PI * f[k] * i + 0.7 * k);
        x[i] = env * s + 0.02 * frand();
    }
}

typedef struct {
    double sigma2;  /* signal variance */
    double h_spec;  /* spectral entropy, bit */
    double k_eff;   /* 2^h_spec: effective occupied bins */
    int m;          /* one-sided bins */
    double h_rate;  /* Gaussian entropy-rate estimate, bit/sample */
} Source;

static Source analyze(const double* x, int n) {
    static double re[MAXN], im[MAXN];
    Source s;
    s.m = n / 2;
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;
    s.sigma2 = 0.0;
    for (int i = 0; i < n; i++) { double d = x[i] - mean; s.sigma2 += d * d; }
    s.sigma2 /= n;

    for (int i = 0; i < n; i++) { re[i] = x[i]; im[i] = 0.0; }
    fft(re, im, n, 0);

    double sum = 0.0, smax = 0.0;
    for (int k = 0; k < n; k++) {
        double p = re[k] * re[k] + im[k] * im[k];
        if (k < s.m) sum += p;
        if (p > smax) smax = p;
    }
    double h = 0.0;
    for (int k = 0; k < s.m; k++) {
        double p = (re[k] * re[k] + im[k] * im[k]) / sum;
        if (p > 0.0) h -= p * log2d(p);
    }
    s.h_spec = h;
    s.k_eff = pow(2.0, h);

    /* Gaussian entropy rate from the two-sided periodogram; the floor
       keeps empty bins from sending the estimate to -inf. */
    double floor = smax * 1e-12 / n;
    double acc = 0.0;
    for (int k = 0; k < n; k++) {
        double sk = (re[k] * re[k] + im[k] * im[k]) / n;
        if (sk < floor) sk = floor;
        acc += log2d(sk);
    }
    s.h_rate = 0.5 * log2d(2.0 * PI * E) + 0.5 * acc / n;
    return s;
}

/* uniform B-bit symmetric quantization of m values sharing one peak
   scale; returns mean squared error */
static double quant_block(const double* v, int m, int bits, double* out) {
    double peak = 1e-30;
    for (int i = 0; i < m; i++) if (fabs(v[i]) > peak) peak = fabs(v[i]);
    double qmax = (double)((1u << (bits - 1)) - 1);
    double mse = 0.0;
    for (int i = 0; i < m; i++) {
        long q = lround(v[i] / peak * qmax);
        if (q > (long)qmax) q = (long)qmax;
        if (q < -(long)qmax) q = -(long)qmax;
        double dq = (double)q / qmax * peak;
        out[i] = dq;
        double d = dq - v[i];
        mse += d * d;
    }
    return mse / m;
}

/* quantize the one-sided spectrum (n/2+1 complex coeffs, one shared
   peak) and reconstruct; returns per-sample MSE, sets *nbits */
static double spec_quant(const double* x, int n, int bits, long* nbits) {
    static double re[MAXN], im[MAXN], packed[MAXN + 2], qp[MAXN + 2];
    int m = n / 2 + 1;
    for (int i = 0; i < n; i++) { re[i] = x[i]; im[i] = 0.0; }
    fft(re, im, n, 0);
    for (int k = 0; k < m; k++) { packed[2 * k] = re[k]; packed[2 * k + 1] = im[k]; }
    quant_block(packed, 2 * m, bits, qp);
    for (int k = 0; k < m; k++) { re[k] = qp[2 * k]; im[k] = qp[2 * k + 1]; }
    im[0] = 0.0; im[n / 2] = 0.0; /* DC and Nyquist stay real */
    for (int k = m; k < n; k++) { re[k] = re[n - k]; im[k] = -im[n - k]; }
    fft(re, im, n, 1);
    double mse = 0.0;
    for (int i = 0; i < n; i++) { double d = re[i] - x[i]; mse += d * d; }
    *nbits = 2L * m * bits + 16; /* coefficients + peak header */
    return mse / n;
}

/* time-domain uniform quantization at the same bit-width */
static double time_quant(const double* x, int n, int bits, long* nbits) {
    static double xh[MAXN];
    double mse = quant_block(x, n, bits, xh);
    *nbits = (long)n * bits + 16;
    return mse;
}

static const int CAND[] = { 2, 3, 4, 5, 6, 8, 10, 12, 16 };
#define NCAND (int)(sizeof(CAND) / sizeof(CAND[0]))

typedef struct {
    int chosen;
    double h_spec;
    double snr_spec; /* SNR at the chosen bit-width */
    double snr_uni;  /* same-size uniform-quantization SNR */
    long bytes;
} Result;

static Result run_length(int n, uint32_t seed) {
    static double x[MAXN];
    Result r;
    gen_material(x, n, seed);
    Source s = analyze(x, n);
    double bound = s.h_spec; /* lower bound L, bit/coeff */
    printf("n=%d  sigma2=%.4f  H_spec=%.2f bit  K_eff=%.1f/%d (%.1f%%)  h_rate=%.2f bit/sample\n",
           n, s.sigma2, s.h_spec, s.k_eff, s.m, 100.0 * s.k_eff / s.m, s.h_rate);
    printf("lower bound L = H_spec = %.2f bit/coeff;  rule: max B with B <= 2L = %.2f\n",
           bound, 2.0 * bound);
    printf("  B | pkg_bytes | spec_snr_dB | uni_snr_dB |  B/L\n");
    r.chosen = -1;
    r.h_spec = s.h_spec;
    r.snr_spec = 0.0; r.snr_uni = 0.0; r.bytes = 0;
    for (int i = 0; i < NCAND; i++) {
        int b = CAND[i];
        long sb, tb;
        double ds = spec_quant(x, n, b, &sb);
        double dt = time_quant(x, n, b, &tb);
        double snrs = 10.0 * log10(s.sigma2 / ds);
        double snrt = 10.0 * log10(s.sigma2 / dt);
        printf(" %2d | %9ld | %11.2f | %10.2f | %.2f%s\n",
               b, sb / 8, snrs, snrt, b / bound,
               b <= 2.0 * bound ? "" : "  (over 2x)");
        if (b <= 2.0 * bound) {
            r.chosen = b;
            r.snr_spec = snrs;
            r.snr_uni = snrt;
            r.bytes = sb / 8;
        }
    }
    printf("chosen B=%d  pkg=%ld bytes  spec_snr=%.2f dB  uni_snr=%.2f dB\n\n",
           r.chosen, r.bytes, r.snr_spec, r.snr_uni);
    return r;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* fft round trip */
    static double re[MAXN], im[MAXN], orig[MAXN];
    for (int i = 0; i < 256; i++) { re[i] = (i % 7) - 3.0; im[i] = 0.0; orig[i] = re[i]; }
    fft(re, im, 256, 0);
    fft(re, im, 256, 1);
    double maxd = 0.0;
    for (int i = 0; i < 256; i++) if (fabs(re[i] - orig[i]) > maxd) maxd = fabs(re[i] - orig[i]);
    check("fft_roundtrip", maxd < 1e-9);

    /* spectral entropy ordering: tone < material < white noise */
    static double t[MAXN], w[MAXN];
    for (int i = 0; i < 1024; i++) t[i] = sin(2.0 * PI * 0.0617 * i);
    rng_seed(99);
    for (int i = 0; i < 1024; i++) w[i] = frand();
    static double m[MAXN];
    gen_material(m, 1024, 7);
    Source st = analyze(t, 1024), sw = analyze(w, 1024), sm = analyze(m, 1024);
    printf("H_spec: tone=%.2f material=%.2f noise=%.2f\n", st.h_spec, sm.h_spec, sw.h_spec);
    check("entropy_order", st.h_spec < sm.h_spec && sm.h_spec < sw.h_spec);

    /* two lengths */
    Result r1 = run_length(1024, 7);
    Result r2 = run_length(8192, 11);

    /* selection must respect the 2x-bound rule */
    check("chosen_within_2x_n1024", r1.chosen > 0 && r1.chosen <= 2.0 * r1.h_spec + 1e-9);
    check("chosen_within_2x_n8192", r2.chosen > 0 && r2.chosen <= 2.0 * r2.h_spec + 1e-9);

    /* 16-bit spectral quantization must still reconstruct decently,
       and same-size uniform quantization must beat it on this peaky
       material (the gap the report has to explain) */
    long nb;
    static double x[MAXN];
    gen_material(x, 1024, 7);
    Source s = analyze(x, 1024);
    double d16s = spec_quant(x, 1024, 16, &nb);
    double d16t = time_quant(x, 1024, 16, &nb);
    double snr16s = 10.0 * log10(s.sigma2 / d16s);
    double snr16t = 10.0 * log10(s.sigma2 / d16t);
    printf("snr16: spec=%.2f uni=%.2f\n", snr16s, snr16t);
    check("spec16_sane", snr16s > 30.0);
    check("uni_beats_spec", snr16t > snr16s + 6.0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
