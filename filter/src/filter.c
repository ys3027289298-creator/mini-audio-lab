
/* filter: 4th-order bandpass (Butterworth HP + LP biquad cascade).
 *
 * Passband flatness relations (derivation in README.md):
 *   single 2nd-order stage with Q_s > 1/sqrt(2) resonates:
 *     G_peak(Q_s) = 2 Q_s^2 / sqrt(4 Q_s^2 - 1),  G_peak(1/sqrt(2)) = 1
 *     LP peaks at f0*sqrt(1 - 1/(2 Q_s^2)), HP at f0/sqrt(1 - 1/(2 Q_s^2)).
 *   Setting stage Q = band Q (the bug) makes the HP stage resonate at fL
 *   and lifts the low end of the passband at high Q. Fix: stage Q is
 *   always 1/sqrt(2); band Q only sets the cutoffs:
 *     fL = fc*(sqrt(1+1/(4Q^2)) - 1/(2Q)), fH = fc*(sqrt(1+1/(4Q^2)) + 1/(2Q))
 *   composite passband peak (at fc, with bilinear frequency warping):
 *     G_pb(Q) = wH^2 / (wL^2 + wH^2),  w = tan(pi f / fs)
 *
 * Self-test: validates the stage relation, reproduces the low-end lift,
 * runs 24 (fc, Q) groups end-to-end (impulse -> filter -> FFT) against
 * G_pb and flags any group whose measured flatness deviates, then chains
 * the limiter and checks that high-Q overshoot never engages it.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../../fft/src/fft.h"
#include "../../limiter/src/limiter.h"

#define FS 48000.0f
#define NFFT 32768
#define BLK 128
#define TONE_N 4864
#define PI 3.14159265358979
#define SQRT1_2 0.70710678118654752f

typedef struct { float b0, b1, b2, a1, a2, z1, z2; } Biquad;

static Biquad biquad_lp(float f0, float q, float fs) {
    float w = 2.0f * (float)PI * f0 / fs;
    float alpha = sinf(w) / (2.0f * q);
    float c = cosf(w);
    float a0 = 1.0f + alpha;
    Biquad b;
    b.b0 = (1.0f - c) * 0.5f / a0;
    b.b1 = (1.0f - c) / a0;
    b.b2 = b.b0;
    b.a1 = -2.0f * c / a0;
    b.a2 = (1.0f - alpha) / a0;
    b.z1 = b.z2 = 0.0f;
    return b;
}

static Biquad biquad_hp(float f0, float q, float fs) {
    float w = 2.0f * (float)PI * f0 / fs;
    float alpha = sinf(w) / (2.0f * q);
    float c = cosf(w);
    float a0 = 1.0f + alpha;
    Biquad b;
    b.b0 = (1.0f + c) * 0.5f / a0;
    b.b1 = -(1.0f + c) / a0;
    b.b2 = b.b0;
    b.a1 = -2.0f * c / a0;
    b.a2 = (1.0f - alpha) / a0;
    b.z1 = b.z2 = 0.0f;
    return b;
}

static float biquad_step(Biquad* b, float x) {
    float y = b->b0 * x + b->z1;
    b->z1 = b->b1 * x - b->a1 * y + b->z2;
    b->z2 = b->b2 * x - b->a2 * y;
    return y;
}

/* analytic magnitude of one biquad at frequency f */
static float biquad_mag(const Biquad* b, float f, float fs) {
    double w = 2.0 * PI * f / fs;
    double c1 = cos(w), s1 = sin(w), c2 = cos(2.0 * w), s2 = sin(2.0 * w);
    double nr = b->b0 + b->b1 * c1 + b->b2 * c2;
    double ni = -(b->b1 * s1 + b->b2 * s2);
    double dr = 1.0 + b->a1 * c1 + b->a2 * c2;
    double di = -(b->a1 * s1 + b->a2 * s2);
    return (float)(sqrt(nr * nr + ni * ni) / sqrt(dr * dr + di * di));
}

typedef struct { Biquad hp, lp; float fL, fH; } Bandpass;

/* buggy: stage Q = band Q (resonates, lifts passband edges)
 * fixed: stage Q = 1/sqrt(2) (maximally flat Butterworth stages) */
static void bp_init(Bandpass* bp, float fc, float q, float fs, int buggy) {
    double s = sqrt(1.0 + 1.0 / (4.0 * (double)q * q));
    double t = 1.0 / (2.0 * (double)q);
    bp->fL = (float)(fc * (s - t));
    bp->fH = (float)(fc * (s + t));
    float qs = buggy ? q : SQRT1_2;
    bp->hp = biquad_hp(bp->fL, qs, fs);
    bp->lp = biquad_lp(bp->fH, qs, fs);
}

static float bp_step(Bandpass* bp, float x) {
    return biquad_step(&bp->lp, biquad_step(&bp->hp, x));
}

/* predicted composite passband peak from the flatness relation */
static double g_pb(double fL, double fH, double fs) {
    double wL = tan(PI * fL / fs), wH = tan(PI * fH / fs);
    return wH * wH / (wL * wL + wH * wH);
}

/* end-to-end magnitude response: impulse -> filter -> FFT */
static FftC g_buf[NFFT];
static float g_mag[NFFT / 2];
static float g_mag2[NFFT / 2];

static void measure_response(Bandpass* bp, float* mag) {
    memset(g_buf, 0, sizeof(g_buf));
    for (int i = 0; i < NFFT; i++)
        g_buf[i].re = bp_step(bp, i == 0 ? 1.0f : 0.0f);
    fft(g_buf, NFFT, 0);
    for (int k = 0; k < NFFT / 2; k++)
        mag[k] = sqrtf(g_buf[k].re * g_buf[k].re + g_buf[k].im * g_buf[k].im);
}

static float mag_at(const float* mag, float f) {
    int k = (int)(f * NFFT / FS + 0.5f);
    return mag[k];
}

/* filter -> limiter chain on a tone burst; returns limited sample count */
static int run_tone_chain(Bandpass* bp, float fc, float amp, float ceiling,
                          float* peak_in, float* peak_out) {
    Limiter l;
    lim_init(&l, 1.0f, 50.0f, FS);
    int limited = 0;
    float pin = 0.0f, pout = 0.0f;
    for (int s = 0; s < TONE_N; s += BLK) {
        float in[BLK], filt[BLK], out[BLK];
        Report rep;
        for (int k = 0; k < BLK; k++)
            in[k] = amp * sinf(2.0f * (float)PI * fc * (s + k) / FS);
        for (int k = 0; k < BLK; k++)
            filt[k] = bp_step(bp, in[k]);
        lim_process(&l, filt, out, BLK, ceiling, &rep);
        limited += rep.limited_samples;
        if (rep.peak_in > pin) pin = rep.peak_in;
        if (rep.peak_out > pout) pout = rep.peak_out;
    }
    *peak_in = pin;
    *peak_out = pout;
    return limited;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static float db(float g) { return 20.0f * log10f(g); }

int main(void) {
    /* 1. stage relation: G_peak(Q_s) = 2 Q_s^2 / sqrt(4 Q_s^2 - 1) */
    int rel_ok = 1;
    const float qrel[5] = { 1.0f, 2.0f, 4.0f, 8.0f, 16.0f };
    for (int i = 0; i < 5; i++) {
        float q = qrel[i];
        Biquad b = biquad_lp(1000.0f, q, FS);
        float best = 0.0f;
        for (int k = 0; k < 8000; k++) {
            float f = 500.0f * powf(4.0f, (float)k / 7999.0f);
            float m = biquad_mag(&b, f, FS);
            if (m > best) best = m;
        }
        float pred = 2.0f * q * q / sqrtf(4.0f * q * q - 1.0f);
        if (fabsf(db(best) - db(pred)) > 0.05f) rel_ok = 0;
    }
    check("stage_peak_relation", rel_ok);

    /* 2. bug demo: stage Q = band Q lifts the low end at high Q */
    Bandpass buggy, fixed;
    bp_init(&buggy, 1000.0f, 8.0f, FS, 1);
    bp_init(&fixed, 1000.0f, 8.0f, FS, 0);
    float fhp = buggy.fL / sqrtf(1.0f - 1.0f / (2.0f * 8.0f * 8.0f));
    measure_response(&buggy, g_mag);
    measure_response(&fixed, g_mag2);
    float gb = mag_at(g_mag, fhp), gf = mag_at(g_mag2, fhp);
    printf("low_end f=%.1fHz buggy=%.2fdB fixed=%.2fdB\n", fhp, db(gb), db(gf));
    check("buggy_low_end_lift", db(gb) - db(gf) > 6.0f);
    check("fixed_no_lift", gf <= 1.02f);

    /* 3. 24 (fc, Q) groups: measured passband flatness vs G_pb relation */
    const float fcs[6] = { 200.0f, 500.0f, 1000.0f, 4000.0f, 8000.0f, 12000.0f };
    const float qs[4] = { 0.7f, 2.0f, 6.0f, 16.0f };
    int flagged = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            Bandpass bp;
            bp_init(&bp, fcs[i], qs[j], FS, 0);
            measure_response(&bp, g_mag);
            float gmax = 0.0f;
            for (int k = 0; k < NFFT / 2; k++) {
                float f = (float)k * FS / NFFT;
                if (f >= bp.fL && f <= bp.fH && g_mag[k] > gmax) gmax = g_mag[k];
            }
            float pred = (float)g_pb(bp.fL, bp.fH, FS);
            float dev = fabsf(db(gmax) - db(pred));
            int ok = dev <= 0.35f && gmax <= 1.02f;
            printf("group fc=%-6g Q=%-4g flat_meas=%6.2fdB flat_pred=%6.2fdB %s\n",
                   fcs[i], qs[j], db(gmax), db(pred), ok ? "ok" : "FLAG");
            if (!ok) flagged++;
        }
    }
    check("groups_match_relation", flagged == 0);

    /* 4. limiter chain: high-Q overshoot must not engage the limiter */
    int trips = 0;
    float worst_in = 0.0f, worst_out = 0.0f;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            Bandpass bp;
            bp_init(&bp, fcs[i], qs[j], FS, 0);
            float pin, pout;
            int limited = run_tone_chain(&bp, fcs[i], 0.7f, 0.9f, &pin, &pout);
            if (pin > worst_in) worst_in = pin;
            if (pout > worst_out) worst_out = pout;
            if (limited > 0) {
                trips++;
                printf("FLAG limiter engaged fc=%g Q=%g limited=%d\n", fcs[i], qs[j], limited);
            }
        }
    }
    printf("limiter_worst_peak_in=%.4f out=%.4f ceiling=0.9\n", worst_in, worst_out);
    check("limiter_never_engaged", trips == 0);

    /* 5. sanity: the buggy coefficients would trip the same chain */
    float pin, pout;
    bp_init(&buggy, 1000.0f, 16.0f, FS, 1);
    int limited = run_tone_chain(&buggy, 1000.0f, 0.7f, 0.9f, &pin, &pout);
    printf("buggy_chain peak_in=%.2f limited=%d\n", pin, limited);
    check("buggy_would_trip_limiter", limited > 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
