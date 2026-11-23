
/* spatial: constant-power panning and inverse-distance attenuation. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

typedef struct { float left, right; } Stereo;

/* Equal power: theta = (x + 1) * pi / 4, L = cos(theta), R = sin(theta).
 * L^2 + R^2 = 1 for every x in [-1, 1]; the endpoint limits are
 * x -> -1: (L, R) -> (1, 0), x -> +1: (L, R) -> (0, 1).
 * The float angle path leaves ~1e-7 residue on the silent side, so the two
 * endpoints are pinned to exact +0.0f / 1.0f. Interior x is untouched and
 * still follows the curve. pan_set_exact_endpoints(0) restores the old
 * formula byte for byte; it is on by default. */
static int exact_endpoints = 1;

static void pan_set_exact_endpoints(int on) { exact_endpoints = on != 0; }

static int pan(float x, Stereo* s) {
    if (x < -1.0f || x > 1.0f) return -1;
    if (exact_endpoints) {
        if (x <= -1.0f) { s->left = 1.0f; s->right = 0.0f; return 0; }
        if (x >= 1.0f)  { s->left = 0.0f; s->right = 1.0f; return 0; }
    }
    float a = (x + 1.0f) * 0.25f * 3.14159265f;   /* 0..pi/2 */
    s->left = cosf(a);
    s->right = sinf(a);
    return 0;
}

static float attenuate(float dist, float ref, float rolloff) {
    if (dist < 0.0f || ref <= 0.0f || rolloff < 0.0f) return -1.0f;
    if (dist <= ref) return 1.0f;
    return ref / (ref + rolloff * (dist - ref));
}

static unsigned fbits(float f) {
    unsigned u;
    memcpy(&u, &f, sizeof u);
    return u;
}
static int bit_eq(float a, float b) { return fbits(a) == fbits(b); }
static int pos_zero(float f) { return fbits(f) == 0u; }
static double power(const Stereo* s) { return (double)s->left * s->left + (double)s->right * s->right; }

/* Old formula, kept verbatim in the test as the compatibility reference. */
static void pan_legacy(float x, Stereo* s) {
    float a = (x + 1.0f) * 0.25f * 3.14159265f;
    s->left = cosf(a);
    s->right = sinf(a);
}

/* Compare default-mode pan() against the old formula, bit for bit. */
static int same_as_legacy(int n) {
    for (int i = 0; i <= n; i++) {
        float x = -1.0f + 2.0f * (float)i / (float)n;
        Stereo a, b;
        if (pan(x, &a) != 0) return 0;
        pan_legacy(x, &b);
        if (!bit_eq(a.left, b.left) || !bit_eq(a.right, b.right)) return 0;
    }
    return 1;
}

/* Strict mode must change only the two endpoints. */
static int strict_sweep(int n) {
    for (int i = 0; i <= n; i++) {
        float x = -1.0f + 2.0f * (float)i / (float)n;
        Stereo a, b;
        if (pan(x, &a) != 0) return 0;
        if (x <= -1.0f) {
            if (!bit_eq(a.left, 1.0f) || !pos_zero(a.right)) return 0;
        } else if (x >= 1.0f) {
            if (!pos_zero(a.left) || !bit_eq(a.right, 1.0f)) return 0;
        } else {
            pan_legacy(x, &b);
            if (!bit_eq(a.left, b.left) || !bit_eq(a.right, b.right)) return 0;
        }
    }
    return 1;
}

static double amp_db(double a, double b) {
    if (a == b) return 0.0;
    if (a == 0.0) return -INFINITY;
    if (b == 0.0) return INFINITY;
    return 20.0 * log10(a / b);
}
static double pow_db(double a, double b) { return 10.0 * log10(a / b); }

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Stereo s, ref;

    /* original assertions, still green with endpoints pinned */
    check("center", pan(0.0f, &s) == 0 && fabsf(s.left - s.right) < 1e-6f);
    check("center_power", fabsf(s.left * s.left + s.right * s.right - 1.0f) < 1e-5f);

    pan(-1.0f, &s);
    check("hard_left", s.left > 0.99f && s.right < 1e-6f);

    pan(1.0f, &s);
    check("hard_right", s.right > 0.99f && s.left < 1e-6f);

    check("pan_bad", pan(1.5f, &s) == -1);

    check("near_full", fabsf(attenuate(1.0f, 2.0f, 1.0f) - 1.0f) < 1e-6f);
    check("farther_quieter", attenuate(5.0f, 2.0f, 1.0f) < attenuate(3.0f, 2.0f, 1.0f));
    check("bad_dist", attenuate(-1.0f, 2.0f, 1.0f) < 0.0f);

    /* feature off: bit-identical to the old implementation */
    pan_set_exact_endpoints(0);
    pan(1.0f, &s);
    pan_legacy(1.0f, &ref);
    check("legacy_off_bits", bit_eq(s.left, ref.left) && bit_eq(s.right, ref.right));
    check("legacy_off_dense_bits", same_as_legacy(100000));

    /* feature on: endpoints exact, interior curve unchanged */
    pan_set_exact_endpoints(1);
    Stereo l, r;
    pan(-1.0f, &l);
    pan(1.0f, &r);
    check("strict_left_endpoint", bit_eq(l.left, 1.0f) && pos_zero(l.right));
    check("strict_right_endpoint", pos_zero(r.left) && bit_eq(r.right, 1.0f));
    check("strict_endpoint_power", power(&l) == 1.0 && power(&r) == 1.0);
    check("strict_pan_bad", pan(-1.5f, &s) == -1 && pan(1.5f, &s) == -1);

    /* curve sample points: energy stays 1, interior bits match legacy */
    static const float xs[] = {-1.0f, -0.75f, -0.5f, -0.25f, 0.0f,
                               0.25f, 0.5f, 0.75f, 1.0f};
    const int nx = (int)(sizeof xs / sizeof xs[0]);
    int interior_same = 1, curve_power = 1;
    for (int i = 0; i < nx; i++) {
        double x = (double)xs[i];
        double lin_l = (1.0 - x) * 0.5, lin_r = (1.0 + x) * 0.5;
        double lin_e = lin_l * lin_l + lin_r * lin_r;
        pan(xs[i], &s);
        pan_legacy(xs[i], &ref);
        if (xs[i] > -1.0f && xs[i] < 1.0f &&
            (!bit_eq(s.left, ref.left) || !bit_eq(s.right, ref.right))) interior_same = 0;
        if (fabs(power(&s) - 1.0) > 1e-6) curve_power = 0;
        printf("curve x=%+.2f ep_L=%.9f ep_R=%.9f ep_E=%.12f "
               "lin_L=%.9f lin_R=%.9f lin_E=%.9f "
               "dL_dB=%+.4f dR_dB=%+.4f dE_dB=%+.4f\n",
               x, (double)s.left, (double)s.right, power(&s),
               lin_l, lin_r, lin_e,
               amp_db((double)s.left, lin_l), amp_db((double)s.right, lin_r),
               pow_db(power(&s), lin_e));
    }
    check("strict_curve_bits", interior_same);
    check("strict_curve_power", curve_power);
    check("strict_sweep_only_endpoints", strict_sweep(100000));

    /* endpoint limit: residual falls toward zero as x -> -1 */
    static const float ds[] = {1e-1f, 1e-2f, 1e-3f, 1e-4f, 1e-5f, 1e-6f, 1e-7f};
    double prev = 1.0;
    int limit_ok = 1;
    for (int i = 0; i < (int)(sizeof ds / sizeof ds[0]); i++) {
        pan(-1.0f + ds[i], &s);
        double res = (double)s.right;
        if (!(res > 0.0 && res < prev)) limit_ok = 0;
        if (fabs(power(&s) - 1.0) > 1e-6) limit_ok = 0;
        prev = res;
    }
    pan(-1.0f, &s);
    check("endpoint_limit", limit_ok && pos_zero(s.right) && bit_eq(s.left, 1.0f));

    /* stereo waveform at the endpoints: loud side keeps input bits, silent side is zero */
    enum { NW = 64 };
    float in[NW], ol[NW], or_[NW];
    for (int i = 0; i < NW; i++) in[i] = (float)((i * 37) % 101 - 50) / 50.0f;
    int loud_bits = 1, silent_zero = 1;
    pan(-1.0f, &s);
    for (int i = 0; i < NW; i++) {
        ol[i] = in[i] * s.left;
        or_[i] = in[i] * s.right;
        if (!bit_eq(ol[i], in[i])) loud_bits = 0;
        if (or_[i] != 0.0f) silent_zero = 0;
    }
    check("wave_left_loud_bits", loud_bits);
    check("wave_left_silent_zero", silent_zero);

    loud_bits = 1; silent_zero = 1;
    pan(1.0f, &s);
    for (int i = 0; i < NW; i++) {
        ol[i] = in[i] * s.left;
        or_[i] = in[i] * s.right;
        if (!bit_eq(or_[i], in[i])) loud_bits = 0;
        if (ol[i] != 0.0f) silent_zero = 0;
    }
    check("wave_right_loud_bits", loud_bits);
    check("wave_right_silent_zero", silent_zero);

    /* numeric gap to linear panning */
    pan(0.0f, &s);
    check("vs_linear_center_channel", fabs(amp_db((double)s.left, 0.5) - 3.0103) < 1e-3);
    pan(0.5f, &s);
    check("vs_linear_half_left", fabs(amp_db((double)s.left, 0.25) - 3.6980) < 1e-3);
    check("vs_linear_half_right", fabs(amp_db((double)s.right, 0.75) - 1.8111) < 1e-3);
    check("vs_linear_half_power", fabs(pow_db(power(&s), 0.625) - 2.0412) < 1e-3);

    /* turning the feature off restores the old formula byte for byte */
    pan_set_exact_endpoints(0);
    check("strict_off_restores_bits", same_as_legacy(100000));
    pan_set_exact_endpoints(1);

    /* default state is the fixed one. Do this last so the endpoint checks
       above run through the same public entry point callers use. */
    {
        Stereo dl, dr;
        pan(-1.0f, &dl);
        pan(1.0f, &dr);
        check("strict_on_by_default",
              pos_zero(dl.right) && bit_eq(dl.left, 1.0f) &&
              pos_zero(dr.left) && bit_eq(dr.right, 1.0f));
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
