
/* hop: hop-size selection.
 *
 * Frame shift h trades precision against compute:
 *   L(h) ~ alpha * h^2  (linear-interpolation RMSE, same interp as resampler)
 *   C(h) ~ T / h        (number of frames, each one NFFT-point FFT)
 * The composite cost is strictly convex (see README for the proof), so the
 * minimizer is located in closed form instead of enumerating candidate hops.
 * This program builds the two functions, measures them at three bracketing
 * points only, and checks the analytic minimizer against the measurements.
 */
#include <stdio.h>
#include <math.h>

#define PI   3.14159265358979
#define T    1024      /* signal length, samples                         */
#define NFFT 64        /* frame size used by the fft module              */
#define F0   64        /* frame budget: frames allowed at reference hop  */
#define TAU  0.03      /* RMSE tolerance on a unit-amplitude signal (3%) */
#define HMIN 4         /* allowed hop interval: 50%..~94% overlap for    */
#define HMAX 32        /* NFFT=64, h <= N/2 keeps the overlap-add valid  */

static double signal_at(int i) {
    double w = 2.0 * PI * 0.01;          /* f = 0.01 cycles/sample */
    return cos(w * (double)i);
}

/* Reconstruction from frame anchors at positions j*h, linear interpolation
 * between adjacent anchors (same convention as resampler_linear). */
static double reconstruct(const double* x, int i, int h) {
    int j = i / h;
    int u = j * h;
    int v = u + h;
    if (v > T - 1) v = T - 1;
    double t = (double)(i - u) / (double)h;
    return (1.0 - t) * x[u] + t * x[v];
}

/* Frames are the anchors j*h with j*h < T, hence K = ceil(T/h). */
static int frame_count(int h) { return (T + h - 1) / h; }

static int frame_count_loop(int h) {
    int k = 0;
    for (int j = 0; j * h < T; j++) k++;
    return k;
}

typedef struct {
    int h;
    double L;      /* measured precision loss: RMSE */
    int    K;      /* measured frame count          */
    double J;      /* measured composite cost       */
} Meas;

static Meas measure(int h, double tau, int f0) {
    double x[T];
    for (int i = 0; i < T; i++) x[i] = signal_at(i);
    double se = 0.0;
    for (int i = 0; i < T; i++) {
        double e = x[i] - reconstruct(x, i, h);
        se += e * e;
    }
    Meas m;
    m.h = h;
    /* measured RMSE; model says L(h) = alpha h^2 + O(h^4), see README */
    m.L = sqrt(se / (double)T);
    m.K = frame_count_loop(h);
    double p = m.L / tau;                 /* normalized precision term */
    double c = (double)m.K / (double)f0;  /* normalized compute term   */
    m.J = p * p + c * c;
    return m;
}

/* Model cost Jm(h) = a*h^4 + b/h^2 and its first derivative. */
static double Jm(double h, double a, double b) {
    return a * h * h * h * h + b / (h * h);
}
static double dJm(double h, double a, double b) {
    return 4.0 * a * h * h * h - 2.0 * b / (h * h * h);
}

/* Closed-form minimizer clipped to the allowed interval — no enumeration.
 * Strict convexity guarantees the integer optimum is one of floor(h*),
 * ceil(h*) (clipped to the interval), so at most two O(1) evaluations. */
static int optimal_hop(double a, double b, int hmin, int hmax) {
    double hs = pow(b / (2.0 * a), 1.0 / 6.0);
    int lo = (int)floor(hs);
    if (lo < hmin) lo = hmin;
    if (lo > hmax) lo = hmax;
    int hi = lo + 1;
    if (hi > hmax) hi = hmax;
    if (hi <= lo) return lo;
    return Jm((double)lo, a, b) <= Jm((double)hi, a, b) ? lo : hi;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* ---- model coefficients ------------------------------------------------
     * Leading interpolation error for a band-limited signal, anchors every h:
     *   e(t) ~= h^2/2 t(1-t) x''(u).  Averaging t over [0,1] and the signal
     *   curvature over time gives
     *   RMSE = sqrt(E[(x'')^2] / 120) * h^2 = alpha * h^2.
     * For x = cos(w i), E[(x'')^2] = w^4/2, hence alpha = w^2/sqrt(240).       */
    double w = 2.0 * PI * 0.01;
    double m2 = w * w;
    double alpha = m2 / sqrt(240.0);
    double a = (alpha / TAU) * (alpha / TAU);
    double b = ((double)T / (double)F0) * ((double)T / (double)F0);
    double hstar = pow(b / (2.0 * a), 1.0 / 6.0);

    printf("model: L(h) ~= %.6e h^2,  C(h) = ceil(%d/h)\n", alpha, T);
    printf("cost : J(h) = %.6e h^4 + %.4f / h^2\n", a, b);
    printf("analytic minimizer h* = (b/2a)^(1/6) = %.4f, interval [%d,%d]\n",
           hstar, HMIN, HMAX);
    printf("\n%-4s %-10s %-10s %-6s %-10s %-10s\n",
           "h", "L_meas", "alpha*h2", "K", "J_meas", "J_model");

    int pts[3] = { 9, 11, 13 };           /* bracket of h* ~= 11.0 */
    Meas mm[3];
    for (int q = 0; q < 3; q++) {
        mm[q] = measure(pts[q], TAU, F0);
        printf("%-4d %-10.6f %-10.6f %-6d %-10.5f %-10.5f\n",
               mm[q].h, mm[q].L, alpha * pts[q] * pts[q], mm[q].K,
               mm[q].J, Jm((double)pts[q], a, b));
    }

    int hop = optimal_hop(a, b, HMIN, HMAX);
    printf("\nconstrained integer minimizer = %d\n", hop);

    /* ---- checks ----------------------------------------------------------- */
    /* 1. frame-count formula K = ceil(T/h) matches an actual anchor loop */
    int k_ok = 1;
    for (int q = 0; q < 3; q++)
        k_ok &= (mm[q].K == frame_count(mm[q].h));
    check("frame_count_ceil", k_ok);

    /* 2. quadratic precision law at all three points (h^6 residual < 12%) */
    int q_ok = 1;
    for (int q = 0; q < 3; q++) {
        double ratio = mm[q].L / (alpha * mm[q].h * mm[q].h);
        q_ok &= fabs(ratio - 1.0) < 0.12;
    }
    check("quadratic_loss_law", q_ok);

    /* 3. three-point measured valley sits at the middle (bracketing) point */
    check("three_point_valley",
          mm[0].J > mm[1].J && mm[2].J > mm[1].J);

    /* 4. the measured valley point is the closed-form pick, no enumeration */
    check("pick_matches_valley", hop == mm[1].h);

    /* 5. analytic derivative signs: decreasing left of h*, increasing right */
    check("derivative_signs",
          dJm(9.0, a, b) < 0.0 && dJm(13.0, a, b) > 0.0);

    /* 6. strict convexity of the model cost on the whole interval */
    int convex_ok = 1;
    for (int h = HMIN; h <= HMAX; h++) {
        double d2 = 12.0 * a * h * h + 6.0 * b / pow((double)h, 4.0);
        convex_ok &= (d2 > 0.0);
    }
    check("strictly_convex", convex_ok);

    /* 7. interval clipping: a minimizer outside the interval hits the edge */
    check("clip_left_edge",  optimal_hop(a, b, 16, 32) == 16);
    check("clip_right_edge", optimal_hop(a, b, 4, 8) == 8);

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
