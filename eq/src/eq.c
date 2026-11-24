
/* eq: 8-band peaking EQ. Neighbor band gains interact through cross-coupling
   terms in the biquad coefficients; per-band suppression boxes bound the
   response deformation of bands squeezed by their neighbors, and a 30-group
   sweep of adjacent-triplet gain peaks flags configs whose intermodulation
   distortion exceeds a limit. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SR 48000.0f
#define NBANDS 8
#define NGRID 97
#define FMIN 30.0f
#define FMAX 16000.0f
#define TRIPLETS (NBANDS - 2)                 /* adjacent 3-band windows */
#define NLEVELS 5
#define NGROUPS (TRIPLETS * NLEVELS)          /* 6 triplets x 5 peaks = 30 */
#define IMD_LIMIT_DB 1.6f                     /* forbidden above this */

static const float BAND_F[NBANDS] = {
    60.0f, 120.0f, 240.0f, 480.0f, 960.0f, 1920.0f, 3840.0f, 7680.0f
};
static const float PEAK_DB[NLEVELS] = { 6.0f, 9.0f, 12.0f, 15.0f, 18.0f };

typedef struct { float b0, b1, b2, a1, a2; } Biquad;

typedef struct {
    float q;
    float couple;        /* cross-coupling strength between band gains */
    float couple_width;  /* coupling kernel width in octaves */
} EqCfg;

typedef struct { float lo, hi; } Box;         /* suppression box bounds, dB */

typedef struct {
    int band0;           /* first band of the adjacent triplet */
    float peak_db;       /* gain peak applied to the triplet */
    float imd_db;        /* worst deviation from the uncoupled ideal */
    int forbidden;
} Group;

static float grid_f(int i) {
    return FMIN * powf(FMAX / FMIN, (float)i / (float)(NGRID - 1));
}

/* RBJ peaking biquad */
static Biquad peaking(float f0, float q, float gain_db) {
    float A = powf(10.0f, gain_db / 40.0f);
    float w0 = 2.0f * (float)M_PI * f0 / SR;
    float alpha = sinf(w0) / (2.0f * q);
    float cw = cosf(w0);
    float a0 = 1.0f + alpha / A;
    Biquad b;
    b.b0 = (1.0f + alpha * A) / a0;
    b.b1 = (-2.0f * cw) / a0;
    b.b2 = (1.0f - alpha * A) / a0;
    b.a1 = b.b1;
    b.a2 = (1.0f - alpha / A) / a0;
    return b;
}

/* cross-coupling terms: each band's effective gain is eaten into by its
   neighbors, weighted by a kernel over octave distance (gain competition):
   geff[i] = g[i] - couple * sum_j exp(-oct(i,j)/width) * g[j] */
static void effective_gains(const EqCfg* c, const float* g, float* geff) {
    for (int i = 0; i < NBANDS; i++) {
        float e = g[i];
        for (int j = 0; j < NBANDS; j++) {
            if (j == i) continue;
            float oct = fabsf(log2f(BAND_F[i] / BAND_F[j]));
            e -= c->couple * expf(-oct / c->couple_width) * g[j];
        }
        geff[i] = e;
    }
}

static float biquad_db(const Biquad* b, float f) {
    float w = 2.0f * (float)M_PI * f / SR;
    float c1 = cosf(w), s1 = sinf(w);
    float c2 = cosf(2.0f * w), s2 = sinf(2.0f * w);
    float nr = b->b0 + b->b1 * c1 + b->b2 * c2;
    float ni = -(b->b1 * s1 + b->b2 * s2);
    float dr = 1.0f + b->a1 * c1 + b->a2 * c2;
    float di = -(b->a1 * s1 + b->a2 * s2);
    return 10.0f * log10f((nr * nr + ni * ni) / (dr * dr + di * di));
}

static float cascade_db(const EqCfg* c, const float* gains, float f) {
    float geff[NBANDS];
    effective_gains(c, gains, geff);
    float db = 0.0f;
    for (int i = 0; i < NBANDS; i++) {
        Biquad b = peaking(BAND_F[i], c->q, geff[i]);
        db += biquad_db(&b, f);
    }
    return db;
}

/* ideal: bands fully independent, coupling disabled */
static float ideal_db(const EqCfg* c, const float* gains, float f) {
    EqCfg nc = *c;
    nc.couple = 0.0f;
    return cascade_db(&nc, gains, f);
}

/* deviation of the coupled cascade from the independent-band ideal */
static void compute_deviation(const EqCfg* c, const float* gains, float* dev) {
    for (int i = 0; i < NGRID; i++) {
        float f = grid_f(i);
        dev[i] = cascade_db(c, gains, f) - ideal_db(c, gains, f);
    }
}

/* intermodulation distortion: worst-case deviation over the band plot */
static float imd_db(const EqCfg* c, const float* gains) {
    float dev[NGRID], m = 0.0f;
    compute_deviation(c, gains, dev);
    for (int i = 0; i < NGRID; i++)
        if (fabsf(dev[i]) > m) m = fabsf(dev[i]);
    return m;
}

/* suppression box per band: min/max deviation over the band's half-octave
   span, i.e. the bounds of the deformation inflicted by the neighbors */
static void suppression_boxes(const EqCfg* c, const float* gains, Box* boxes) {
    float dev[NGRID];
    compute_deviation(c, gains, dev);
    for (int b = 0; b < NBANDS; b++) {
        float flo = BAND_F[b] / sqrtf(2.0f);
        float fhi = BAND_F[b] * sqrtf(2.0f);
        float lo = 1e30f, hi = -1e30f;
        for (int i = 0; i < NGRID; i++) {
            float f = grid_f(i);
            if (f < flo || f > fhi) continue;
            if (dev[i] < lo) lo = dev[i];
            if (dev[i] > hi) hi = dev[i];
        }
        boxes[b].lo = lo;
        boxes[b].hi = hi;
    }
}

/* 30 groups: 6 adjacent triplets x 5 gain-peak levels */
static int sweep(const EqCfg* c, Group* groups) {
    int n = 0;
    for (int t = 0; t < TRIPLETS; t++) {
        for (int l = 0; l < NLEVELS; l++) {
            float g[NBANDS] = { 0.0f };
            g[t] = g[t + 1] = g[t + 2] = PEAK_DB[l];
            Group* gr = &groups[n++];
            gr->band0 = t;
            gr->peak_db = PEAK_DB[l];
            gr->imd_db = imd_db(c, g);
            gr->forbidden = gr->imd_db > IMD_LIMIT_DB;
        }
    }
    return n;
}

#define ROWS 24
#define COLS NGRID

static int row_of(float db, float vmin, float vmax) {
    int r = (int)((vmax - db) / (vmax - vmin) * (float)(ROWS - 1) + 0.5f);
    if (r < 0) r = 0;
    if (r > ROWS - 1) r = ROWS - 1;
    return r;
}

/* band plot with suppression boxes: '#' coupled response, '.' uncoupled
   ideal, '-'/'|' the per-band deformation bounds */
static void plot(const EqCfg* c, const float* gains, const Box* boxes) {
    static char canvas[ROWS][COLS + 1];
    float act[NGRID], idl[NGRID];
    float vmin = 1e30f, vmax = -1e30f;
    for (int i = 0; i < NGRID; i++) {
        float f = grid_f(i);
        act[i] = cascade_db(c, gains, f);
        idl[i] = ideal_db(c, gains, f);
        if (act[i] < vmin) vmin = act[i];
        if (act[i] > vmax) vmax = act[i];
        if (idl[i] < vmin) vmin = idl[i];
        if (idl[i] > vmax) vmax = idl[i];
    }
    for (int b = 0; b < NBANDS; b++) {
        if (boxes[b].lo < vmin) vmin = boxes[b].lo;
        if (boxes[b].hi > vmax) vmax = boxes[b].hi;
    }
    vmin -= 0.5f;
    vmax += 0.5f;

    for (int r = 0; r < ROWS; r++) {
        memset(canvas[r], ' ', COLS);
        canvas[r][COLS] = '\0';
    }
    /* suppression boxes first so the curves draw over them */
    for (int b = 0; b < NBANDS; b++) {
        float flo = BAND_F[b] / sqrtf(2.0f);
        float fhi = BAND_F[b] * sqrtf(2.0f);
        int c0 = COLS - 1, c1 = 0;
        for (int i = 0; i < NGRID; i++) {
            float f = grid_f(i);
            if (f >= flo && i < c0) c0 = i;
            if (f <= fhi && i > c1) c1 = i;
        }
        int rlo = row_of(boxes[b].lo, vmin, vmax);
        int rhi = row_of(boxes[b].hi, vmin, vmax);
        for (int col = c0; col <= c1; col++) {
            canvas[rlo][col] = '-';
            canvas[rhi][col] = '-';
        }
        for (int r = rhi; r <= rlo; r++) {
            canvas[r][c0] = '|';
            canvas[r][c1] = '|';
        }
    }
    for (int i = 0; i < NGRID; i++) {
        canvas[row_of(idl[i], vmin, vmax)][i] = '.';
        canvas[row_of(act[i], vmin, vmax)][i] = '#';
    }
    printf("band plot (dB, %.0f Hz..%.0f Hz): # actual, . ideal, box = suppression bounds\n",
           FMIN, FMAX);
    for (int r = 0; r < ROWS; r++) {
        float db = vmax - (vmax - vmin) * (float)r / (float)(ROWS - 1);
        printf("%+7.2f |%s|\n", db, canvas[r]);
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    EqCfg cfg = { 1.4f, 0.35f, 0.5f };
    float g[NBANDS];

    memset(g, 0, sizeof(g));
    check("flat_unity", fabsf(cascade_db(&cfg, g, 1000.0f)) < 1e-3f);
    check("flat_no_imd", imd_db(&cfg, g) < 1e-3f);

    memset(g, 0, sizeof(g));
    g[3] = 12.0f;
    check("single_band_peak", fabsf(cascade_db(&cfg, g, BAND_F[3]) - 12.0f) < 0.5f);

    /* triplet gain peak: coupling off -> no intermodulation, on -> distortion */
    EqCfg nc = cfg;
    nc.couple = 0.0f;
    memset(g, 0, sizeof(g));
    g[2] = g[3] = g[4] = 18.0f;
    check("no_coupling_no_imd", imd_db(&nc, g) < 1e-3f);
    check("coupling_imd_present", imd_db(&cfg, g) > 0.5f);

    float g6[NBANDS] = { 0.0f };
    g6[2] = g6[3] = g6[4] = 6.0f;
    check("imd_grows_with_gain", imd_db(&cfg, g) > imd_db(&cfg, g6));

    /* 30-group sweep over adjacent triplets x gain peaks */
    Group groups[NGROUPS];
    int n = sweep(&cfg, groups);
    check("groups_30", n == 30);

    int worst = 0, nforbidden = 0;
    for (int i = 0; i < n; i++) {
        if (groups[i].imd_db > groups[worst].imd_db) worst = i;
        if (groups[i].forbidden) nforbidden++;
    }
    check("worst_is_max_peak", groups[worst].peak_db == PEAK_DB[NLEVELS - 1]);
    check("forbidden_some_not_all", nforbidden > 0 && nforbidden < NGROUPS);

    float gw[NBANDS] = { 0.0f };
    int t0 = groups[worst].band0;
    gw[t0] = gw[t0 + 1] = gw[t0 + 2] = groups[worst].peak_db;
    Box boxes[NBANDS];
    suppression_boxes(&cfg, gw, boxes);
    int ordered = 1;
    for (int b = 0; b < NBANDS; b++)
        if (!(boxes[b].lo <= boxes[b].hi + 1e-6f)) ordered = 0;
    check("boxes_ordered", ordered);

    printf("imd_limit_db=%.2f groups=%d forbidden=%d\n", (double)IMD_LIMIT_DB, n, nforbidden);
    printf("%-6s %-12s %-8s %-8s %s\n", "group", "bands_hz", "peak_db", "imd_db", "status");
    for (int i = 0; i < n; i++) {
        Group* gr = &groups[i];
        printf("%-6d %.0f-%-8.0f %-8.1f %-8.3f %s\n", i,
               (double)BAND_F[gr->band0], (double)BAND_F[gr->band0 + 2],
               (double)gr->peak_db, (double)gr->imd_db,
               gr->forbidden ? "FORBIDDEN" : "ok");
    }
    printf("worst_group=%d bands=%.0f-%.0f peak_db=%.1f imd_db=%.3f\n", worst,
           (double)BAND_F[t0], (double)BAND_F[t0 + 2],
           (double)groups[worst].peak_db, (double)groups[worst].imd_db);

    printf("suppression boxes (worst group):\n");
    for (int b = 0; b < NBANDS; b++)
        printf("band %d (%.0f Hz): [%+.3f, %+.3f] dB\n", b,
               (double)BAND_F[b], (double)boxes[b].lo, (double)boxes[b].hi);
    plot(&cfg, gw, boxes);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
