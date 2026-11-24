
/* reflect: image-source reflection paths with direction-aware stereo split.
 *
 * Each valid reflection path yields an arrival azimuth at the listener.
 * Path energy is split between left/right ears by a constant-power pan of
 * sin(azimuth) (same law as spatial.c) plus a far-ear head-shadow factor
 * driven by the interaural path difference d_ear * sin(azimuth).
 */
#include <stdio.h>
#include <math.h>

#define MAX_WALLS 8
#define MAX_ORDER 3
#define MAX_PATHS 512

#define ABSORP 0.8   /* amplitude kept per bounce */
#define D_EAR 0.18   /* interaural distance, m */
#define S0 0.10      /* head-shadow length scale, m */
#define PI 3.14159265358979323846
#define EPS 1e-9

typedef struct { double x, y; } V2;
typedef struct { V2 a, b; } Wall;
typedef struct { V2 p, fwd; } Listener;

typedef struct {
    double az;      /* arrival azimuth, rad, + = left of facing */
    double energy;
    int order;
} Path;

typedef struct {
    Wall walls[MAX_WALLS];
    int nwalls;
    V2 source;
    Listener lis;
    Path paths[MAX_PATHS];
    int npaths;
} Scene;

static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static double vdot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static double vcross(V2 a, V2 b) { return a.x * b.y - a.y * b.x; }
static double vlen(V2 a) { return sqrt(vdot(a, a)); }
static V2 vnorm(V2 a) { double l = vlen(a); V2 r = { a.x / l, a.y / l }; return r; }

/* mirror p across the wall's infinite line */
static V2 reflect_point(V2 p, Wall w) {
    V2 d = vnorm(vsub(w.b, w.a));
    V2 rel = vsub(p, w.a);
    double t = vdot(rel, d);
    V2 foot = { w.a.x + d.x * t, w.a.y + d.y * t };
    V2 r = { 2.0 * foot.x - p.x, 2.0 * foot.y - p.y };
    return r;
}

/* intersection of segment a->b with the wall segment; 1 if within both */
static int cross_point(Wall w, V2 a, V2 b, V2* hit) {
    V2 r = vsub(b, a), s = vsub(w.b, w.a);
    double den = vcross(r, s);
    if (fabs(den) < EPS) return 0;
    V2 d = vsub(w.a, a);
    double t = vcross(d, s) / den;   /* along a->b */
    double u = vcross(d, r) / den;   /* along wall */
    if (t < -EPS || t > 1.0 + EPS || u < -EPS || u > 1.0 + EPS) return 0;
    hit->x = a.x + r.x * t;
    hit->y = a.y + r.y * t;
    return 1;
}

/* proper crossing strictly inside both segments => leg is obstructed */
static int leg_blocked(Scene* sc, V2 a, V2 b) {
    V2 r = vsub(b, a);
    for (int i = 0; i < sc->nwalls; i++) {
        V2 s = vsub(sc->walls[i].b, sc->walls[i].a);
        double den = vcross(r, s);
        if (fabs(den) < EPS) continue;
        V2 d = vsub(sc->walls[i].a, a);
        double t = vcross(d, s) / den;
        double u = vcross(d, r) / den;
        if (t > 1e-6 && t < 1.0 - 1e-6 && u > 1e-6 && u < 1.0 - 1e-6)
            return 1;
    }
    return 0;
}

/* validate one wall sequence and record the path */
static void trace(Scene* sc, const int* seq, int n) {
    V2 img[MAX_ORDER], hit[MAX_ORDER];
    img[0] = reflect_point(sc->source, sc->walls[seq[0]]);
    for (int k = 1; k < n; k++)
        img[k] = reflect_point(img[k - 1], sc->walls[seq[k]]);

    V2 cur = sc->lis.p;
    for (int k = n - 1; k >= 0; k--) {
        if (!cross_point(sc->walls[seq[k]], img[k], cur, &hit[k])) return;
        cur = hit[k];
    }
    V2 from = sc->source;
    for (int k = 0; k < n; k++) {
        if (leg_blocked(sc, from, hit[k])) return;
        from = hit[k];
    }
    if (leg_blocked(sc, from, sc->lis.p)) return;
    if (sc->npaths >= MAX_PATHS) return;

    double plen = vlen(vsub(sc->lis.p, img[n - 1]));
    double amp = pow(ABSORP, n) / fmax(plen, 1.0);
    V2 arr = vnorm(vsub(hit[n - 1], sc->lis.p));   /* apparent source dir */
    Path* p = &sc->paths[sc->npaths++];
    p->az = atan2(vcross(sc->lis.fwd, arr), vdot(sc->lis.fwd, arr));
    p->energy = amp * amp;
    p->order = n;
}

static void gen(Scene* sc, int* seq, int depth) {
    if (depth > 0) trace(sc, seq, depth);
    if (depth == MAX_ORDER) return;
    for (int w = 0; w < sc->nwalls; w++) {
        if (depth > 0 && seq[depth - 1] == w) continue;  /* degenerate re-bounce */
        seq[depth] = w;
        gen(sc, seq, depth + 1);
    }
}

/* constant-power pan of sin(az) + far-ear shadow from interaural distance */
static void split_lr(double az, double e, double* el, double* er) {
    double s = sin(az);
    double x = -s;                              /* spatial.c: -1 left .. +1 right */
    double a = (x + 1.0) * 0.25 * PI;
    double gl = cos(a), gr = sin(a);
    double h = 1.0 / (1.0 + D_EAR * fabs(s) / S0);
    if (s > 0.0) gr *= h; else gl *= h;         /* az > 0: from left, right ear far */
    *el += e * gl * gl;
    *er += e * gr * gr;
}

static void stereo_sum(Scene* sc, double* el, double* er, double* az_out) {
    double l = 0.0, r = 0.0, sx = 0.0, sy = 0.0;
    for (int i = 0; i < sc->npaths; i++) {
        split_lr(sc->paths[i].az, sc->paths[i].energy, &l, &r);
        sx += sc->paths[i].energy * cos(sc->paths[i].az);
        sy += sc->paths[i].energy * sin(sc->paths[i].az);
    }
    *el = l;
    *er = r;
    *az_out = atan2(sy, sx);
}

static void run_scene(const char* name, Scene* sc) {
    int seq[MAX_ORDER];
    sc->npaths = 0;
    gen(sc, seq, 0);
    double el, er, az;
    stereo_sum(sc, &el, &er, &az);
    printf("RESULT %s paths=%d lr_ratio=%.3f azimuth_deg=%+.1f\n",
           name, sc->npaths, er > 0.0 ? el / er : 0.0, az * 180.0 / PI);
}

/* per-order azimuths for the corridor: where does direction get averaged out? */
static void corridor_diag(Scene* sc) {
    for (int o = 1; o <= MAX_ORDER; o++) {
        double el = 0.0, er = 0.0, sx = 0.0, sy = 0.0;
        printf("DIAG corridor order=%d:", o);
        for (int i = 0; i < sc->npaths; i++) {
            Path* p = &sc->paths[i];
            if (p->order != o) continue;
            printf(" az=%+.1f(E=%.2e)", p->az * 180.0 / PI, p->energy);
            split_lr(p->az, p->energy, &el, &er);
            sx += p->energy * cos(p->az);
            sy += p->energy * sin(p->az);
        }
        printf(" | order_lr=%.3f order_az=%+.1f\n",
               er > 0.0 ? el / er : 0.0, atan2(sy, sx) * 180.0 / PI);
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static Scene make_scene(V2 src, V2 lis) {
    Scene sc;
    sc.nwalls = 0;
    sc.npaths = 0;
    sc.source = src;
    sc.lis.p = lis;
    sc.lis.fwd = vnorm(vsub(src, lis));   /* face the direct sound */
    return sc;
}

int main(void) {
    /* geometry self-checks */
    Wall wx = { { 0.0, 0.0 }, { 10.0, 0.0 } };
    V2 m = reflect_point((V2){ 3.0, 4.0 }, wx);
    check("reflect_x", fabs(m.x - 3.0) < 1e-9 && fabs(m.y + 4.0) < 1e-9);
    V2 hit;
    check("cross_hit", cross_point(wx, (V2){ 5.0, 2.0 }, (V2){ 5.0, -2.0 }, &hit)
          && fabs(hit.x - 5.0) < 1e-9 && fabs(hit.y) < 1e-9);
    check("cross_miss", !cross_point(wx, (V2){ 12.0, 2.0 }, (V2){ 12.0, -2.0 }, &hit));

    /* symmetric corridor: source on centerline must give L/R == 1 */
    {
        Scene sc = make_scene((V2){ 5.0, 0.0 }, (V2){ 25.0, 0.0 });
        sc.walls[sc.nwalls++] = (Wall){ { 0.0, 1.5 }, { 40.0, 1.5 } };
        sc.walls[sc.nwalls++] = (Wall){ { 0.0, -1.5 }, { 40.0, -1.5 } };
        int seq[MAX_ORDER];
        gen(&sc, seq, 0);
        double el, er, az;
        stereo_sum(&sc, &el, &er, &az);
        check("symmetric_lr", fabs(el / er - 1.0) < 1e-6);
    }

    /* experiment 1: narrow corridor, 3 m wide, shot off-center */
    Scene corridor = make_scene((V2){ 5.0, 0.6 }, (V2){ 25.0, 0.0 });
    corridor.walls[corridor.nwalls++] = (Wall){ { 0.0, 1.5 }, { 40.0, 1.5 } };
    corridor.walls[corridor.nwalls++] = (Wall){ { 0.0, -1.5 }, { 40.0, -1.5 } };
    run_scene("corridor", &corridor);
    corridor_diag(&corridor);

    /* experiment 2: open field, walls 100 m away */
    Scene field = make_scene((V2){ -30.0, 5.0 }, (V2){ 0.0, 0.0 });
    field.walls[field.nwalls++] = (Wall){ { -100.0, -100.0 }, { 100.0, -100.0 } };
    field.walls[field.nwalls++] = (Wall){ { 100.0, -100.0 }, { 100.0, 100.0 } };
    field.walls[field.nwalls++] = (Wall){ { 100.0, 100.0 }, { -100.0, 100.0 } };
    field.walls[field.nwalls++] = (Wall){ { -100.0, 100.0 }, { -100.0, -100.0 } };
    run_scene("open_field", &field);

    /* experiment 3: corner, two walls meeting at 90 deg */
    Scene corner = make_scene((V2){ 4.0, 14.0 }, (V2){ 12.0, 9.0 });
    corner.walls[corner.nwalls++] = (Wall){ { 0.0, 0.0 }, { 30.0, 0.0 } };
    corner.walls[corner.nwalls++] = (Wall){ { 0.0, 0.0 }, { 0.0, 30.0 } };
    run_scene("corner", &corner);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
