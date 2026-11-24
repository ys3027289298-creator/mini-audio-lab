
/* loudness: short-time loudness statistics robust to background noise.
 *
 * model: x(t) = a(t)*s(t) + n(t), speech-active fraction p, noise variance sn2.
 * window length T, noise bandwidth B -> window energy is chi-square with
 * nu = 2*B*T degrees of freedom, so the silence-cluster loudness spread is
 *     sigma_L = 4.34/sqrt(B*T)  dB.
 * edge (mixed speech+silence) windows make up w(T) ~ E*T/Ttot of all windows
 * (E = syllable edge rate), growing linearly in T; hence silence contamination
 * of any speech-side quantile decreases monotonically as T shortens, down to
 * the cluster-resolution limit
 *     T_min = (k*4.34/Delta)^2 / B,  Delta = 10*log10(1 + ss2/sn2),  k ~ 3.
 * the quantile must clear the silence+mixed mass: q > (1-p) + w(T)/2; we
 * back-solve q* = 1 - p_hat/2 (centre of the speech-active mass).
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define SR 16000
#define DUR_S 6
#define NSMP (SR * DUR_S)
#define BW (SR / 2)         /* white-noise bandwidth */
#define KSEP 3.0f           /* cluster-separation factor for T_min */
#define WIN_MS 30.0f        /* chosen window length */
#define DEV_THRESH_DB 1.0f  /* allowed quiet/noisy representative deviation */
#define PI 3.14159265358979f

static unsigned long long rng = 0x9e3779b97f4a7c15ULL;
static float frand(void) {
    rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
    return (float)((rng >> 40) & 0xFFFFFF) / 16777216.0f;
}
static float gnoise(void) { float s = 0.0f; for (int i = 0; i < 12; i++) s += frand(); return s - 6.0f; }

/* speech-like envelope: 1.2 s phrases every 1.9 s, 3 syllables/s at 55% duty */
static float envelope(float t, int* in_phrase) {
    float ph = fmodf(t, 1.9f);
    *in_phrase = ph < 1.2f;
    if (!*in_phrase) return 0.0f;
    float per = 1.0f / 3.0f, on = 0.55f * per, r = 0.010f;
    float sp = fmodf(ph, per);
    if (sp >= on) return 0.0f;
    float e = 1.0f;
    if (sp < r) e = sp / r;
    if (on - sp < r && (on - sp) / r < e) e = (on - sp) / r;
    return e;
}

/* same speaker, two rooms: quiet (-60 dBFS) and noisy (-30 dBFS) */
static void synth(float* x, char* act, char* phr, float sigma_n) {
    float hrms = 0.0f;
    for (int h = 1; h <= 5; h++) hrms += 1.0f / (2.0f * h * h);
    float amp = 0.1f / sqrtf(hrms); /* active RMS = 0.1 (-20 dBFS) */
    for (int i = 0; i < NSMP; i++) {
        float t = (float)i / SR;
        int ip;
        float e = envelope(t, &ip);
        float s = 0.0f;
        for (int h = 1; h <= 5; h++) s += sinf(2.0f * PI * 140.0f * h * t) / h;
        x[i] = e * amp * s + sigma_n * gnoise();
        act[i] = e > 0.5f ? 1 : 0;
        phr[i] = (char)ip;
    }
}

typedef struct {
    int nw;
    float mean_db, median_db, qstar, qstar_db, p_hat, noise_db;
    float gap_of_sil, gap_of_all; /* intra-phrase gaps read as silence */
    float contam;                 /* not-fully-active share of windows >= qstar_db */
} Stats;

static int cmpf(const void* a, const void* b) {
    float d = *(const float*)a - *(const float*)b;
    return d < 0 ? -1 : d > 0;
}

static float quant(const float* s, int n, float q) {
    float pos = q * (n - 1);
    int i = (int)pos;
    float f = pos - (float)i;
    return i + 1 < n ? s[i] * (1.0f - f) + s[i + 1] * f : s[n - 1];
}

static Stats analyze(const float* x, const char* act, const char* phr, int wlen) {
    Stats r = { 0 };
    int nw = NSMP / wlen;
    r.nw = nw;
    float* db = malloc((size_t)nw * sizeof *db);
    float* en = malloc((size_t)nw * sizeof *en);
    float* tmp = malloc((size_t)nw * sizeof *tmp);
    for (int w = 0; w < nw; w++) {
        float e = 0.0f;
        for (int i = 0; i < wlen; i++) { float v = x[w * wlen + i]; e += v * v; }
        en[w] = e / wlen;
        db[w] = 10.0f * log10f(en[w] + 1e-12f);
    }
    /* noise variance: mean energy of the lowest decile of windows */
    for (int w = 0; w < nw; w++) tmp[w] = en[w];
    qsort(tmp, (size_t)nw, sizeof *tmp, cmpf);
    float noise = 0.0f;
    for (int w = 0; w < nw / 10; w++) noise += tmp[w];
    noise /= nw / 10;
    r.noise_db = 10.0f * log10f(noise + 1e-12f);
    /* activity threshold: noise level plus k separation sigmas */
    float ethr = noise * (1.0f + KSEP / sqrtf((float)BW * wlen / SR));
    int act_cnt = 0;
    for (int w = 0; w < nw; w++) act_cnt += en[w] > ethr;
    r.p_hat = (float)act_cnt / nw;
    r.qstar = 1.0f - 0.5f * r.p_hat;
    if (r.qstar > 0.99f) r.qstar = 0.99f;
    if (r.qstar < 0.50f) r.qstar = 0.50f;
    for (int w = 0; w < nw; w++) tmp[w] = db[w];
    qsort(tmp, (size_t)nw, sizeof *tmp, cmpf);
    float mean = 0.0f;
    for (int w = 0; w < nw; w++) mean += db[w];
    r.mean_db = mean / nw;
    r.median_db = quant(tmp, nw, 0.5f);
    r.qstar_db = quant(tmp, nw, r.qstar);
    /* gap diagnostics: estimator-silent windows whose samples sit in phrases */
    int sil = 0, gap = 0;
    for (int w = 0; w < nw; w++) {
        if (en[w] > ethr) continue;
        sil++;
        int ip = 0;
        for (int i = 0; i < wlen; i++) ip += phr[w * wlen + i];
        gap += ip > wlen / 2;
    }
    r.gap_of_sil = sil ? (float)gap / sil : 0.0f;
    r.gap_of_all = (float)gap / nw;
    /* contamination: windows at/above qstar_db that are not fully active */
    int above = 0, dirty = 0;
    for (int w = 0; w < nw; w++) {
        if (db[w] < r.qstar_db) continue;
        above++;
        int a = 0;
        for (int i = 0; i < wlen; i++) a += act[w * wlen + i];
        dirty += a < wlen;
    }
    r.contam = above ? (float)dirty / above : 0.0f;
    free(db); free(en); free(tmp);
    return r;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    static float xq[NSMP], xn[NSMP];
    static char act[NSMP], phr[NSMP];
    synth(xq, act, phr, 1.0e-3f);   /* quiet room */
    synth(xn, act, phr, 3.1623e-2f); /* noisy room, same speech */

    /* T_min from the noisy-room cluster separation Delta */
    float delta = 10.0f * log10f(1.0f + 1.0e-2f / 1.0e-3f);
    float tmin_ms = KSEP * 4.34f / delta;
    tmin_ms = tmin_ms * tmin_ms / BW * 1.0e3f;
    printf("delta_db=%.1f tmin_ms=%.2f win_ms=%.0f\n", delta, tmin_ms, WIN_MS);

    int wlen = (int)(SR * WIN_MS / 1000.0f);
    Stats q = analyze(xq, act, phr, wlen);
    Stats n = analyze(xn, act, phr, wlen);
    printf("quiet mean=%.2f median=%.2f qstar=%.3f qstar_db=%.2f p_hat=%.3f noise_db=%.1f\n",
           q.mean_db, q.median_db, q.qstar, q.qstar_db, q.p_hat, q.noise_db);
    printf("noisy mean=%.2f median=%.2f qstar=%.3f qstar_db=%.2f p_hat=%.3f noise_db=%.1f\n",
           n.mean_db, n.median_db, n.qstar, n.qstar_db, n.p_hat, n.noise_db);

    float dev_mean = fabsf(q.mean_db - n.mean_db);
    float dev_med = fabsf(q.median_db - n.median_db);
    float dev_qs = fabsf(q.qstar_db - n.qstar_db);
    printf("dev_mean=%.2f dev_median=%.2f dev_qstar=%.2f thresh=%.2f\n",
           dev_mean, dev_med, dev_qs, DEV_THRESH_DB);

    /* window sweep on the noisy clip: silence contamination of the q*
     * statistic must decrease monotonically as T shortens toward T_min */
    static const int sweep_ms[] = { 160, 80, 40, 20, 10 };
    int mono = 1;
    float prev = 1.0f;
    for (int i = 0; i < 5; i++) {
        int wl = SR * sweep_ms[i] / 1000;
        Stats s = analyze(xn, act, phr, wl);
        float sig_l = 4.34f / sqrtf((float)BW * wl / SR);
        printf("sweep win_ms=%d contam=%.4f sigma_l=%.2f\n", sweep_ms[i], s.contam, sig_l);
        if (s.contam > prev + 1e-6f) mono = 0;
        prev = s.contam;
    }

    if (dev_qs > DEV_THRESH_DB) {
        /* noise variance estimate mistook speech gaps for silence */
        printf("gap_share quiet: of_sil=%.3f of_all=%.3f\n", q.gap_of_sil, q.gap_of_all);
        printf("gap_share noisy: of_sil=%.3f of_all=%.3f\n", n.gap_of_sil, n.gap_of_all);
    }

    check("qstar_deviation_ok", dev_qs <= DEV_THRESH_DB);
    check("mean_deviation_large", dev_mean > 5.0f);
    check("monotonic_contamination", mono);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
