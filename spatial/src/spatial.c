
/* spatial: constant-power panning and inverse-distance attenuation. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "stream.h"

typedef struct { float left, right; } Stereo;

static int pan(float x, Stereo* s) {
    if (x < -1.0f || x > 1.0f) return -1;
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

/*
 * Per-frame processing shared by both paths: one mono frame -> one stereo
 * frame, constant-power pan plus distance attenuation.
 */
static void spatial_frame(float in, float pan_x, float att, Stereo* out) {
    Stereo s = {0.0f, 0.0f};
    pan(pan_x, &s); /* pan_x is validated by the stream setup */
    out->left = in * s.left * att;
    out->right = in * s.right * att;
}

/*
 * Streaming spatializer. The spatial stage is length-neutral per pair; the
 * only length decision is how odd frame counts are handled:
 *  - legacy: each block rounds up to ceil(n/2) pairs, so a dangling last
 *    frame occupies a stereo slot by itself; the next block's first frame
 *    then lands in the wrong slot - the cumulative shift clicks at joins;
 *  - unified: the dangling frame is carried into the next block and paired
 *    with its first frame, so pair counts follow stream_total_pairs() and
 *    no extra stereo slot is ever invented.
 * The two L/R arrays are written separately; interleaved output is the
 * caller's job. Returns the number of pairs written.
 */
typedef struct {
    float pan_x;
    float att;
    int carry;          /* unified: pending frame flag */
    float held;         /* unified: pending frame value */
    int64_t pairs;
} SpatialStream;

MA_MAYBE_UNUSED static void spatial_stream_init(SpatialStream* st, float pan_x, float att) {
    st->pan_x = pan_x;
    st->att = att;
    st->carry = 0;
    st->held = 0.0f;
    st->pairs = 0;
}

/* Unified path: carry the odd frame across blocks, never emit padding.
 * A carried frame is paired (summed) with the next block's first frame,
 * preserving the whole-stream (frame 2p, frame 2p+1) pairing exactly. */
MA_MAYBE_UNUSED static int spatial_stream_block(SpatialStream* st,
                                                const float* in, int n,
                                                Stereo* out) {
    /* Pair count and the outgoing carry both come from stream.h, so this
     * stage never decides a length on its own. */
    int carry_in = st->carry;
    int pairs = (int)stream_block_pairs(n, &st->carry);

    int k = 0;          /* pairs written */
    int i = 0;          /* frames consumed from in[] */
    if (carry_in) {
        Stereo a, b;
        spatial_frame(st->held, st->pan_x, st->att, &a);
        spatial_frame(in[0], st->pan_x, st->att, &b);
        out[k].left = a.left + b.left;
        out[k].right = a.right + b.right;
        k++;
        i = 1;
    }
    for (; k < pairs; k++, i += 2) {
        Stereo a, b;
        spatial_frame(in[i], st->pan_x, st->att, &a);
        spatial_frame(in[i + 1], st->pan_x, st->att, &b);
        out[k].left = a.left + b.left;
        out[k].right = a.right + b.right;
    }
    if (st->carry)
        st->held = in[n - 1]; /* paired with the next block's first frame */
    st->pairs += pairs;
    return pairs;
}

/*
 * Legacy path: exactly ceil(n/2) pairs per block. A dangling last frame is
 * emitted as a full stereo frame on its own (no partner). Reconstructs the
 * pre-unification behavior byte for byte.
 */
MA_MAYBE_UNUSED static int spatial_legacy_block(const SpatialStream* st, const float* in,
                                int n, Stereo* out) {
    int pairs = (n + 1) / 2;
    for (int p = 0; p < pairs; p++) {
        spatial_frame(in[2 * p], st->pan_x, st->att, &out[p]);
        if (2 * p + 1 < n) {
            Stereo r;
            spatial_frame(in[2 * p + 1], st->pan_x, st->att, &r);
            out[p].left += r.left;
            out[p].right += r.right;
        }
    }
    return pairs;
}

#ifndef MA_NO_MAIN
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
#endif /* MA_NO_MAIN */

#ifndef MA_NO_MAIN
int main(void) {
    Stereo s;
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

    /* legacy blocks: ceil(n/2) pairs. With an odd block the dangling last
     * frame becomes a spurious stereo frame (hard-left pan: L keeps it,
     * R holds the same frame's pan result - not a zero pad). */
    float mono[5] = {0.1f, 0.2f, 0.3f, 0.4f, 0.5f};
    Stereo lb[4];
    SpatialStream ss0;
    spatial_stream_init(&ss0, -1.0f, 1.0f);
    int lp = spatial_legacy_block(&ss0, mono, 5, lb);
    check("legacy_pairs", lp == 3);
    check("legacy_spurious_pair",
          fabsf(lb[2].left - 0.5f) < 1e-6f &&
          fabsf(lb[2].right) < 1e-6f); /* sinf(-pi/2) ~ 0 */

    /* unified blocks: odd frame carries over and pairs with next block's
     * first frame; lengths match the stream formulas */
    Stereo ub[4];
    SpatialStream ss;
    spatial_stream_init(&ss, 0.0f, 1.0f);
    int k1 = spatial_stream_block(&ss, mono, 3, ub);
    int k2 = spatial_stream_block(&ss, mono + 3, 2, ub + k1);
    check("unified_carry", k1 == 1 && ss.carry == 1 && k2 == 1 &&
          (int64_t)(k1 + k2) == stream_total_pairs(5));

    /* unified carry must equal the whole-stream pair of the same frames */
    Stereo whole[4];
    SpatialStream sw;
    spatial_stream_init(&sw, 0.0f, 1.0f);
    int kw = spatial_stream_block(&sw, mono, 5, whole);
    check("unified_bitexact", kw == k1 + k2 &&
          memcmp(ub, whole, (size_t)(k1 + k2) * sizeof(Stereo)) == 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

#endif /* MA_NO_MAIN */
