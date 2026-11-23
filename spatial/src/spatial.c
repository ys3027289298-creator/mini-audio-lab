
/* spatial: constant-power panning and inverse-distance attenuation. */
#include <stdio.h>
#include <math.h>
#include <string.h>

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

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

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

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
