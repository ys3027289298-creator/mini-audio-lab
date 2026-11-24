#ifndef MAL_SR_H
#define MAL_SR_H

/* sr: unified sample-rate assumptions for mini-audio-lab.
 *
 * Before this header existed, the three time-aware stages each baked in
 * their own assumption:
 *   - resampler took a bare ratio with no notion of absolute rate
 *   - mixer took ramp length in frames (time meaning depends on rate)
 *   - limiter took ms plus a caller-supplied rate (tests hardcoded 48000)
 * Cross-validating one stage against another was impossible because the
 * same physical duration produced different frame counts.
 *
 * Unified assumption: the pipeline runs at MAL_SR_DEFAULT (48000 Hz);
 * MAL_SR_ALT (44100 Hz) is the second supported rate for comparison.
 * All durations are expressed in milliseconds and converted to frames
 * through mal_ms_to_frames(); all rate conversions go through
 * mal_ratio(). Rounding is round-half-away-from-zero everywhere.
 */

#include <math.h>

#define MAL_SR_DEFAULT 48000.0f
#define MAL_SR_ALT     44100.0f

/* ms -> frames at rate sr, rounded half away from zero. Returns -1 on
 * non-positive ms or sr so callers can propagate their own -1. */
static inline int mal_ms_to_frames(float ms, float sr) {
    if (ms <= 0.0f || sr <= 0.0f) return -1;
    double f = (double)ms * 0.001 * (double)sr;
    return (int)llround(f);
}

/* frames -> ms at rate sr. */
static inline float mal_frames_to_ms(int frames, float sr) {
    if (frames < 0 || sr <= 0.0f) return -1.0f;
    return (float)((double)frames / (double)sr * 1000.0);
}

/* resample ratio (source frames per dest frame) for sr_in -> sr_out. */
static inline double mal_ratio(float sr_in, float sr_out) {
    if (sr_in <= 0.0f || sr_out <= 0.0f) return -1.0;
    return (double)sr_in / (double)sr_out;
}

#endif /* MAL_SR_H */
