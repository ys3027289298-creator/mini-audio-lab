/* meter: 导出/播放共用的统一电平表头。
 *
 * 统一口径(两条路径必须完全一致):
 *   测量点 = 链路输出(post-limiter 样本流)
 *   窗口   = 整段信号(不分块、不截尾块)
 *   单位   = dBFS, 参考 1.0
 *   指标   = peak 与 RMS 同时报告, 不允许一边报峰值一边报 RMS
 * 表头只读样本, 不修改音频; 开关统一不影响输出样本。 */
#ifndef METER_H
#define METER_H

#include <math.h>

typedef struct {
    double sum_sq;  /* 平方和, 双精度累加 */
    float peak;     /* 绝对值峰值 */
    long count;     /* 样本数 */
    int clipped;    /* |x| >= 1.0 的样本数 */
} Meter;

static inline void meter_reset(Meter* m) {
    m->sum_sq = 0.0;
    m->peak = 0.0f;
    m->count = 0;
    m->clipped = 0;
}

static inline void meter_feed(Meter* m, const float* x, int n) {
    for (int i = 0; i < n; i++) {
        float a = fabsf(x[i]);
        if (a > m->peak) m->peak = a;
        if (a >= 1.0f) m->clipped++;
        m->sum_sq += (double)x[i] * (double)x[i];
        m->count++;
    }
}

static inline float meter_peak_dbfs(const Meter* m) {
    return m->peak > 0.0f ? 20.0f * log10f(m->peak) : -INFINITY;
}

static inline float meter_rms_dbfs(const Meter* m) {
    if (m->count <= 0 || m->sum_sq <= 0.0) return -INFINITY;
    return 20.0f * log10f((float)sqrt(m->sum_sq / (double)m->count));
}

#endif
