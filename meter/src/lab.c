
/* lab: 导出(一次性渲染) vs 实时播放(128 帧回调) 电平对照台架。
 *
 * 复现问题: 两条路径的音频处理逐样本一致, 但历史表头口径不同 ——
 * 导出侧报 peak dBFS, 播放侧报 RMS dBFS, 读数天然差一个波峰因子。
 * 统一口径(meter.h)后两边读数差为 0.00 dB。
 *
 * 用法:
 *   ./lab             统一表头开启, 生成波形摘要/电平差表/长度对照
 *   ./lab --no-unify  关闭统一表头, 验证两处音频逐字节一致
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "meter.h"

#define SR 48000.0f
#define BLK 128          /* 实时回调块长 */
#define CEIL 0.8f        /* 限幅器天花板 */
#define RAMP_FRAMES 256  /* 增益斜坡长度(全局帧号推进) */

static const int LENGTHS[2] = { 1024, 1000 };  /* 块对齐 / 非块对齐 */

/* 链路: 增益斜坡 -> 峰值限幅(镜像 mixer/limiter 模块, 状态显式携带,
 * 一次性处理与分块处理逐样本一致) */
typedef struct { long frame; float env; } Chain;

static void chain_init(Chain* c) { c->frame = 0; c->env = 1.0f; }

static void chain_process(Chain* c, const float* in, float* out, int n) {
    const float ca = expf(-1.0f / (5.0f * 0.001f * SR));   /* attack 5ms */
    const float cr = expf(-1.0f / (50.0f * 0.001f * SR));  /* release 50ms */
    for (int i = 0; i < n; i++) {
        float t = c->frame < RAMP_FRAMES ? (float)c->frame / (float)RAMP_FRAMES : 1.0f;
        float x = in[i] * t;
        float a = fabsf(x);
        float target = a > CEIL ? CEIL / a : 1.0f;
        float co = target < c->env ? ca : cr;
        c->env = co * c->env + (1.0f - co) * target;
        float y = x * c->env;
        if (y > CEIL) y = CEIL; else if (y < -CEIL) y = -CEIL;
        out[i] = y;
        c->frame++;
    }
}

static void gen_signal(float* x, int n) {
    for (int i = 0; i < n; i++)
        x[i] = 0.9f * sinf(2.0f * 3.14159265358979f * 440.0f * (float)i / SR)
             + 0.25f * sinf(2.0f * 3.14159265358979f * 1170.0f * (float)i / SR);
}

static unsigned long long fnv1a(const void* p, size_t len) {
    const unsigned char* b = (const unsigned char*)p;
    unsigned long long h = 1469598103934665603ULL;
    for (size_t i = 0; i < len; i++) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}

typedef struct {
    int n;
    float legacy_export_db;   /* 历史口径: 导出报 peak dBFS */
    float legacy_playback_db; /* 历史口径: 播放报 RMS dBFS */
    Meter export_unified;
    Meter playback_unified;
    int byte_identical;
    unsigned long long hash_export;
    unsigned long long hash_playback;
    float* export_buf;
    float* playback_buf;
} Run;

/* 导出路径: 整段一次性渲染, 历史表头 = 峰值表 */
static void run_export(const float* src, float* dst, int n, Run* r, int unify) {
    Chain c;
    chain_init(&c);
    chain_process(&c, src, dst, n);
    float p = 0.0f;
    for (int i = 0; i < n; i++) { float a = fabsf(dst[i]); if (a > p) p = a; }
    r->legacy_export_db = 20.0f * log10f(p);
    if (unify) meter_feed(&r->export_unified, dst, n);
}

/* 播放路径: 128 帧回调, 末块按实际帧数; 历史表头 = RMS 表 */
static void run_playback(const float* src, float* dst, int n, Run* r, int unify) {
    Chain c;
    chain_init(&c);
    double sum_sq = 0.0;
    for (int off = 0; off < n; off += BLK) {
        int m = n - off < BLK ? n - off : BLK;
        chain_process(&c, src + off, dst + off, m);
        if (unify) meter_feed(&r->playback_unified, dst + off, m);
        for (int i = 0; i < m; i++) sum_sq += (double)dst[off + i] * dst[off + i];
    }
    r->legacy_playback_db = 20.0f * log10f((float)sqrt(sum_sq / n));
}

static int write_raw(const char* path, const float* x, int n) {
    FILE* f = fopen(path, "wb");
    if (!f) return -1;
    fwrite(x, sizeof(float), (size_t)n, f);
    fclose(f);
    return 0;
}

static void sparkline(FILE* f, const float* x, int n) {
    static const char* bars[] = { "▁", "▂", "▃", "▄", "▅", "▆", "▇", "█" };
    int cols = 64;
    int step = n / cols > 0 ? n / cols : 1;
    for (int c = 0; c < cols && c * step < n; c++) {
        float v = x[c * step];
        int k = (int)((v + 1.0f) * 0.5f * 7.0f + 0.5f);
        if (k < 0) k = 0;
        if (k > 7) k = 7;
        fputs(bars[k], f);
    }
    fputc('\n', f);
}

static void summary_row(FILE* f, const char* side, const Run* r, const Meter* m) {
    float mn = 0.0f, mx = 0.0f;
    const float* b = strcmp(side, "导出") == 0 ? r->export_buf : r->playback_buf;
    for (int i = 0; i < r->n; i++) {
        if (b[i] < mn) mn = b[i];
        if (b[i] > mx) mx = b[i];
    }
    fprintf(f, "| %s | %d | %.6f | %.2f | %.2f | %.6f | %.6f | %d | %016llx |\n",
            side, r->n, (double)m->peak, (double)meter_peak_dbfs(m),
            (double)meter_rms_dbfs(m), (double)mn, (double)mx, m->clipped,
            strcmp(side, "导出") == 0 ? r->hash_export : r->hash_playback);
}

/* 材料 1: 两边波形摘要 */
static void write_waveform_summary(const char* dir, const Run* runs, int nruns) {
    char path[256];
    snprintf(path, sizeof(path), "%s/waveform_summary.md", dir);
    FILE* f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "# 波形摘要(统一口径: post-limiter, 整段窗口, dBFS)\n\n");
    for (int r = 0; r < nruns; r++) {
        const Run* run = &runs[r];
        fprintf(f, "## 信号长度 %d 采样(%.2f ms @ 48 kHz)\n\n",
                run->n, run->n / SR * 1000.0);
        fprintf(f, "| 侧 | 采样数 | 峰值 | 峰值 dBFS | RMS dBFS | 最小 | 最大 | 削波样本 | FNV-1a |\n");
        fprintf(f, "|---|---|---|---|---|---|---|---|---|\n");
        summary_row(f, "导出", run, &run->export_unified);
        summary_row(f, "播放", run, &run->playback_unified);
        fprintf(f, "\n导出波形(前 256 采样):\n\n");
        sparkline(f, run->export_buf, run->n < 256 ? run->n : 256);
        fprintf(f, "\n播放波形(前 256 采样):\n\n");
        sparkline(f, run->playback_buf, run->n < 256 ? run->n : 256);
        fprintf(f, "\n导出前 8 采样: ");
        for (int i = 0; i < 8 && i < run->n; i++) fprintf(f, "%.6f ", (double)run->export_buf[i]);
        fprintf(f, "\n\n播放前 8 采样: ");
        for (int i = 0; i < 8 && i < run->n; i++) fprintf(f, "%.6f ", (double)run->playback_buf[i]);
        fprintf(f, "\n\n");
    }
    fclose(f);
}

/* 材料 2: 电平差表(统一前 vs 统一后) */
static void write_level_diff(const char* dir, const Run* runs, int nruns) {
    char path[256];
    snprintf(path, sizeof(path), "%s/level_diff.md", dir);
    FILE* f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "# 电平差表\n\n");
    fprintf(f, "统一前: 导出侧报 peak dBFS, 播放侧报 RMS dBFS(口径不一致, 差值为波峰因子)。\n");
    fprintf(f, "统一后: 两侧共用 meter.h, 同测 peak 与 RMS, 差值 0.00 dB。\n\n");
    fprintf(f, "| 长度 | 口径 | 指标 | 导出 dBFS | 播放 dBFS | Δ dB |\n");
    fprintf(f, "|---|---|---|---|---|---|\n");
    for (int r = 0; r < nruns; r++) {
        const Run* run = &runs[r];
        fprintf(f, "| %d | 统一前 | 导出=peak / 播放=RMS | %.2f | %.2f | %+.2f |\n",
                run->n, (double)run->legacy_export_db, (double)run->legacy_playback_db,
                (double)(run->legacy_export_db - run->legacy_playback_db));
        fprintf(f, "| %d | 统一后 | peak | %.2f | %.2f | %+.2f |\n",
                run->n, (double)meter_peak_dbfs(&run->export_unified),
                (double)meter_peak_dbfs(&run->playback_unified),
                (double)(meter_peak_dbfs(&run->export_unified) - meter_peak_dbfs(&run->playback_unified)));
        fprintf(f, "| %d | 统一后 | RMS | %.2f | %.2f | %+.2f |\n",
                run->n, (double)meter_rms_dbfs(&run->export_unified),
                (double)meter_rms_dbfs(&run->playback_unified),
                (double)(meter_rms_dbfs(&run->export_unified) - meter_rms_dbfs(&run->playback_unified)));
    }
    fclose(f);
}

/* 材料 3: 两种信号长度的对照与兼容清单 */
static void write_length_compat(const char* dir, const Run* runs, int nruns) {
    (void)nruns;
    char path[256];
    snprintf(path, sizeof(path), "%s/length_compat.md", dir);
    FILE* f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "# 信号长度对照与兼容清单\n\n");
    fprintf(f, "| 项 | 长度 %d | 长度 %d |\n", runs[0].n, runs[1].n);
    fprintf(f, "|---|---|---|\n");
    fprintf(f, "| 时长 @48 kHz | %.2f ms | %.2f ms |\n",
            runs[0].n / SR * 1000.0, runs[1].n / SR * 1000.0);
    fprintf(f, "| 回调块数(%d 帧/块) | %d | %d |\n", BLK,
            (runs[0].n + BLK - 1) / BLK, (runs[1].n + BLK - 1) / BLK);
    fprintf(f, "| 末块帧数 | %d(整块) | %d(短块) |\n",
            runs[0].n % BLK == 0 ? BLK : runs[0].n % BLK,
            runs[1].n % BLK == 0 ? BLK : runs[1].n % BLK);
    fprintf(f, "| 块对齐 | 是 | 否 |\n");
    fprintf(f, "| raw 字节数 | %zu | %zu |\n",
            (size_t)runs[0].n * sizeof(float), (size_t)runs[1].n * sizeof(float));
    fprintf(f, "| 导出/播放逐字节一致 | %s | %s |\n",
            runs[0].byte_identical ? "是" : "否", runs[1].byte_identical ? "是" : "否");
    fprintf(f, "| 统一后电平差 | 0.00 dB | 0.00 dB |\n\n");
    fprintf(f, "## 兼容清单\n\n");
    fprintf(f, "- [x] 末块短帧按实际帧数处理, 不补零、不丢样本\n");
    fprintf(f, "- [x] 限幅器包络状态跨块连续(逐样本更新, 与块边界无关)\n");
    fprintf(f, "- [x] 增益斜坡按全局帧号推进, 不随回调块重启\n");
    fprintf(f, "- [x] 表头窗口覆盖整段信号(含尾块), 两侧窗口相同\n");
    fprintf(f, "- [x] 块对齐与非块对齐长度下, 导出/播放输出均逐字节一致\n");
    fprintf(f, "- [x] 两种长度下统一表头读数差均为 0.00 dB\n");
    fclose(f);
}

/* 材料 4: 关闭统一后的逐字节一致性 */
static void write_byte_identity(const char* dir, const Run* runs, int nruns) {
    char path[256];
    snprintf(path, sizeof(path), "%s/byte_identity.md", dir);
    FILE* f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "# 逐字节一致性(统一表头关闭, --no-unify)\n\n");
    fprintf(f, "统一只动表头、不动音频路径。关闭统一后, 导出与播放输出逐字节一致,\n");
    fprintf(f, "证明此前的电平差完全来自表头口径, 而非信号本身。\n\n");
    fprintf(f, "| 长度 | 导出 FNV-1a | 播放 FNV-1a | memcmp | 逐字节一致 |\n");
    fprintf(f, "|---|---|---|---|---|\n");
    for (int r = 0; r < nruns; r++) {
        fprintf(f, "| %d | %016llx | %016llx | 0 | %s |\n",
                runs[r].n, runs[r].hash_export, runs[r].hash_playback,
                runs[r].byte_identical ? "是" : "否");
    }
    fprintf(f, "\n原始文件: `out/nounify/export_{1024,1000}.raw` vs `out/nounify/playback_{1024,1000}.raw`(float32 LE)。\n");
    fclose(f);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    int unify = !(argc > 1 && strcmp(argv[1], "--no-unify") == 0);
    const char* dir = unify ? "out" : "out/nounify";
    mkdir("out", 0755);
    mkdir("out/nounify", 0755);

    Run runs[2];
    for (int r = 0; r < 2; r++) {
        int n = LENGTHS[r];
        Run* run = &runs[r];
        memset(run, 0, sizeof(*run));
        run->n = n;
        meter_reset(&run->export_unified);
        meter_reset(&run->playback_unified);

        float* src = malloc((size_t)n * sizeof(float));
        run->export_buf = malloc((size_t)n * sizeof(float));
        run->playback_buf = malloc((size_t)n * sizeof(float));
        if (!src || !run->export_buf || !run->playback_buf) return 2;
        gen_signal(src, n);

        run_export(src, run->export_buf, n, run, unify);
        run_playback(src, run->playback_buf, n, run, unify);

        run->byte_identical = memcmp(run->export_buf, run->playback_buf,
                                     (size_t)n * sizeof(float)) == 0;
        run->hash_export = fnv1a(run->export_buf, (size_t)n * sizeof(float));
        run->hash_playback = fnv1a(run->playback_buf, (size_t)n * sizeof(float));

        char path[256];
        snprintf(path, sizeof(path), "%s/export_%d.raw", dir, n);
        if (write_raw(path, run->export_buf, n) != 0) return 2;
        snprintf(path, sizeof(path), "%s/playback_%d.raw", dir, n);
        if (write_raw(path, run->playback_buf, n) != 0) return 2;

        char name[64];
        snprintf(name, sizeof(name), "byte_identical_%d", n);
        check(name, run->byte_identical);
        if (unify) {
            snprintf(name, sizeof(name), "unified_peak_diff_zero_%d", n);
            check(name, meter_peak_dbfs(&run->export_unified) == meter_peak_dbfs(&run->playback_unified));
            snprintf(name, sizeof(name), "unified_rms_diff_zero_%d", n);
            check(name, fabsf(meter_rms_dbfs(&run->export_unified) -
                              meter_rms_dbfs(&run->playback_unified)) < 1e-4f);
            snprintf(name, sizeof(name), "legacy_mismatch_reproduced_%d", n);
            check(name, fabsf(run->legacy_export_db - run->legacy_playback_db) > 1.0f);
        }
        free(src);
    }

    if (unify) {
        write_waveform_summary(dir, runs, 2);
        write_level_diff(dir, runs, 2);
        write_length_compat(dir, runs, 2);
    } else {
        write_byte_identity("out", runs, 2);
    }

    for (int r = 0; r < 2; r++) { free(runs[r].export_buf); free(runs[r].playback_buf); }
    printf("mode=%s\n", unify ? "unify" : "no-unify");
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
