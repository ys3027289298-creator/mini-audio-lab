# meter

导出文件与实时播放的电平统一口径。

- `src/meter.h` — 两侧共用的统一表头(测量点 post-limiter、整段窗口、dBFS、peak+RMS 同报)。
- `src/lab.c` — 对照台架: 导出(一次性渲染) vs 播放(128 帧回调), 两种信号长度(1024 块对齐 / 1000 短尾块)。

`make test` 跑自检并验证逐字节一致; `make report` 生成四样材料到 `out/`:

1. `out/waveform_summary.md` — 两边波形摘要
2. `out/level_diff.md` — 电平差表(统一前/后)
3. `out/length_compat.md` — 两种信号长度对照与兼容清单
4. `out/byte_identity.md` — 关闭统一(`--no-unify`)后两处逐字节一致
