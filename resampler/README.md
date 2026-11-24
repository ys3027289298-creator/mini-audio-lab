# resampler

`make test` 跑自检。

采样率假设：旧接口 `resample_linear(in, n, ratio, out, cap)` 只吃比率、无绝对采样率，保持不变。
统一入口 `resample_rates(in, n, sr_in, sr_out, out, cap)` 用 `common/sr.h` 的 `mal_ratio()` 推导比率。
对照数据见 `report/sr_report.txt`（`cd ../report && make test` 重新生成）。
