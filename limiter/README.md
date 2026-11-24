# limiter

`make test` 跑自检。

采样率假设：旧接口 `lim_init(..., sr)` 由调用方传采样率，保持不变。
统一入口 `lim_init_default(...)` 使用 `common/sr.h` 的 `MAL_SR_DEFAULT`（48000 Hz）。
对照数据见 `report/sr_report.txt`（`cd ../report && make test` 重新生成）。
