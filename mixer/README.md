# mixer

`make test` 跑自检。

采样率假设：旧接口 `mix_ramp(..., ramp_frames, ...)` 以帧数表示斜坡，保持不变。
统一入口 `mix_ramp_ms(..., ramp_ms, sr, ...)` 用 `common/sr.h` 的 `mal_ms_to_frames()` 把毫秒换算成帧。
对照数据见 `report/sr_report.txt`（`cd ../report && make test` 重新生成）。
