# resampler

`make test` 跑自检。

- `src/resampler.c`：线性重采样与输出长度推导。
- `src/aliasing_bound.c`：由滤波器阶数推混叠残留解析上界、反解满足目标残留的最小倍率、三档实测验证。推导与结果见 [ALIASING.md](ALIASING.md)。
