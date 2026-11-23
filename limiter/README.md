# limiter

峰值限幅：一阶包络按启动/释放时间常数逼近逐样本目标增益
（`ceiling / |x|`），超限样本乘包络后再硬限幅到 ceiling。

`make test` 跑自检，并打印四项验收数据：

- `[envelope_curve]` — 阶跃信号（0.2 → 2.0 → 0.2）下逐点包络曲线，附带闭式解和误差
- `[gain_reduction_table]` — 20 ms 窗口的峰值输入/输出、最小包络、最大增益衰减 dB
- `[time_constants]` — 5 ms 启动与 50 ms 释放的实测 vs 解析 t63 / t90 / t99 耗时
- `[bypass_identity]` — `bypass=1` 时输出与输入逐字节相同（memcmp + FNV-1a）

包络递推为 `env[n+1] = c*env[n] + (1-c)*target`，`c = exp(-1/(tau*sr))`；
恒定目标下的闭式为 `env[n] = target + (env0-target)*exp(-n/(tau*sr))`。
