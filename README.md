# mini-audio-lab

游戏音频基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `resampler/` | 线性重采样、输出长度推导 |
| `mixer/` | 增益斜坡混音、削波检测 |
| `fft/` | 基 2 快速傅里叶变换、幅度、Parseval 校验 |
| `limiter/` | 峰值限幅、启动释放包络、增益衰减报告 |
| `spatial/` | 等功率声像、距离衰减 |
| `common/` | 统一采样率假设（`sr.h`：默认 48000、备选 44100，ms↔帧↔比率换算） |
| `report/` | 交叉验证报告：长度对照表、各环波形摘要、两种采样率耗时、旧接口兼容清单 |

采样率假设已统一到 `common/sr.h`：resampler（比率）、mixer（斜坡帧数）、limiter（ms+采样率）
三环此前各自为政，现都经同一组换算函数推导。`cd report && make test` 生成对照数据，
`make check-compat` 校验五个模块自测输出与统一前逐字节相同（golden 文件在 `report/golden/`）。
