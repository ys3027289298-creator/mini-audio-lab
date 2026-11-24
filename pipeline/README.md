# pipeline

流式演示管线（resampler → mixer → spatial）与耗时基准，一份源码构建两个二进制：

| 二进制 | 宏 | 行为 |
|---|---|---|
| `pipeline_legacy` | `MA_UNIFY=0` | 改造前的逐块长度推导，逐字节复现旧行为 |
| `pipeline_unified` | `MA_UNIFY=1` | 长度全部来自 `stream/stream.h` 的统一公式 |

```
make test                      # 构建两种模式并跑演示
./pipeline_unified demo <dir> # 写各环节 raw float32 与 run.json
./pipeline_unified bench       # 采样率翻倍耗时基准（CSV）
```

旧函数签名一字未改；新流式函数带 `MA_MAYBE_UNUSED`，在 legacy 构建中不产生警告。
演示内置整流参考实现，两种模式都对三个环节做字节级自检（`all_byte_exact`）。
