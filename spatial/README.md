# spatial

`make test` 跑自检。

`pan(x)` 端点严格归零：`x = -1` 输出 `L = 1.0f, R = +0.0f`，`x = +1` 输出 `L = +0.0f, R = 1.0f`，静音侧是正零而不是 `sinf`/`cosf` 的浮点残差。
内部 `|x| < 1` 仍走等功率曲线，逐位不变。`pan_set_exact_endpoints(0)` 可整体退回旧公式，输出与原实现逐字节相同。
