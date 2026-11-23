# spatial

等功率声像与**统一距离衰减**。

`make test` 跑自检（含证据生成）；`make report` 只生成材料到 `out/`。

## 统一衰减

旧实现按远近两套衰减参数在硬边界切换，切换处电平突变。现在只有一组参数，
整条曲线 C1（值、一阶导都连续）：

```
d <= ref              : g = 1
ref < d <= ref + knee : g = ref / (ref + rolloff * h(d - ref))
d > ref + knee        : g = ref / (ref + rolloff * (d - ref))

h(x) = x^2 * (2*knee - x) / knee^2,  0 <= x <= knee
```

- `h(0) = 0, h'(0) = 0`：近场平台与膝部在 `d = ref` 处值、斜率同时相等。
- `h(knee) = knee, h'(knee) = 1`：膝部与远场在 `d = ref + knee` 处值、斜率同时相等。
- `h' >= 0` 于 `[0, knee]`，因此 `g` 全域单调不增。
- `rolloff = 1` 时远场就是 `ref / d`，与画面中目标的投影尺寸同一规律。

默认配置 `ref = 1 m, rolloff = 1, knee = 0.5 m`。

## 错误

`dist == 0`（含 `-0.0`）返回 `ATTEN_ERR_DIST_ZERO`，`dist < 0` 返回
`ATTEN_ERR_DIST_NEG`，NaN/inf 返回 `ATTEN_ERR_NONFINITE`；参数非法分别返回
`ATTEN_ERR_REF / ROLLOFF / KNEE`。所有错误路径把输出增益写为 `0.0`，**不会**
静默返回满电平。

## 证据材料（`make report` → `out/`）

| 文件 | 内容 |
|---|---|
| `distance_level.csv` | 距离-电平：旧模型/统一模型的增益与 dB，及画面预期 `min(1, ref/d)` |
| `first_diff.csv` | 一阶差分：相邻距离的增益差与 dB 差，含接缝与旧切换点标注 |
| `jump_counts.csv` | 突跳采样计数：音频采样/视频帧/距离扫描三条序列，多阈值下超限样本数 |
| `event_deviation.csv` | 与画面事件对齐后的偏差：逐帧、事件标记、与画面线索的 dB 偏差 |
| `summary.md` | C1 证明、接缝实测值、关键数字汇总 |
| `error_cases.csv` | 零距离等错误路径的返回码与输出 |

`spatial.c` 的自检会同时检查定义域、单调性、接缝单侧极限，以及零距离不被当作满电平。
