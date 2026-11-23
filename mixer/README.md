# mixer

`make test` 跑自检；`make evidence` 跑证据批次（七份材料，缺一视为未完成）。

## 增益结构：按分贝域求和推总线余量

混音总线的问题两头都有：单路上来太轻听不见，多路一叠就爆。现在总线按
**分贝域**求和来估余量，再据此给每路分配：

| 量 | 公式 | 含义 |
|---|---|---|
| 相干和（最坏情况） | `S_peak = 20·log10( Σ 10^(L_i/20) )` | 所有路同时到达峰值时的和 |
| 非相干和（节目材料） | `S_rms = 10·log10( Σ 10^(L_i/10) )` | 各路不相关、功率相加 |
| 总线余量 | `H = target_dbfs − S_peak` | 负值表示超限，需要衰减 |
| 每路份额 | `share_i = 10^(L_i/20) / Σ 10^(L_j/20)` | 该路占相干和的比例 |
| 统一衰减 | `trim = H`（上抬有 `+24 dB` 上限） | 所有路共用，保持混音平衡 |
| 每路分配 | `alloc_i = L_i + trim` | 该路衰减后在总线里的峰值 |
| 每路余量 | `margin_i = 0 dBFS − alloc_i` | 离满刻度还剩多少 |

用**相干和**而不是功率和定余量，是为了保证任何相位关系下总线都不超过
`target_dbfs`（自检 `random_grouped_never_clips` 用 200 组随机相位验证）；
`rms_margin_db` 另行给出，用于看节目材料的实际余量有多少。

由于对所有路用**同一个** `trim`，比例关系不变（`balance_preserved`）：单路
源拿到全部预算（-18 dBFS 的单路被抬到 -1 dBFS，不再"不响"），N 路则分摊
（8 路 -6 dBFS 叠到 +12 dBFS 被统一压回 -1 dBFS，不再"爆"）。

```c
GainPlan plan;
gain_plan_init(&plan, peak_dbfs, nch, -1.0f);   /* 每路实测峰值电平 */
gain_bus_sum(ch, nch, frames, &plan, 0, out, &meter, spans, cap, &list);
/* 或一步到位：自动测峰值 -> 建计划 -> 求和 */
mix_bus_auto(ch, nch, frames, -1.0f, out, &plan, &meter, spans, cap, &list);
```

`gain_plan_init` 对非法输入返回 -1：路数越界、电平为 NaN/+inf、目标电平
高于 0 dBFS 或低于 -120 dBFS。静音路（`-inf`）不参与分配（`silent[i]=1`）。

## 关闭分组 = 回到现在

`MIX_GAIN_GROUP=0`（`make test NOGROUP=1`）时结构整体编译掉，
`gain_bus_sum(plan=NULL)` 退化为逐路增益 1.0 的直接求和：

- 自检程序 stdout 与**改动前逐字节一致**（`make test NOGROUP=1` 对照原始二进制）。
- API 级 `group_off_legacy_parity`：`gain_bus_sum(plan=NULL)` 与旧 `mix_ramp`
  在三种斜坡下的输出 `memcmp` 完全相同、削波计数相同。
- 证据里的 `[group_off_parity]` 逐斜坡复核，并打印 `group_off_byte_identical=1`。

`make evidence NOGROUP=1` 会拒绝运行，避免用被编译掉的实现出材料。

## 七份材料

`make evidence` 依次输出，`materials_complete=1` 仅在七项都非空且校验通过时成立：

| 材料 | 前缀 | 内容 |
|---|---|---|
| 总线余量 | `[bus_headroom]` | 相干/功率和、总线余量、目标、单路情形 |
| 各路余量分配 | `[channel_allocation]` | 混合电平原案：in/share/trim/alloc/margin 逐路列表 |
| 电平表 | `[level_meter]` | 统一前后总线的 peak/RMS/削波样本/区间数，及逐路电平 |
| 削波区间清单 | `[clip_spans]` | 半开区间 `[start,end)`、峰值、正负极性；覆盖样本数与电平表吻合 |
| 两档路数耗时 | `[timing]` | 2 路与 32 路，分组/关闭两配置的 ns/样本与比值 |
| 关闭一致性 | `[group_off_parity]` | 三种斜坡下与旧实现逐字节一致 |
| 覆盖清单 | `[coverage]` | 10 项检查（含单路可闻、多路下压、平衡保持、静音路处理） |

削波区间来自 `gain_bus_sum` 对**限幅前**和值的判定，区间总数 `ClipList.count`
与缓冲区不足标记 `overflow` 分开报告，缓冲区满时可扩容重取。
