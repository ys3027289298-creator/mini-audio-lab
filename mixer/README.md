# mixer

`make test` 跑自检。

## 削波定位

削波不再只报“发生过”：`mix_ramp_locate` 把削波采样聚成极大连续区间，半开区间
`[start, end)`；`MixStat` 同时给出波形摘要（peak/rms/dc/crest、正负极性计数、
首末削波采样）。

```c
ClipSpan spans[8];
MixStat st;
mix_ramp_locate(a, b, 1.0f, 1.0f, out, n, ramp_frames, &st, spans, 8);
/* st.span_count 为真实区间总数；spans 缓冲不足时 st.span_overflow 置 1，
   st.span_written 为实际写出数量，调用者据此扩容重取。 */
```

`mix_ramp` 保留原签名，内部走同一路径但不写区间（`span_written == 0`）。

## 关闭定位

`make test NOLOCATE=1` 时以 `-DMIX_CLIP_LOCATE=0` 构建，回到原始“只报次数”
实现，程序输出与改动前逐字节一致。证据跑批 `make evidence` 在该配置下不可用。

## 证据

`make evidence` 运行注入缺陷跑批，输出四段材料：`[clip_spans]` 削波区间清单、
`[wave_summary]` 波形摘要、`[injection_capture]` 注入后捕获率（含
`control_no_location_rate`：仅报次数的检测器捕获率为 0；`control_shift_rejected`：
区间偏移 1 采样会被判错）、`[coverage]` 覆盖清单。
