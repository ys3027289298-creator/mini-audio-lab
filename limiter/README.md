# limiter

`make test` 跑自检。

`make compare` 跑新旧释放策略对照：同一段素材（0.4 s 超标大声段 + 1.6 s
小声段）上叠加旧（单一 300 ms 释放）与新（按电平分包络：大声 60 ms 快放、
小声 400 ms 慢放，分裂点 0.25）两条增益曲线，另跑一段纯脉冲查增益过冲。
图在 `out/gain_compare.png`，频段时间对照表在 `out/segments.md`。
