# dither

TPDF 抖动 + 16 位量化，附一份可复跑的证据包。

```bash
make test                     # 21 项 C 自检
./evidence/run_all.sh         # 导出 WAV、独立校验、生成图表与指标 JSON
```

- `REPORT.md` —— 结论文档：信噪比、噪声频谱、与无抖动的对照、逐字节复现性。
- `evidence/verify_dither.py` —— 只读导出的 WAV 独立复测，不链接被测代码。
- `DITHER=off ./dither_export tone60 out.wav` —— 关闭抖动，输出与纯舍入逐字节一致。
