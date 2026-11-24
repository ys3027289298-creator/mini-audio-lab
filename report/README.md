# report

交叉验证报告：把 resampler/mixer/limiter 的真实 DSP 代码链接进一个程序，
在统一采样率假设（`common/sr.h`）下输出四组对照数据。

- `make test` — 生成报告（最新结果已存 `sr_report.txt`）
  1. 长度对照表：ms→帧（48000/44100）、重采样输出长度推导
  2. 各环波形摘要：n/peak/rms/mean + 波形字节 FNV-1a 哈希，两种采样率各一份
  3. 两种采样率下各环耗时（0.25 s 缓冲 × 2000 次）
  4. 旧接口兼容清单：逐条调用旧 API 并验证返回值
- `make check-compat` — 重新构建五个模块的自测，与 `golden/` 下统一前的
  输出逐字节比较（cmp），全部相同才算通过。
