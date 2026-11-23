# fft

- `make test` 跑自检：`fft`（幅度/Parseval/往返）与 `spectral`（Hann 谱平均、设参公式、1/√K 收敛，14 项）。
- 谱平均帧数不靠猜：推导、低 SNR 偏离与设参依据见 [SPECTRAL_AVERAGING.md](SPECTRAL_AVERAGING.md)。
- 蒙特卡洛验证与图表：`tools/spectral_mc.py`、`tools/make_figures.py`（输出在 `../artifacts/spectral_v4/`）。

## 设参接口

```c
unsigned k = spectral_frames_needed(/*j_max=*/0.10f, /*snr_db=*/0.0f);
/* K = ceil(1.5 * C(rho) / J_max^2), rho = (N/3)*10^(SNR/10), C=(1+2rho)/(1+rho)^2 */
```

保守取 Hann 50% 重叠膨胀因子 1.5（实测真实 IF 随 ρ 落在 1.03–1.33，1.5 恒为安全上界）。
