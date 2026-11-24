# 后台预算自适应降级模拟 (budget_sim)

同一段 2 秒、48 kHz 立体声"音乐"信号依次穿过与仓库五个 C 模块一一对应的处理环节，
模拟应用切到后台、处理预算被压缩时，各环节按**频谱掩蔽贡献排序**动态降低精度的过程。

## 运行

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install numpy matplotlib
python -I budget_sim.py            # 输出默认写到 ../artifacts/budget_sim/
python -I budget_sim.py --out DIR  # 自定义输出目录
```

## 处理链与精度档位

| 环节 (对应 C 模块) | 高精度 | 中精度 | 低精度（降级手段） |
|---|---|---|---|
| mixer (mixer.c) | f64，连续增益斜坡 | f32，8 级增益量化 | f16 累加，4 级硬件增益阶梯 |
| resampler (resampler.c) | 16 抽头加窗 sinc | 线性插值 | 最近邻（产生镜像/混叠频率） |
| spatial (spatial.c) | f64 直接三角函数 | 128 点 pan LUT | 3 点 pan LUT（台阶式声像） |
| limiter (limiter.c) | 逐样本包络 | 16 样本保持、10 bit 系数量化 | 64 样本保持、3 bit 系数量化 |
| fft (fft.c) | 1024 点 STFT、全频带 | 512 点、截到 10 kHz | 128 点、截到 4 kHz |

每档有抽象计算成本（见脚本 `COST` 表），三档预算为 32 / 20 / 10。

## 掩蔽贡献如何估算

1. 全高精度参考输出做 2048 点 STFT，幅度归一化到峰值 90 dB SPL；
2. 按 Bark 刻度做 Johnston 扩散函数，功率域求和并加 9 dB 掩蔽偏移，
   与 Terhardt 绝对听阈（ATH，限幅在 20 Hz–16 kHz 有效区间）取最大值，
   得到逐帧逐频点的掩蔽阈值；
3. 每个环节**单独**降到低精度（其余保持高精度）跑一遍，取其误差谱，
   误差超过掩蔽阈值部分的均值（dB）即该环节的**掩蔽贡献**——
   可闻误差越大，贡献越大，越该被保护、越晚降级；
4. 贪心分配：预算超限时，按贡献从小到大逐环节降档，直到成本不超预算。

另实现一个**符号取反的错误估算器**：把"落在阈值之下、被掩蔽的误差"
（阈值与误差的差值均值）当成贡献。这等价于把掩蔽效应算反——
被掩蔽得越好的环节反而被认为越重要。

## 反转判定

最低预算档再用翻转后的排序分配一次。若正确模型中贡献最高（最不该降级）
的环节在翻转模型里被降级，判定为 `INVERTED_MASKING_DETECTED`，
并在控制台、`contributions.csv` 与 `fig2` 中把两套贡献值与降级顺序并列。

当前信号/参数下的结果：

- 正确降级顺序（先降→最后降）：
  spatial → mixer → limiter → fft → **resampler（最该保护）**
- 翻转顺序：**resampler（最先降）** → spatial → fft → mixer → limiter
- 判定：`INVERTED_MASKING_DETECTED`——最低档把重采样器降到 low，
  其在 1.1 kHz 与 8.5–11.7 kHz 掩蔽薄弱区的镜像误差本是最大的可闻失真来源。

## 输出文件 (`artifacts/budget_sim/`)

| 文件 | 内容 |
|---|---|
| `report.json` | 全部数值：成本表、两套贡献与顺序、三档精度/谱差、反转判定 |
| `contributions.csv` | 各环节掩蔽贡献值 × 两套降级顺序（并列） |
| `tier_precision.csv` | 三档预算下每环节实际精度与成本 |
| `spectral_diff.csv` | 三档输出在 8 个对数频带相对参考的谱差 |
| `fig1_masking_and_errors.png` | 掩蔽阈值曲线 + 五环节单环节误差谱 |
| `fig2_contribution_ranking.png` | 正确/翻转两套贡献排序与降级序号并列 |
| `fig3_tier_precision.png` | 三档 × 五环节实际精度热力表 |
| `fig4_spectral_diff_tiers.png` | 三档输出的频带谱差曲线 |
| `wav/` | 参考、5 个单环节探针、3 个预算档、翻转档共 10 个 16-bit WAV |

## 模型边界

- 掩蔽阈值只建模了同时掩蔽（频域扩散 + ATH），不含时域前向/后向掩蔽；
- "贡献"是对固定测试信号的测量值，换信号（尤其频谱形状不同）排序可能变化，
  这正是真实系统需要随帧重估的原因；
- 成本单位是相对抽象值，只用于驱动贪心分配，不代表真实 CPU 周期。
