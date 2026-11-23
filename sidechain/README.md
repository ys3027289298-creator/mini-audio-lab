# sidechain

旁链（key 输入）触发检测器：attack/release 包络跟随 + 阈值比较 + hold 状态机。

## 构建与自检

```
make test
```

无参数运行内置自检（静默不误触发、步进信号单次触发、hold 期间不掉门、
释放后可重新触发、参数校验）。

## 流式模式

```
./sidechain in.f32 sr thresh_db attack_ms release_ms hold_ms env_out.f32 gate_out.u8
```

输入为单声道 raw float32；输出逐采样包络（float32）和门状态（uint8），
stdout 打印触发事件数与 engaged 采样数。

## 阈值分析（`analysis/`）

复现实验（需要 numpy / scipy / matplotlib）：

```
../.venv/bin/python -I analysis/gen_signal.py data
../.venv/bin/python -I analysis/analyze.py data results
```

或直接 `./run.sh`（含构建自检）。产物在 `results/`：

| 文件 | 内容 |
|---|---|
| `counts.json` | 分布参数、BIC、三档计数、11 点扫描、事件级结果 |
| `report.md` | 分布估计、理论曲线 vs 三档实测、超 30% 偏差归因 |
| `distribution.png` | 噪声/事件帧电平分布与拟合 |
| `curves.png` | FPR/FNR 理论曲线与实测点（三档 + 全扫描） |

方法要点：

- 触发电平按 10 ms 帧取包络最大值（阈值实际作用的量），dBFS 域分
  噪声帧/事件帧分别拟合单正态与 2 分量 GMM，BIC 选模。
- 理论曲线：FPR(T)=1−F₀(T)（噪声分布上尾），FNR(T)=F₁(T)（事件分布下尾）。
- 实测同时统计原始逐帧越限和 hold 状态机输出；理论只建模逐帧电平，
  hold 对噪声帧的额外占用单独报告。
- 偏差超 30% 必须归因（小样本精确 Wilson 区间 / 尾部结构性失配），
  归因不上脚本以非零码退出。
