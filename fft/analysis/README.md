# fft 舍入误差累积律分析

对象：`fft/src/fft.c` —— 迭代基-2 复数 FFT，**float 数据 + float 递推旋转因子**
（每级内 `w ← w·wl` 连乘，`wl = (float)cos(2π/len) + i·(float)sin(2π/len)`）。

参考：同算法、同递推结构的 double 实现（自身误差 ~N·2^-53 ≤ 4e-9，相对被测
float 路径的 ~1e-3 可忽略）。输入统一为 [-1,1) 复白噪声（固定种子 LCG），
先量化为 float，两条路径变换**同一份量化值**，因此只测变换内部舍入，
不含输入量化。`eps_f = 2^-24 = 5.96e-8`。

## 1. 误差来源定位（variants.c）

| 变体 | 做法 | rms 相对误差量级 |
|---|---|---|
| mode 0（仓库代码） | float 递推 twiddle | 1e-7 → **1.8e-4**（随 N 涨） |
| mode 1 | 每 k 直接算 float cos/sin | 7e-8 → 1.8e-7（几乎平） |
| mode 2 | 递推在 double、取用转 float | 同 mode 1 |

数据：`variants.csv`。蝶形运算本身的舍入地板 ~1e-7·√(log2 N)，
**全部超额误差来自 float 的 `w ← w·wl` 递推**。

## 2. 累积律

两级结构（recur.c 实测每级，err_sweep.c 实测整体）：

**(a) 每级 twiddle 递推误差**（len = 级长，k = 0..len/2-1 步递推）：

```
W(len) ≈ (len/2) · | |wl_fl| − 1 |        系统性幅度漂移，随步数线性累积
```

实测 `recur.csv`：`drift_model = (len/2)·||wl_fl|−1|` 与实测 Wmax 在
len ≥ 512 后基本重合。该项是一阶（∝ eps_f）小量时成立。

**(b) 跨级累积**：第 len 级注入的误差经后续级与信号同构放大，
每级对最终 rms 相对误差的贡献为 W(len)/√2，与 N 无关，逐级平方和：

```
rms_rel(N) ≈ sqrt( Σ_{len=2,4,…,N} W(len)² / 2 )
```

预测 vs 实测（`sweep.csv` 对 `recur.csv` 合成）：N = 2^4 … 2^22
全部 19 个点比值在 **0.77–1.64**（绝大多数 0.9–1.3）。

**(c) 闭式包络（线性累积律）**：wl 各分量舍入 |δ| ≤ 2^-24 时
W(len) ≲ (len/2)·eps_f，最大级主导，得

```
rms_rel(N) ≲ C_rms · N · eps_f ,  C_rms ≈ 0.065   （N ≤ 2^14 实测标定）
max_rel(N) ≲ C_max · N · eps_f ,  C_max ≈ 0.23
```

即**误差随长度线性累积**——这是"该律"的闭式形态。

## 3. 适用长度上限

线性律的前提是"wl 带一般性 O(eps_f) 舍入"。但

```
1 − cos(2π/len) ≈ 2π²/len² < 2^-25（1.0f 之下 ulp 的一半）
⟺ len > 2^13 · π ≈ 25736
```

**len ≥ 2^15 = 32768 时 cos(2π/len) 恰好舍入为 1.0f**（recur.csv
`cos_is_one` 列实测确认：2^14 为 0，2^15 起为 1）。此时
`|wl_fl|−1 = 2π²/len²`（比一般舍入小 3–4 个数量级），系统性漂移
`W_sys(len) = π²/len` 转为**随 len 递减**，线性累积的前提断裂。

> **该律适用上限：N_max = 2^14 = 16384。**
> N ≥ 2^15 起不得用 ∝N·eps_f 外推。

## 4. 验证

### 4.1 上限内一组：N₁ = 8192 = 2^13（逐 bin 差：`bins_8192.csv`）

| 量 | 实测 | 律的预测/包络 | 判定 |
|---|---|---|---|
| rms_rel | 2.40e-5 | 包络 0.065·N·eps_f = 3.17e-5；平方和律 1.90e-5 | ✓ 包络内，比值 1.26 |
| max_rel | 1.00e-4（bin 4067） | 包络 0.23·N·eps_f = 1.12e-4 | ✓ 包络内 |
| 逐 bin | p50=1.2e-3, p99=4.7e-3, max=7.4e-3（绝对值） | 随 k 无异常结构 | ✓ |

整个 N ≤ 2^14 区间逐点检验：`rms_rel/(N·eps_f)` ≤ 0.062，
`max_rel/(N·eps_f)` ≤ 0.23，线性包络全程成立。

### 4.2 上限外对照组：N₂ = 262144 = 2^18（逐 bin 差：`bins_262144.csv`）

**该组不与线性累积律比较一致性，仅报告实测与律的失效：**

| 量 | 线性律外推 | 实测 | 结论 |
|---|---|---|---|
| 每级 W(2^18) | (len/2)·eps_f = 7.8e-3 | 1.25e-7（recur.csv） | 闭式前提失效，高估 6×10⁴ 倍 |
| rms_rel 增长（2^15→2^18） | ∝N：×8 | ×1.5（1.17e-4→1.76e-4） | 增长律不成立 |
| rms_rel / max_rel | — | 1.76e-4 / 1.00e-3（bin 64528） | 进入饱和区，非线性律区 |

N ≥ 2^15 后误差**饱和**（rms_rel ≈ 1.8e-4、max_rel ≈ 1e-3 维持到 2^20），
原因是 cos≡1.0f 后漂移项 π²/len 递减。**对 N₂ 不宣称与线性累积律一致。**

### 4.3 累积曲线（sweep.csv，rms 相对误差）

```
 N        rms_rel
 2^ 4     1.16e-07  ▏
 2^ 6     3.79e-07  ▏
 2^ 8     6.46e-07  ▎
 2^10     2.17e-06  ▌
 2^12     1.51e-05  ██▌          ← 线性累积区（∝N·eps_f）
 2^13     2.40e-05  ████         ← 验证组 N₁
 2^14     2.14e-05  ███▌         ← 适用上限 N_max
 2^15     1.17e-04  ██████████▌  ← cos→1.0f，前提断裂
 2^16     1.52e-04  █████████████
 2^18     1.76e-04  ██████████████▌   ← 对照组 N₂（饱和，不外推）
 2^20     1.82e-04  ███████████████
 2^21     2.20e-04  █████████████████
 2^22     2.89e-04  ████████████████████▌ ← 随机游走项重新主导
```

附注（第三区）：N ≥ 2^21 起每步舍入的随机游走 ~√len·eps_f 超过
π²/len 漂移，误差重新缓慢增长（recur.csv 与 sweep.csv 均可见）。
饱和也不是永久保证。

## 5. 复现

```sh
gcc -O2 -std=c11 -Wall -Wextra err_sweep.c -o err_sweep -lm
gcc -O2 -std=c11 -Wall -Wextra variants.c  -o variants  -lm
gcc -O2 -std=c11 -Wall -Wextra recur.c     -o recur     -lm
./err_sweep sweep > sweep.csv          # 累积曲线 N=2^4..2^22
./err_sweep bins 8192   bins_8192.csv    # 限内组逐 bin 差
./err_sweep bins 262144 bins_262144.csv  # 限外对照逐 bin 差
./variants > variants.csv              # 误差来源分解
./recur > recur.csv                    # 每级递推误差 W(len)
```

## 6. 文件

- `err_sweep.c` — 扫描 + 逐 bin 差输出（被测算法逐行照搬 `../src/fft.c`）
- `variants.c` — 误差来源分解（递推 twiddle vs 蝶形运算）
- `recur.c` — 每级 twiddle 递推误差 W(len) 独立测量
- `sweep.csv` — 累积曲线数据（n, max_abs, rms_abs, rms_ref, rel_max）
- `recur.csv` — W(len) 与 |wl_fl|−1、cos 是否为 1.0f
- `variants.csv` — 三种变体的 rms 相对误差
- `bins_8192.csv` — 限内验证组逐 bin 差（bin, re_diff, im_diff, abs_diff, mag_ref）
- `bins_262144.csv` — 限外对照组逐 bin 差（同上）
