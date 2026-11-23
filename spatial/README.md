# spatial

`make test` 跑自检（`spatial` = 等功率声像/距离衰减，`early_reflection` = 预延迟早期反射）。

## early_reflection

### 预延迟

反射抽头整体抬到预延迟之后，和干声拉开：

```
t_i = pre_ms/1000 + dist_i / c        [s]
D_i = t_i * fs                        [samples]
n_i = floor(D_i),  mu_i = D_i - n_i
```

默认 `fs = 48000 Hz`、`c = 343 m/s`、`pre_ms = 5 ms`（240 样本）。反射相对直达声的距离
`dist_i`（米）换算成抽头时刻，再用分数延迟插值落到非整数样本上。

反向换算（抽头时刻 -> 距离）：

```
dist_i = ((n_i + mu_i)/fs - pre_ms/1000) * c
```

`D_i <= 0`（含 `dist = 0, pre_ms = 0`）报无效；反射的第一个子抽头落到偏移 0 或负偏移
（即与直达声重合或提前）也报无效。有效构建下所有抽头偏移 `>= 1`；默认配置的
`first_offset = 588`，与 240 样本预延迟之间留有间隔。

### 分数延迟插值

奇数阶 Lagrange 插值器，阶数 `N`，`n = (N+1)/2`，抽头 `j = -(n-1)..n`：

```
h_j(mu) = prod_{i != j} (mu - i)/(j - i)
y[k]    = sum_j h_j * x[k - n_i - j]
```

成立两个恒等式（自检里按 0.01 步长扫全部 mu 验证）：

```
sum_j h_j     = 1      -> 直流增益 1
sum_j j*h_j   = mu     -> 一阶矩 = mu
```

所以每个反射的抽头质心 `sum (n_i+j) h_j / sum h_j = n_i + mu_i = D_i`，与解析延迟零偏差
（浮点残差 ~1e-8 样本）。

### 插值误差上界

输入频谱限制在 `fmax` 以内、幅度不超过 `amp` 时，分数延迟滤波器的频响误差满足

```
|E(e^jw)| = w^(N+1)/(N+1)! * prod_j |mu - j|,   w = 2*pi*fmax/fs
err <= amp * w^(N+1)/(N+1)! * prod_j |mu - j| <= amp * w^(N+1)/(N+1)! * P_N
```

`P_N = max_mu prod_j |mu-j| = prod_j |0.5 - j|`（在 `mu = 1/2` 取到）：

| 阶数 N | 抽头数 | P_N |
|---|---|---|
| 1 | 2 | 0.25 |
| 3 | 4 | 0.5625 |
| 5 | 6 | 3.515625 |

`amp = 1` 时的上界（fs = 48 kHz）：

| fmax | N=1 | N=3 | N=5 |
|---|---|---|---|
| 1 kHz | 2.14e-3 | 6.88e-6 | 2.46e-8 |
| 4 kHz | 3.43e-2 | 1.76e-3 | 1.01e-4 |
| 8 kHz | 1.37e-1 | 2.82e-2 | 6.44e-3 |

误差随阶数按 `w^(N+1)/(N+1)!` 下降（8 kHz 处 1 阶 -> 3 阶约 4.86 倍，3 阶 -> 5 阶约 4.38 倍）。

相位延迟偏差由上界推出：若 `|R - e^{-j w mu}| <= e`，则
`|tau_p - mu| <= asin(e/(1-e))/w`。自检在 4 kHz 实测：N=1 为 3.9e-3 样本，N=3 为 2.3e-4，
N=5 为 1.4e-5，均在界内。

### 自检输出

- `tap{r}_dist_m / _delay_ms / _delay_samples / _int / _frac / _dist_back_m`：抽头换算表
- `tap{r}_centroid_samples / _centroid_dev`：质心与解析延迟的偏差
- `tap{r}_ir{k}`：每个反射的分数延迟子抽头（偏移: 权重），拼成脉冲响应
- `ir_len`、`ir_first_offset`、`predelay_samples`：IR 长度和预延迟间隔
- `interp_bound_o{N}` / `interp_err_o{N}`：各阶误差上界与 7.9 kHz 带限测试信号实测误差
- `err_bound_o{N}_f{F}`：上界表
- `zero_delay_invalid`、`neg_dist_invalid`、`tap_on_direct_invalid` 等无效输入用例
