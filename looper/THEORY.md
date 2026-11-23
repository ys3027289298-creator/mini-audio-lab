# 循环接缝的相位连续性条件与泄漏量估计

记号：采样率 $f_s$，循环缓冲长度 $N$（样本），信号 $x[i]\;(i=0,\dots,N-1)$，
循环播放时索引按模 $N$ 回绕，接缝指 $x[N-1]\to x[0]$ 的回绕点。

## 1. 两个必须区分的接缝量

讨论"咔哒声"时常见的混淆来自把两个不同的量都叫作"跳变"：

| 量 | 定义 | 含义 |
|---|---|---|
| 相位残差 $\Delta$ | $x[0]-x_{\text{synth}}(N)$：与"信号若不回绕、继续走到第 $N$ 个样本"之差 | 真正的相位断裂量；为零当且仅当循环相位连续 |
| 接缝步长 $s$ | $x[0]-x[N-1]$：回绕瞬间实际发生的样本跳变 | 咔哒声的直接载体；与 $\Delta$ 只差一个平滑单样本步长 |

对测试音 $x[i]=\cos(2\pi q\,i/N)$（周期 $P=N/q$，频率 $f_s/P$），
$x_{\text{synth}}(N)=\cos(2\pi q)=1$，故

$$\Delta=1-\cos(2\pi q),\qquad
s=x[0]-x[N-1]=\Delta+O\!\Big(\frac{2\pi q}{N}\Big).$$

$O(2\pi q/N)$ 就是信号自身一个样本的平滑运动量（测试参数下约 0.144）。
因此整周期循环时 $s$ 并不为零，但它与缓冲内部一阶差分同量级——**判据不是
$s=0$，而是 $s$ 相对内部差分是否是离群值**。自检中：整周期 $s=0.0108$，
内部最大差分 $0.1471$（同量级）；半 bin 失配 $s=1.9896$，是内部的 13.5 倍（离群）。

## 2. 基频整周期条件

循环一阶差分 $d[i]=x[i]-x[(i-1)\bmod N]$。对单音 $x[i]=A\cos(2\pi f_0 i/f_s+\varphi)$，
回绕后相位连续要求

$$2\pi f_0 N/f_s=2\pi q,\ q\in\mathbb Z^+
\;\Longleftrightarrow\; N=qP,\ P=f_s/f_0
\;\Longleftrightarrow\; f_0=q\,\frac{f_s}{N}.$$

即**循环长度必须是基频周期的整数倍**，等价于基频恰好落在 DFT 第 $q$ 个 bin 上。
该条件下 $\Delta=0$，$d[0]$ 与内部差分严格同分布，循环信号在 $N$ 点 DFT 中
只占 bin $q$ 与 $N-q$。多谐信号 $x[i]=\sum_h A_h\cos(2\pi h f_0 i/f_s+\varphi_h)$
条件相同：一次裁掉整数个基频周期，所有谐波同时闭合。

不满足时写 $q=q_0+\varepsilon$，$q_0=\operatorname{round}(q)$，
$\varepsilon\in(-\tfrac12,\tfrac12]$ 为分数 bin 偏移。取 $\varphi=0$：

$$\Delta=A(1-\cos 2\pi\varepsilon)=2A\sin^2(\pi\varepsilon).$$

由 $|\sin u|\le |u|$ 得到紧二次界与最坏相位线性界：

$$|\Delta|\le 2\pi^2 A\varepsilon^2,\qquad |\Delta|\le 2\pi A|\varepsilon|.$$

自检参数 $q=23.5,\ \varepsilon=0.5,\ A=1$：理论 $\Delta=2$，
紧界 $2\pi^2\varepsilon^2=4.935$、线性界 $2\pi\varepsilon=3.142$，
两个界都被数值校验（`residual_inside_bounds`）。

## 3. 泄漏量估计（严格离散版本）

### 3.1 关键事实：缓冲内 DFT 没有泄漏

对 $N$ 点缓冲做 DFT，分数 bin 余弦的能量只落在两个 bin 上（$q=23.5$ 时
恰为 bin 23、24，各占一半），**缓冲内谱是完全稀疏的**。"频谱泄漏"不在静态
DFT 里——它在回绕动作发生时才产生。因此正确的分析对象是循环一阶差分 $d[i]$：
缓冲内部 $d[i]$ 只是音调的平滑导数，而 $d[0]=s$ 携带接缝跳变。

### 3.2 DFT 导数恒等式

循环差分对应 DFT 域乘以 $1-e^{-j2\pi k/N}$：

$$D[k]=(1-e^{-j2\pi k/N})X[k],\qquad
|D[k]|=2\big|\sin(\pi k/N)\big|\,|X[k]|.$$

另有直流望远镜恒等式 $\sum_i d[i]=0$（自检 `click_identity_dc`）。
去掉两个音调 bin 后，$D$ 的能量就是接缝脉冲铺出的宽带谱。自检在
bin 30–80 上对每个 bin 验证导数恒等式（`derivative_identity`，相对误差
$10^{-8}$ 量级），并对照两种循环：

| | 整周期 $q=24$ | 半 bin $q=23.5$ |
|---|---|---|
| 接缝样本差分能量 $d[0]^2$ | $1.17\times10^{-4}$ | $3.96$ |
| 内部最大单样本差分能量 | $2.16\times10^{-2}$ | $2.08\times10^{-2}$ |
| bin 30–80 宽带底噪 | $1.5\times10^{-25}$ | $3.95\times10^{2}$ |

宽带底噪相差 $2.6\times10^{27}$ 倍（`click_floor_present`）：接缝脉冲的谱
近乎白色，这就是静音背景下"咔"声格外可闻的原因。接缝单样本占全部差分能量的
27%，且是任何内部样本的 190 倍以上（`seam_energy_share`）。

### 3.3 连续模型：锯齿修正项与闭式泄漏

工程估计用连续周期延拓视角。把非整周期信号写成整周期分量加修正项，
修正项的回绕跳变为 $-\Delta$，形状是斜率 $\Delta/N$ 的斜坡（锯齿）$c[i]=\Delta i/N$。
其 DFT（$k\ne0$）为

$$C[k]=-\frac{\Delta}{1-e^{-j2\pi k/N}},
\qquad |C[k]|=\frac{|\Delta|}{2|\sin(\pi k/N)|}
\;\approx\;\frac{|\Delta|N}{2\pi k}\quad(k\ll N).$$

低频按 $1/k$ 缓降、高频趋于常数 $|\Delta|$：平坦宽带底噪，与 3.2 的离散
观测一致。由 Parseval，锯齿总能量

$$E_c=\sum_{i=0}^{N-1}\!\Big(\frac{\Delta i}{N}\Big)^2
=\Delta^2\frac{(N-1)(2N-1)}{6N}\approx\frac{\Delta^2 N}{3}.$$

相对基频能量 $E_s=A^2N/2$：

$$\frac{E_c}{E_s}\approx\frac{2\Delta^2}{3A^2}
\;\lesssim\;\frac{8\pi^2}{3}\varepsilon^2.$$

该比值与 $N$ 无关：加长缓冲不会让咔哒声变小，因为 $\Delta$ 本身不变。

## 4. 两种样本级修复与连续性证明

### 4.1 周期裁剪（`loop_trim_period` / `loop_trim_fractional`）

估计基频周期 $\hat P$（自相关，搜索区间给定），保留
$N'=\lfloor N/\hat P\rfloor\hat P$ 个样本（整数周期）或按真实分数周期取整。
样本级裁剪只能落在整数索引上，残余相位残差受单样本运动界约束：

$$|\Delta'|\lesssim \frac{2\pi q_0}{N'}.$$

自检中分数周期裁剪到 1002 样本（23 个周期），$\Delta'=4.7\times10^{-4}$，
仅为该界（0.144）的 0.33%；接缝步长与内部差分之比 0.106，泄漏能量占比
$7.8\times10^{-5}$。代价：损失至多一个周期的样本，且要求信号近似周期。

### 4.2 接缝钉值的 $\cos^2$ 交叉淡化（`loop_crossfade`）

不假设周期性。令 $\delta=x[0]-x[N-1]$，在回绕接缝两侧各取 $M$ 个样本，窗函数

$$w_j=\cos^2\!\frac{\pi j}{2(M-1)},\qquad j=0,\dots,M-1,$$

满足 $w_0=1$、$w_{M-1}=0$，输出

$$y[j]=x[j]-\tfrac{w_j}{2}\delta,\qquad
y[N-M+j]=x[N-M+j]+\tfrac{1-w_j}{2}\delta.$$

**C0 严格为零**：接缝两端 $y[0]=x[0]-\delta/2$、
$y[N-1]=x[N-1]+\delta/2$ 恒等（自检 `xfade_c0_exact`，残差 $<10^{-12}$）。
**远端零修正**：$w_{M-1}=0$、$1-w_0=0$，样本 $M-1$ 与 $N-M$ 原样保留，
中间样本 $M,\dots,N-M-1$ 完全不触碰（`xfade_leg_rejoin`、
`xfade_interior_untouched`）。**C1 不失控**：淡化区内每样本最大修正量
$\le\delta/[2(M-1)]$，取 $M\approx P$ 时接缝差分回到信号自身尺度
（`xfade_no_click`，接缝/内部差分之比 $\le2$）。宽带底噪相对未处理下降约
8.7 倍（`xfade_leak_down`）；残余来自淡化区自身的缓变，而非脉冲。

### 4.3 关闭处理

`enable=0` 走 `loop_passthrough`，逐样本拷贝，自检用 `memcmp` 证明输出与
输入位一致（`off_passthrough`）；任何非零启用值在直通函数中返回 `-1`，
杜绝"半处理"。float API 的全部性质单独走一遍（`float_*` 系列检查），
双精度缓冲仅用于理论证明的数值工作。
