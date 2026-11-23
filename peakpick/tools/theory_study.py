#!/usr/bin/env python3
# 理论预研：四种窗的窗谱形状参数
#  - 主瓣第一零点宽度（bin）
#  - 两等幅谱线在中点处由凹陷变单峰的临界间距（g'' 拐点法，bin）
#  - 最高旁瓣电平（dB）
# 用高倍零填充 FFT 逼近连续 DTFT。
# -I 下用户 site-packages 不在 sys.path，显式补入（不依赖 cwd，防同名模块投毒）。
import sys
sys.path.insert(0, "/home/ys/.local/lib/python3.12/site-packages")
import numpy as np

M = 2 ** 22
WINN = {
    "rect":     np.ones,
    "hann":     np.hanning,
    "hamming":  np.hamming,
    "blackman": np.blackman,
}

def window_spectrum(w):
    n = len(w)
    spec = np.fft.fftshift(np.fft.fft(w, M))
    dc = spec[M // 2]                      # fftshift 后偶数长度 FFT 的 DC 位置
    g = spec[M // 2:] / abs(dc)
    # 频率换算为 bin（相对窗长 N：1 bin = 1/N 周/样本）
    u = np.arange(len(g)) / M * n
    return u, np.abs(g)

def first_null(u, mag):
    # 主瓣第一零点：第一个严格局部极小（Hamming 无精确零点，取主瓣末端谷值）
    loc = np.where((mag[1:-1] < mag[:-2]) & (mag[1:-1] < mag[2:]))[0] + 1
    return u[loc[0]]

def peak_sidelobe(mag):
    # 主瓣外的最大局部极大
    loc = np.where((mag[1:-1] > mag[:-2]) & (mag[1:-1] > mag[2:]))[0] + 1
    edge = np.where((mag[1:-1] < mag[:-2]) & (mag[1:-1] < mag[2:]))[0] + 1
    outside = loc[loc > edge[0]]
    return 20 * np.log10(mag[outside].max())

def crit_spacing(u, g):
    # g(u)=|W(u)|/|W(0)|；中点凹陷消失等价于 g 在半间距处经过凸凹拐点，
    # 即 g'' 由负转正的第一个正频率零点。临界间距 = 2*u_inflect。
    d2 = np.gradient(np.gradient(g, u), u)
    cross = np.where(np.diff(np.sign(d2[1:])) > 0)[0]
    return 2 * u[cross[0] + 1]

print(f"{'window':9s} {'first_null(bin)':>15s} {'crit_d(bin)':>12s} {'PSL(dB)':>9s}")
for name, fn in WINN.items():
    w = fn(4096)
    u, mag = window_spectrum(w)
    print(f"{name:9s} {2*first_null(u, mag):15.3f} {crit_spacing(u, mag):12.3f} "
          f"{peak_sidelobe(mag):9.2f}")
