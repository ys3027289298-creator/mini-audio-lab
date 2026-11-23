"""peakpick 共用分析引擎（噪声-free 确定性试验）。

信号模型
  fs=48000，两个等幅正弦 f0、f1=f0+df，加长度 L 的分析窗后做 ZP 倍零填充 FFT。
峰值检测（与 src/peakpick.c 同一规格）
  幅值(dB) 的严格局部极大，且高于「全谱最大值 - GATE_DB」门限。
  - 误检 FP：ROI 内被报告、但距任一真实频率超过 0.5 bin 的峰
  - 漏检 FN：真实频率 ±0.5 bin 内没有被报告峰
  - 两峰合并时局部极大只剩 1 个 -> 必然产生 1 个 FN
窗定义（对称型，C/Python 必须逐式一致）
  rect     w[n]=1
  hann     0.50 - 0.50 cos(2*pi*n/(L-1))
  hamming  0.54 - 0.46 cos(2*pi*n/(L-1))
  blackman 0.42 - 0.50 cos(2*pi*n/(L-1)) + 0.08 cos(4*pi*n/(L-1))
"""
import sys
sys.path.insert(0, "/home/ys/.local/lib/python3.12/site-packages")
import numpy as np

FS = 48000.0
K0 = 12                     # 低频音的固定 bin 索引（保证全 L 网格内远离 DC/Nyquist）
GATE_DB = 40.0              # 相对主峰的报告门限（-40 dB）
ZP = 4                      # 零填充倍数
FRAC_OFFSETS = np.array([0.0, 1/8, 2/8, 3/8, 4/8, 5/8, 6/8, 7/8])
PHASES = (0.0, np.pi / 2)   # 第二个音相对第一个的相位
WIN_NAMES = ("rect", "hann", "hamming", "blackman")


def window(name, L):
    n = np.arange(L)
    if name == "rect":
        return np.ones(L)
    a = 2.0 * np.pi * n / (L - 1)
    if name == "hann":
        return 0.50 - 0.50 * np.cos(a)
    if name == "hamming":
        return 0.54 - 0.46 * np.cos(a)
    if name == "blackman":
        return 0.42 - 0.50 * np.cos(a) + 0.08 * np.cos(2.0 * a)
    raise ValueError(name)


def two_tone(L, df, frac, phase):
    """生成一个双音试验片段。

    f0 固定在 bin 索引 K0+frac 上（而非固定绝对频率），这样两音的归一化
    几何只取决于 d_bin = df*L/fs，与 L 无关，可行性谓词随 L 严格单调。
    frac: f0 相对整数 bin 的分数偏移。
    """
    n = np.arange(L)
    bin_hz = FS / L
    f0 = (K0 + frac) * bin_hz
    f1 = f0 + df
    x = np.cos(2 * np.pi * f0 * n / FS) + np.cos(
        2 * np.pi * f1 * n / FS + phase)
    return f0, f1, x


def mag_spectrum(x, w, zp=ZP):
    L = len(w)
    nfft = 1 << int(np.ceil(np.log2(zp * L)))
    spec = np.fft.rfft(x * w, nfft)
    mag = np.abs(spec)
    freqs = np.fft.rfftfreq(nfft, 1.0 / FS)
    return freqs, mag


def pick_peaks(freqs, mag, gate_db=GATE_DB):
    """返回被报告峰的频率列表：严格局部极大 + 相对门限。"""
    if mag.max() <= 0:
        return np.empty(0)
    floor = mag.max() * 10.0 ** (-gate_db / 20.0)
    is_max = (mag[1:-1] > mag[:-2]) & (mag[1:-1] > mag[2:])
    idx = np.where(is_max & (mag[1:-1] >= floor))[0] + 1
    return freqs[idx]


def count_errors(peaks, truths, bin_hz):
    """按 0.5 bin 容差配对，返回 (FP, FN)。"""
    used = np.zeros(len(peaks), dtype=bool)
    fn = 0
    for t in truths:
        hit = np.where((~used) & (np.abs(peaks - t) <= 0.5 * bin_hz))[0]
        if len(hit) == 0:
            fn += 1
        else:
            used[hit[0]] = True
    fp = int((~used).sum())
    return fp, fn


def run_cell(win, L, df, zp=ZP):
    """一个 (窗, L, 间距) 单元：16 个试验，返回逐试验 (fp, fn) 列表。"""
    w = window(win, L)
    bin_hz = FS / L
    out = []
    for frac in FRAC_OFFSETS:
        for phase in PHASES:
            f0, f1, x = two_tone(L, df, frac, phase)
            freqs, mag = mag_spectrum(x, w, zp)
            peaks = pick_peaks(freqs, mag)
            # ROI：两个真实音左右各 8 bin，排除 DC/边缘与远端旁瓣的统计口径争议
            lo, hi = f0 - 8 * bin_hz, f1 + 8 * bin_hz
            peaks = peaks[(peaks >= lo) & (peaks <= hi)]
            out.append(count_errors(peaks, (f0, f1), bin_hz))
    return out


def cell_totals(win, L, df, zp=ZP):
    """单元汇总：总误检、总漏检、最差单试验值、全通过标志。"""
    trials = run_cell(win, L, df, zp)
    fps = [t[0] for t in trials]
    fns = [t[1] for t in trials]
    return {
        "fp_sum": sum(fps), "fn_sum": sum(fns),
        "fp_max": max(fps), "fn_max": max(fns),
        "clean": sum(fps) == 0 and sum(fns) == 0,
        "n_trials": len(trials),
    }


def measured_lcrit(win, df, mode="full", L_lo=32, L_hi=8192):
    """二分求最小整数窗长 L。

    mode="res" : 全部 16 个试验漏检为 0（分辨力临界，误检不管）
    mode="full": 全部试验误检、漏检均为 0（工程可行临界）
    单调性前提：clean 谓词随 L 单调（网格已验证 False*->True* 结构）。
    返回 None 表示搜索上界内仍不可行（如低旁瓣门限下的 rect/hann）。
    """
    def ok(L):
        t = cell_totals(win, L, df)
        return t["fn_sum"] == 0 if mode == "res" else t["clean"]

    if not ok(L_hi):
        return None
    lo, hi = L_lo, L_hi
    while lo < hi:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def gate_sensitivity(L=512, gates=(-20.0, -30.0, -40.0, -50.0, -60.0)):
    """单音（整数 bin K0 上）、各窗型在不同报告门限下的假峰数（ROI ±8 bin）。"""
    n = np.arange(L)
    bin_hz = FS / L
    f0 = K0 * bin_hz
    x = np.cos(2 * np.pi * f0 * n / FS)
    out = {}
    for win in WIN_NAMES:
        freqs, mag = mag_spectrum(x, window(win, L))
        row = {}
        for g in gates:
            pk = pick_peaks(freqs, mag, gate_db=-g)
            pk = pk[np.abs(pk - f0) <= 8 * bin_hz]
            row[int(g)] = int(sum(abs(p - f0) > 0.5 * bin_hz for p in pk))
        out[win] = row
    return out
