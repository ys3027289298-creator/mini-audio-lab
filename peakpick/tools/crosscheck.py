#!/usr/bin/env python3
"""C/Python 交叉校验：C 可执行文件的 raw_* 原始计数 == Python 引擎同工况计数。

C 与 Python 是两套独立 FFT/峰值检测实现；原始计数逐项相等才能证明
README 里的误检/漏检结论不是单一实现的产物。
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/ys/.local/lib/python3.12/site-packages")
import numpy as np

from common import FS, window, pick_peaks, WIN_NAMES

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def py_raw():
    """复现 src/peakpick.c 全部工况（固定绝对频率，NFFT 与 C 字面量一致）。"""
    out = {}

    def peaks(win, L, f0, f1, nfft):
        n = np.arange(L)
        x = (np.cos(2 * np.pi * f0 * n / FS)
             + np.cos(2 * np.pi * f1 * n / FS))
        spec = np.fft.rfft(x * window(win, L), nfft)
        return pick_peaks(np.fft.rfftfreq(nfft, 1.0 / FS), np.abs(spec))

    # 试验 1：单音 3000Hz, L=512, NFFT=2048
    L, nfft, f0, abin = 512, 2048, 3000.0, FS / 512
    for win in WIN_NAMES:
        pk = peaks(win, L, f0, f0, nfft)
        pk = pk[np.abs(pk - f0) <= 8 * abin]
        out[f"raw_sidelobe_fp_{win}"] = int(
            sum(abs(p - f0) > 0.5 * abin for p in pk))

    # 试验 2/3：blackman L=256
    abin256 = FS / 256
    for tag, f1 in (("df300", 3300.0), ("df400", 3400.0)):
        pk = peaks("blackman", 256, 3000.0, f1, 1024)
        pk = pk[(pk >= 3000 - 8 * abin256) & (pk <= f1 + 8 * abin256)]
        out[f"raw_blackman_L256_{tag}_peakcount"] = len(pk)

    # 试验 4：blackman L=2048, NFFT=8192
    abin2048 = FS / 2048
    pk = peaks("blackman", 2048, 3000.0, 3400.0, 8192)
    pk = pk[(pk >= 3000 - 8 * abin2048) & (pk <= 3400 + 8 * abin2048)]
    out["raw_blackman_L2048_peakcount"] = len(pk)
    out["raw_blackman_L2048_fp"] = int(
        sum(abs(p - 3000.0) > 0.5 * abin2048
            and abs(p - 3400.0) > 0.5 * abin2048 for p in pk))
    return out


def main():
    expected = py_raw()
    proc = subprocess.run([os.path.join(ROOT, "peakpick")],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        print("C self-test reported failure")
        sys.exit(1)
    c_raw = {}
    for line in proc.stdout.splitlines():
        if line.startswith("raw_") and "=" in line:
            k, v = line.split("=")
            c_raw[k.strip()] = int(v)
    fails = 0
    print(f"{'metric':38s} {'C':>6s} {'Py':>6s}")
    for k, py_v in expected.items():
        c_v = c_raw.get(k)
        ok = c_v == py_v
        print(f"{k:38s} {str(c_v):>6s} {str(py_v):>6s}  {'OK' if ok else 'FAIL'}")
        if not ok:
            fails += 1
    print("crosscheck_ok =", fails == 0)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
