# Minimum / maximum code and the dither-safe peak are the two things a
# game-audio exporter has to get right. These are checked here rather than
# buried in the self-test so they can be inspected on their own.
stops = []

def expect(name, got, want):
    st = "ok " if got == want else "FAIL"
    stops.append(got == want)
    print(f"  [{st}] {name}: {got} (want {want})")

import math
S = 32768.0
def q(x, d=0.0):
    v = x * S + d
    if v + 0.5 >= 32768.0: return 32767
    if v + 0.5 < -32768.0: return -32768
    return math.floor(v + 0.5)

print("16-bit two's-complement rails")
expect("-1.0 -> -32768", q(-1.0), -32768)
expect("+1.0 -> 32767", q(1.0), 32767)
expect("0.0 -> 0", q(0.0), 0)
expect("32766/32768 + up to 1 LSB dither stays in range",
       max(q(32766 / S, d) for d in (-1.0, -0.5, 0.0, 0.5, 0.999)), 32767)
expect("1.0 + 0.999 dither saturates, does not wrap",
       q(1.0, 0.999), 32767)
expect("-1.0 - 0.999 dither saturates, does not wrap",
       q(-1.0, -0.999), -32768)

print("dither amplitude bounds, verified against the exported files")
import wave
import numpy as np
for variant in ("plain", "dither"):
    w = wave.open(f"evidence/wav/tone60_{variant}.wav")
    qq = np.frombuffer(w.readframes(w.getnframes()), dtype="<i2").astype(np.int32)
    w.close()
    x = np.frombuffer(open(f"evidence/wav/tone60_{variant}.wav.f32", "rb").read(),
                      dtype="<f4").astype(np.float64)
    e = (qq / S - x) * S
    print(f"  {variant:6s} err rms {np.sqrt(np.mean(e*e)):.3f} LSB, "
          f"peak {np.abs(e).max():.3f} LSB, codes [{qq.min()}, {qq.max()}]")
    if variant == "plain":
        expect("plain error within +/-0.5 LSB", bool(np.abs(e).max() <= 0.50001), True)
    else:
        expect("TPDF error within +/-1.5 LSB", bool(np.abs(e).max() < 1.50001), True)
        expect("no dither-only clipping at -60 dBFS",
               int(np.sum((qq == 32767) | (qq == -32768))), 0)

# Variance budget: total dither error measured on the export, and the TPDF
# source variance from an independent Python port of the documented PCG32
# (same algorithm as dither_core.c). The two add to the 0.25 LSB^2 figure.
M64 = (1 << 64) - 1

class Pcg32:
    def __init__(self, seed, stream):
        self.state = seed & M64
        self.inc = ((stream << 1) | 1) & M64
        self.next()
    def next(self):
        old = self.state
        self.state = (old * 6364136223846793005 + self.inc) & M64
        xs = ((((old >> 18) ^ old) >> 27)) & 0xffffffff
        rot = (old >> 59) & 31
        return ((xs >> rot) | (xs << ((32 - rot) & 31))) & 0xffffffff
    def uni(self):
        return (self.next() >> 8) * (1.0 / 16777216.0) - 0.5

r = Pcg32(1234, 0x5eed)
tri = np.array([r.uni() + r.uni() for _ in range(200000)])
print(f"  TPDF source variance {tri.var():.4f} LSB^2 (two uniforms => 1/6 = 0.1667), "
      f"range [{tri.min():.3f}, {tri.max():.3f}]")
expect("TPDF source variance matches D^2/6", bool(abs(tri.var() - 1 / 6) < 0.01), True)
expect("TPDF support within (-1, 1) LSB", bool(tri.min() > -1.0 and tri.max() < 1.0), True)
try:
    wd = wave.open("evidence/wav/tone60_dither.wav"); qd = np.frombuffer(
        wd.readframes(wd.getnframes()), dtype="<i2").astype(np.int32); wd.close()
    x = np.frombuffer(open("evidence/wav/tone60_dither.wav.f32", "rb").read(),
                      dtype="<f4").astype(np.float64)
    print(f"  total dither error variance {float(np.var((qd / S - x) * S)):.4f} LSB^2 "
          f"(=> 0.25)")
    expect("total dither error variance matches 0.25",
           bool(abs(float(np.var((qd / S - x) * S)) - 0.25) < 0.01), True)
except FileNotFoundError:
    print("  (exported WAVs absent; skipping export-side variance check)")

print("all_ok=" + ("1" if all(stops) else "0"))
raise SystemExit(0 if all(stops) else 1)
