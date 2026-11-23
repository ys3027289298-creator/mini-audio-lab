#!/usr/bin/env bash
# Reproduce the whole evidence bundle from scratch.
#   ./evidence/run_all.sh
# Rebuilds the binaries, exports the WAVs, runs the independent verifier,
# and verifies the two byte-identity claims:
#   1. dither on, same seed, run twice  -> identical WAV bytes
#   2. DITHER=off, any seed             -> identical to the plain-rounding
#      reference, byte for byte, and identical across seeds.
set -euo pipefail
cd "$(dirname "$0")/.."          # -> dither/

make -s clean
make -s
./dither | tee evidence/out/selftest.txt

mkdir -p evidence/wav evidence/out evidence/svg
FREQ=750   # exactly 64 samples/cycle at 48 kHz: coherent for any 2^n FFT
for mode in tone60 sweep gated; do
    ./dither_export "$mode" "evidence/wav/${mode}_dither.wav" "$FREQ" | tee "evidence/out/${mode}_dither.txt"
    DITHER=off ./dither_export "$mode" "evidence/wav/${mode}_plain.wav" "$FREQ" | tee "evidence/out/${mode}_plain.txt"
done

# Level sweep, straight from the C exporter: -6 .. -100 dBFS, both variants.
mkdir -p evidence/wav/levels
for lvl in -6 -20 -40 -60 -80 -90 -100; do
    TONE_DB=$lvl ./dither_export tone60 "evidence/wav/levels/l${lvl}_dither.wav" "$FREQ" >/dev/null
    TONE_DB=$lvl DITHER=off ./dither_export tone60 "evidence/wav/levels/l${lvl}_plain.wav" "$FREQ" >/dev/null
done

# 1) reproducible dithered bytes
./dither_export tone60 /tmp/repro_a.wav "$FREQ" >/dev/null
./dither_export tone60 /tmp/repro_b.wav "$FREQ" >/dev/null
cmp /tmp/repro_a.wav /tmp/repro_b.wav
echo "repro_same_seed_identical=1" | tee evidence/out/repro.txt
cmp /tmp/repro_a.wav evidence/wav/tone60_dither.wav
echo "repro_matches_bundle=1" | tee -a evidence/out/repro.txt

# 2) dither off == plain, independent of RNG state and seed
DITHER=off DITHER_SEED=1   ./dither_export tone60 /tmp/off_s1.wav "$FREQ" >/dev/null
DITHER=off DITHER_SEED=999 ./dither_export tone60 /tmp/off_s999.wav "$FREQ" >/dev/null
cmp /tmp/off_s1.wav /tmp/off_s999.wav
cmp /tmp/off_s1.wav evidence/wav/tone60_plain.wav
echo "dither_off_seed_independent=1" | tee -a evidence/out/repro.txt
echo "dither_off_matches_bundle=1" | tee -a evidence/out/repro.txt

# dither on with a different seed must NOT be byte-identical (the dither is live)
DITHER=on DITHER_SEED=777 ./dither_export tone60 /tmp/on_s777.wav "$FREQ" >/dev/null
if cmp -s /tmp/on_s777.wav evidence/wav/tone60_dither.wav; then
    echo "different_seed_changes_bytes=0" | tee -a evidence/out/repro.txt; exit 1
fi
echo "different_seed_changes_bytes=1" | tee -a evidence/out/repro.txt

rm -f /tmp/repro_a.wav /tmp/repro_b.wav /tmp/off_s1.wav /tmp/off_s999.wav /tmp/on_s777.wav

python3 evidence/check_rails.py | tee evidence/out/rails.txt
python3 evidence/verify_dither.py evidence/wav evidence/out/metrics.json | tee evidence/out/verify.txt
# fail the pipeline if a verification claim failed
grep -q "all verification claims hold" evidence/out/verify.txt
sha256sum evidence/wav/*.wav evidence/wav/*.f32 | tee evidence/out/sha256.txt
