#!/usr/bin/env bash
# Regenerate every artifact in artifacts/ from a clean build.
# Outputs: demo signals (both modes), length table, spectrograms, timing.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ART="$ROOT/artifacts"
LEG="$ART/demo/legacy"
UNI="$ART/demo/unified"
TM_LEG="$ART/timing/legacy"
TM_UNI="$ART/timing/unified"

echo "== module self-tests =="
for m in resampler mixer spatial; do
  make -C "$ROOT/$m" clean test >/dev/null
done

echo "== pipeline build (MA_UNIFY=0 and =1) =="
make -C "$ROOT/pipeline" clean all >/dev/null

mkdir -p "$LEG" "$UNI" "$TM_LEG" "$TM_UNI"
"$ROOT/pipeline/pipeline_legacy" demo "$LEG"
"$ROOT/pipeline/pipeline_unified" demo "$UNI"

echo "== timing: sample-rate doubling bench (100 iters, median) =="
"$ROOT/pipeline/pipeline_legacy" bench | tee "$TM_LEG/timing.csv"
"$ROOT/pipeline/pipeline_unified" bench | tee "$TM_UNI/timing.csv"

echo "== independent python verification (-I) =="
PYTHONPATH="$ROOT/tools" python3 -I "$ROOT/tools/analyze.py" "$LEG" "$UNI" "$ART"

PYTHONPATH="$ROOT/tools" python3 -I "$ROOT/tools/make_timing_table.py" \
  "$TM_LEG/timing.csv" "$TM_UNI/timing.csv" "$ART"

echo "all artifacts regenerated in $ART"
