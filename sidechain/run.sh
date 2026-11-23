#!/usr/bin/env bash
# Reproduce the full sidechain threshold analysis end to end.
set -euo pipefail
cd "$(dirname "$0")"
PY=../.venv/bin/python

make test
"$PY" -I analysis/gen_signal.py data
"$PY" -I analysis/analyze.py data results
echo "run.sh: OK"
