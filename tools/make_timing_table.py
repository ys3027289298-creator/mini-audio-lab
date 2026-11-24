#!/usr/bin/env python3
"""Merge the two bench CSVs into timing_table.md / timing_table.csv."""
import csv
import os
import sys


def read(path):
    rows = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows[(r["stage"], r["sr_in"], r["sr_out"])] = r
    return rows


def main():
    leg = read(sys.argv[1])
    uni = read(sys.argv[2])
    out = sys.argv[3]
    os.makedirs(out, exist_ok=True)

    order = [("mixer", "8000", "16000"), ("mixer", "16000", "32000"),
             ("resampler", "8000", "16000"), ("resampler", "16000", "32000"),
             ("spatial", "8000", "16000"), ("spatial", "16000", "32000")]

    with open(os.path.join(out, "timing_table.md"), "w") as md, \
         open(os.path.join(out, "timing_table.csv"), "w", newline="") as cf:
        md.write("# Stage timing, sample-rate doubling (ratio=0.5), "
                 "12 s workload, 100 iterations, median\n\n")
        md.write("| stage | sr_in | sr_out | legacy median (ms) | "
                 "unified median (ms) | delta (ms) | "
                 "unified throughput (Msps) |\n")
        md.write("|---|---|---|---|---|---|---|\n")
        w = csv.writer(cf)
        w.writerow(["stage", "sr_in", "sr_out", "legacy median (ms)",
                    "unified median (ms)", "delta (ms)",
                    "unified throughput (Msps)"])
        for key in order:
            l, u = leg[key], uni[key]
            lm, um = float(l["median_ms"]), float(u["median_ms"])
            msps = float(u["throughput_Msps"])
            md.write(f"| {key[0]} | {key[1]} | {key[2]} | {lm:.4f} | "
                     f"{um:.4f} | {um - lm:+.4f} | {msps:.1f} |\n")
            w.writerow([key[0], key[1], key[2], f"{lm:.4f}", f"{um:.4f}",
                        f"{um - lm:+.4f}", f"{msps:.1f}"])


if __name__ == "__main__":
    main()
