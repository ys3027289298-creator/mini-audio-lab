#!/usr/bin/env python3
"""Plot seam_wave.csv / seam_diff.csv / leakage_spec.csv into seam_report.svg.

Pure standard library (no numpy/matplotlib): a minimal SVG line-chart writer.
Run from the looper/ directory: python3 -I tools/plot.py
"""
import csv
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

COLORS = {"raw": "#c0392b", "trim": "#2980b9", "crossfade": "#27ae60"}
LABELS = {"raw": "raw loop (click)", "trim": "integer-period trim",
          "crossfade": "equal-gain crossfade"}


def load(name, xcol, ycol, case_col="case"):
    data = {}
    with open(os.path.join(HERE, name), newline="") as f:
        for row in csv.DictReader(f):
            data.setdefault(row[case_col], []).append(
                (float(row[xcol]), float(row[ycol])))
    for case in data:
        data[case].sort()
    return data


def svg_path(points, x0, y0, w, h, xmin, xmax, ymin, ymax):
    def sx(v):
        return x0 + (v - xmin) / (xmax - xmin) * w

    def sy(v):
        return y0 + h - (v - ymin) / (ymax - ymin) * h

    return " ".join(
        ("M" if i == 0 else "L") + f"{sx(x):.2f},{sy(y):.2f}"
        for i, (x, y) in enumerate(points))


def chart(parts, x0, y0, w, h, title, xlabel, ylabel, log_y=False,
          vline=None, vlabel=""):
    out = []
    xs = [p[0] for case in parts.values() for p in case]
    ys = [p[1] for case in parts.values() for p in case]
    xmin, xmax = min(xs), max(xs)
    if log_y:
        pos = [max(abs(v), 1e-12) for v in ys]
        import math
        ymin, ymax = math.log10(min(pos)), math.log10(max(pos))
        plotted = {c: [(x, math.log10(max(abs(y), 1e-12))) for x, y in pts]
                   for c, pts in parts.items()}
    else:
        ymin, ymax = min(ys), max(ys)
        pad = (ymax - ymin) * 0.08 or 1.0
        ymin -= pad
        ymax += pad
        plotted = parts
    # frame and grid
    out.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" '
               f'fill="white" stroke="#888"/>')
    for g in range(1, 4):
        gy = y0 + h * g / 4
        out.append(f'<line x1="{x0}" y1="{gy:.1f}" x2="{x0+w}" y2="{gy:.1f}" '
                   f'stroke="#eee"/>')
    if vline is not None and xmin <= vline <= xmax:
        vx = x0 + (vline - xmin) / (xmax - xmin) * w
        out.append(f'<line x1="{vx:.1f}" y1="{y0}" x2="{vx:.1f}" '
                   f'y2="{y0+h}" stroke="#333" stroke-dasharray="4,3"/>')
        out.append(f'<text x="{vx+4:.1f}" y="{y0+12}" font-size="11" '
                   f'fill="#333">{vlabel}</text>')
    for case, pts in plotted.items():
        out.append(f'<path d="{svg_path(pts, x0, y0, w, h, xmin, xmax, ymin, ymax)}" '
                   f'fill="none" stroke="{COLORS[case]}" stroke-width="1.5"/>')
    out.append(f'<text x="{x0}" y="{y0-8}" font-size="14" '
               f'font-weight="bold">{title}</text>')
    out.append(f'<text x="{x0+w/2}" y="{y0+h+22}" font-size="12" '
               f'text-anchor="middle">{xlabel}</text>')
    out.append(f'<text x="{x0-38}" y="{y0+h/2}" font-size="12" '
               f'text-anchor="middle" transform="rotate(-90 {x0-38},{y0+h/2})">'
               f'{ylabel}</text>')
    if log_y:
        out.append(f'<text x="{x0+w-6}" y="{y0+12}" font-size="10" '
                   f'text-anchor="end" fill="#666">log10 magnitude</text>')
    return out


def main():
    wave = load("seam_wave.csv", "index", "value")
    diff = load("seam_diff.csv", "index", "diff")
    spec = load("leakage_spec.csv", "fn", "magnitude")

    W, H = 900, 1080
    x0, w = 70, 780
    ch, gap = 300, 60
    parts = [
        chart(wave, x0, 50, w, ch, "1. waveform around the seam (index 1023 -> 1024)",
              "sample index (wrapped)", "amplitude", vline=1023.5,
              vlabel="seam"),
        chart(diff, x0, 50 + ch + gap, w, ch,
              "2. cyclic first difference (a spike here is the click)",
              "sample index (wrapped)", "x[i] - x[i-1]", vline=1023.5,
              vlabel="seam"),
        chart(spec, x0, 50 + 2 * (ch + gap), w, ch,
              "3. DFT magnitude: flat broadband floor = leakage from the seam jump",
              "normalized frequency k/N", "magnitude", log_y=True),
    ]
    legend = []
    lx = x0
    for case in ("raw", "trim", "crossfade"):
        legend.append(
            f'<line x1="{lx}" y1="20" x2="{lx+22}" y2="20" '
            f'stroke="{COLORS[case]}" stroke-width="2.5"/>'
            f'<text x="{lx+27}" y="24" font-size="12">{LABELS[case]}</text>')
        lx += 200
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'viewBox="0 0 {W} {H}" font-family="sans-serif">'
           f'<rect width="{W}" height="{H}" fill="white"/>'
           + "".join(legend)
           + "".join("".join(p) for p in parts) + "</svg>")
    out = os.path.join(HERE, "seam_report.svg")
    with open(out, "w") as f:
        f.write(svg)
    print("wrote", out)


if __name__ == "__main__":
    main()
