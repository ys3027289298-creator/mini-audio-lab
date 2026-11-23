#!/usr/bin/env python3
"""Render the measured anti-aliasing frequency response (CSV -> SVG).

Stdlib only.  Usage: plot_freq.py doc/freq_response.csv doc/freq_response.svg
Each ratio in the CSV gets one panel: red = filter off (aliased),
blue = filter on, dashed markers = passband/stopband edges.
"""
import csv
import sys

Y0, Y1 = -100.0, 5.0  # dB range


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def main(csv_path, svg_path):
    groups = {}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            groups.setdefault(float(row["ratio"]), []).append(
                (float(row["freq"]), float(row["db_on"]), float(row["db_off"])))

    W, PH = 780, 330
    ML, MR, MT, MB = 70, 24, 46, 54
    GW, GH = W - ML - MR, PH - MT - MB
    ratios = sorted(groups)
    H = PH * len(ratios)

    def xmap(f):
        return ML + GW * f / 0.5

    def ymap(db):
        return MT + GH * (Y1 - db) / (Y1 - Y0)

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"'
        f' font-family="sans-serif" font-size="12">',
        f'<rect width="{W}" height="{H}" fill="white"/>',
    ]
    for pi, ratio in enumerate(ratios):
        yoff = pi * PH

        def X(f):
            return xmap(f)

        def Y(db):
            return yoff + ymap(db)

        out.append(f'<rect x="{ML}" y="{yoff + MT}" width="{GW}" height="{GH}"'
                   f' fill="#fafafa" stroke="#333"/>')
        f = 0.0
        while f <= 0.5001:
            x = X(f)
            out.append(f'<line x1="{x:.1f}" y1="{yoff + MT}" x2="{x:.1f}"'
                       f' y2="{yoff + MT + GH}" stroke="#e2e2e2"/>')
            out.append(f'<text x="{x:.1f}" y="{yoff + MT + GH + 16}"'
                       f' text-anchor="middle" fill="#333">{f:.2f}</text>')
            f += 0.05
        db = Y0
        while db <= Y1:
            y = Y(db)
            out.append(f'<line x1="{ML}" y1="{y:.1f}" x2="{ML + GW}" y2="{y:.1f}"'
                       f' stroke="#e2e2e2"/>')
            out.append(f'<text x="{ML - 6}" y="{y + 4:.1f}" text-anchor="end"'
                       f' fill="#333">{db:.0f}</text>')
            db += 20

        fp_, fs_ = 0.4 / ratio, 0.5 / ratio
        for fe, col in ((fp_, "#2a9d3a"), (fs_, "#888888")):
            out.append(f'<line x1="{X(fe):.1f}" y1="{yoff + MT}" x2="{X(fe):.1f}"'
                       f' y2="{yoff + MT + GH}" stroke="{col}"'
                       f' stroke-dasharray="5,4" stroke-width="1.4"/>')

        data = groups[ratio]
        for col, key in (("#c0392b", 2), ("#1f6fb2", 1)):  # off, on
            pts = " ".join(f"{X(d[0]):.1f},{Y(max(d[key], Y0)):.1f}" for d in data)
            out.append(f'<polyline points="{pts}" fill="none" stroke="{col}"'
                       f' stroke-width="1.6"/>')

        # title, legend, axis labels
        out.append(f'<text x="{ML}" y="{yoff + 22}" font-size="14" fill="#111">'
                   f'降采样比率 ratio = {ratio:g} — 抗混叠频响实测</text>')
        lx = ML + GW - 300
        ly = yoff + MT + 18
        for i, (col, label, dash) in enumerate((
                ("#1f6fb2", "开抗混叠", ""),
                ("#c0392b", "关闭（回卷）", ""),
                ("#2a9d3a", "通带边 f_p", ' stroke-dasharray="5,4"'),
                ("#888888", "阻带边 f_s = 新奈奎斯特", ' stroke-dasharray="5,4"'))):
            yy = ly + i * 18
            out.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 26}" y2="{yy}"'
                       f' stroke="{col}" stroke-width="1.8"{dash}/>')
            out.append(f'<text x="{lx + 32}" y="{yy + 4}" fill="#333">'
                       f'{esc(label)}</text>')
        out.append(f'<text x="{ML + GW / 2:.0f}" y="{yoff + PH - 8}"'
                   f' text-anchor="middle" fill="#111">归一化频率（周期 / 输入采样）</text>')
        out.append(f'<text x="16" y="{yoff + MT + GH / 2:.0f}" text-anchor="middle"'
                   f' fill="#111" transform="rotate(-90 16 {yoff + MT + GH / 2:.0f})">'
                   f'幅度 (dB)</text>')

    out.append("</svg>")
    with open(svg_path, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {svg_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: plot_freq.py <csv> <svg>")
    main(sys.argv[1], sys.argv[2])
