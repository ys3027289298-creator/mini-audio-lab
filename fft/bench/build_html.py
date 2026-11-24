#!/usr/bin/env python3
"""Generate the self-contained interactive trade-off HTML from the CSVs."""
import csv, json, pathlib

here = pathlib.Path(__file__).parent


def read_csv(name):
    with open(here / name) as f:
        return list(csv.DictReader(f))


def numify(rows, ints=()):
    out = []
    for r in rows:
        d = {}
        for k, v in r.items():
            if k in ints:
                d[k] = int(v)
            else:
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
        out.append(d)
    return out


results = numify(read_csv("results_blocking1.csv"), ints=("rate_hz", "block", "navg", "resolved"))
candidates = numify(read_csv("candidates.csv"), ints=("rate_hz", "block"))
timing = numify(read_csv("timing.csv"), ints=("rate_hz", "block"))
fits = (here / "scaling_fit.txt").read_text()

payload = {"results": results, "candidates": candidates, "timing": timing, "fits": fits}
DATA_JSON = json.dumps(payload)

html = """<!doctype html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>分块 FFT 频谱 · 误差-延迟权衡</title>
<style>
  .viz-root {
    color-scheme: light;
    --surface-1:#fcfcfb; --page:#f9f9f7; --ink:#0b0b0b; --ink-2:#52514e; --muted:#898781;
    --grid:#e1e0d9; --axis:#c3c2b7; --border:rgba(11,11,11,.10);
    --blue:#2a78d6; --orange:#eb6834; --aqua:#1baf7a; --violet:#6250d6; --red:#e34948;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) .viz-root {
      color-scheme: dark;
      --surface-1:#1a1a19; --page:#0d0d0d; --ink:#f0efec; --ink-2:#c3c2b7; --muted:#898781;
      --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,.10);
      --blue:#3987e5; --orange:#d95926; --aqua:#199e70; --violet:#9085e9; --red:#e66767;
    }
  }
  :root[data-theme="dark"] .viz-root {
    color-scheme: dark;
    --surface-1:#1a1a19; --page:#0d0d0d; --ink:#f0efec; --ink-2:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,.10);
    --blue:#3987e5; --orange:#d95926; --aqua:#199e70; --violet:#9085e9; --red:#e66767;
  }
  * { box-sizing: border-box; }
  body.viz-root { margin:0; background:var(--page); color:var(--ink);
    font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif; padding:28px; }
  .wrap { max-width:1180px; margin:0 auto; }
  h1 { font-size:22px; margin:0 0 4px; }
  .sub { color:var(--ink-2); margin:0 0 20px; }
  .controls { display:flex; gap:8px; align-items:center; margin-bottom:18px; flex-wrap:wrap; }
  .controls button { font:inherit; padding:6px 14px; border-radius:8px; border:1px solid var(--border);
    background:var(--surface-1); color:var(--ink-2); cursor:pointer; }
  .controls button[aria-pressed="true"] { background:var(--blue); color:#fff; border-color:transparent; }
  .spacer { flex:1; }
  .grid3 { display:grid; grid-template-columns:repeat(3,1fr); gap:16px; }
  @media (max-width:980px){ .grid3{grid-template-columns:1fr;} }
  .card { background:var(--surface-1); border:1px solid var(--border); border-radius:12px; padding:14px 16px; }
  .card h2 { font-size:14px; margin:0 0 2px; color:var(--ink); }
  .card .cap { font-size:12px; color:var(--muted); margin:0 0 8px; }
  svg { width:100%; height:auto; display:block; }
  .axis { stroke:var(--axis); stroke-width:1; }
  .gl { stroke:var(--grid); stroke-width:1; }
  .tick { fill:var(--muted); font-size:10px; }
  .axis-label { fill:var(--ink-2); font-size:11px; }
  .pt { cursor:pointer; }
  .ring { fill:none; stroke-width:2.5; }
  table { width:100%; border-collapse:collapse; font-size:13px; margin-top:6px; }
  th,td { text-align:right; padding:7px 8px; border-bottom:1px solid var(--grid);
    font-variant-numeric:tabular-nums; white-space:nowrap; }
  th:first-child,td:first-child, th:nth-child(2),td:nth-child(2){ text-align:left; }
  th { color:var(--muted); font-weight:600; }
  td { color:var(--ink); }
  .swatch { display:inline-block; width:10px; height:10px; border-radius:3px; margin-right:6px; vertical-align:middle; }
  .fit { white-space:pre-wrap; font-family:ui-monospace,Menlo,Consolas,monospace; font-size:12px;
    color:var(--ink-2); background:var(--surface-1); border:1px solid var(--border);
    border-radius:10px; padding:12px 14px; overflow:auto; }
  #tip { position:fixed; pointer-events:none; background:var(--surface-1); border:1px solid var(--border);
    border-radius:8px; padding:8px 10px; font-size:12px; color:var(--ink);
    box-shadow:0 6px 20px rgba(0,0,0,.18); opacity:0; transition:opacity .1s; z-index:10; }
  #tip b { font-variant-numeric:tabular-nums; }
  .legend { display:flex; gap:16px; font-size:12px; color:var(--ink-2); margin:4px 0 0; flex-wrap:wrap; }
  .yes { color:var(--aqua); font-weight:600; } .no { color:var(--red); font-weight:600; }
</style>
</head>
<body class="viz-root">
<div class="wrap">
  <h1>分块 FFT 频谱：误差–延迟权衡</h1>
  <p class="sub">Hann 窗 + 50% 重叠 Welch 分块估计，对照整段 FFT 参考（density 约定，功率守恒，与 scipy.signal.welch 一致）。</p>
  <div class="controls">
    <button data-rate="48000" aria-pressed="true">48 kHz</button>
    <button data-rate="96000" aria-pressed="false">96 kHz</button>
    <span class="spacer"></span>
    <button id="theme" aria-pressed="false">切换深色</button>
  </div>
  <div class="grid3">
    <div class="card"><h2>误差下降律</h2><p class="cap">分块误差 E(dB) 随块长 B 下降</p><div id="c_err"></div></div>
    <div class="card"><h2>延迟上升律</h2><p class="cap">算法延迟 L = B/fs 随 B 线性上升</p><div id="c_lat"></div></div>
    <div class="card"><h2>权衡曲线</h2><p class="cap">误差 vs 延迟，左下更优；圈=三个候选块长</p><div id="c_trade"></div></div>
  </div>
  <div class="legend">
    <span><span class="swatch" style="background:var(--orange)"></span>误差 (dB)</span>
    <span><span class="swatch" style="background:var(--blue)"></span>延迟 (ms)</span>
    <span><span class="swatch ring" style="border:2.5px solid var(--aqua)"></span>候选 B=128</span>
    <span><span class="swatch ring" style="border:2.5px solid var(--red)"></span>候选 B=512</span>
    <span><span class="swatch ring" style="border:2.5px solid var(--orange)"></span>候选 B=1024</span>
  </div>

  <h2 style="font-size:15px;margin:26px 0 6px">三个候选块长评估</h2>
  <div class="card" style="overflow:auto"><table id="cand"></table></div>
  <h2 style="font-size:15px;margin:22px 0 6px">两档采样率耗时（单次 FFT 与每音频秒 CPU）</h2>
  <div class="card" style="overflow:auto"><table id="time"></table></div>
  <h2 style="font-size:15px;margin:22px 0 6px">拟合的两条互制约标度律</h2>
  <div class="fit" id="fit"></div>
</div>
<div id="tip"></div>

<script>
const DATA = __DATA__;
const CAND = new Set([128, 512, 1024]);
const CAND_COLOR = {128:"var(--aqua)", 512:"var(--red)", 1024:"var(--orange)"};
let rate = 48000;
const tip = document.getElementById("tip");
const css = n => getComputedStyle(document.body).getPropertyValue(n).trim();

function rows(){ return DATA.results.filter(r => r.rate_hz===rate && r.mode==="block")
  .sort((a,b)=>a.block-b.block); }

function logScale(dmin,dmax,r0,r1){ const l0=Math.log10(dmin),l1=Math.log10(dmax);
  return v => r0 + (Math.log10(v)-l0)/(l1-l0)*(r1-r0); }
function linScale(dmin,dmax,r0,r1){ return v => r0+(v-dmin)/(dmax-dmin)*(r1-r0); }

function chart(mount, pts, xacc, yacc, opts){
  const W=360,H=250, M={t:14,r:16,b:42,l:52};
  const iw=W-M.l-M.r, ih=H-M.t-M.b;
  const xs=pts.map(xacc), ys=pts.map(yacc);
  const xmin=Math.min(...xs), xmax=Math.max(...xs), ymin=Math.min(...ys), ymax=Math.max(...ys);
  const sx = opts.xlog ? logScale(xmin*0.9, xmax*1.1, M.l, M.l+iw)
                        : linScale(xmin, xmax, M.l, M.l+iw);
  const sy = linScale(Math.min(0,ymin), ymax*1.08, M.t+ih, M.t);
  let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${opts.aria}">`;
  const xticks = opts.xlog ? [64,128,256,512,1024,2048].filter(v=>v>=xmin&&v<=xmax)
                           : xs.filter((_,i)=>i%2===0);
  xticks.forEach(v=>{ const x=sx(v);
    s+=`<line class="gl" x1="${x}" y1="${M.t}" x2="${x}" y2="${M.t+ih}"/>`;
    s+=`<text class="tick" x="${x}" y="${M.t+ih+15}" text-anchor="middle">${v}</text>`; });
  for(let i=0;i<=4;i++){ const v=ymin+(ymax-ymin)*i/4, y=sy(v);
    s+=`<line class="gl" x1="${M.l}" y1="${y}" x2="${M.l+iw}" y2="${y}"/>`;
    s+=`<text class="tick" x="${M.l-7}" y="${y+3.5}" text-anchor="end">${v.toFixed(opts.ydec??1)}</text>`; }
  s+=`<line class="axis" x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ih}"/>`;
  s+=`<line class="axis" x1="${M.l}" y1="${M.t+ih}" x2="${M.l+iw}" y2="${M.t+ih}"/>`;
  const path = pts.map((p,i)=>`${i?"L":"M"}${sx(xacc(p)).toFixed(1)},${sy(yacc(p)).toFixed(1)}`).join("");
  s+=`<path d="${path}" fill="none" stroke="${opts.color}" stroke-width="2"/>`;
  pts.forEach(p=>{ const x=sx(xacc(p)),y=sy(yacc(p));
    const ring = CAND.has(p.block) ? `<circle class="ring" cx="${x}" cy="${y}" r="9" stroke="${CAND_COLOR[p.block]}"/>` : "";
    s+=`${ring}<circle class="pt" cx="${x}" cy="${y}" r="6" fill="${opts.color}" data-tip='${opts.tip(p)}'/>`; });
  s+=`<text class="axis-label" x="${M.l+iw/2}" y="${H-6}" text-anchor="middle">${opts.xlabel}</text>`;
  s+=`<text class="axis-label" transform="translate(14,${M.t+ih/2}) rotate(-90)" text-anchor="middle">${opts.ylabel}</text>`;
  s+=`</svg>`;
  mount.innerHTML=s;
  mount.querySelectorAll(".pt").forEach(el=>{
    el.addEventListener("mousemove",e=>{ tip.innerHTML=el.dataset.tip; tip.style.opacity=1;
      tip.style.left=(e.clientX+14)+"px"; tip.style.top=(e.clientY+14)+"px"; });
    el.addEventListener("mouseleave",()=>tip.style.opacity=0); });
}

function render(){
  const p = rows();
  const tipE = d=>`B=<b>${d.block}</b><br>延迟 <b>${d.latency_ms.toFixed(2)}</b> ms<br>误差 <b>${d.band_err_db.toFixed(2)}</b> dB<br>峰值误差 <b>${d.peak_err_hz.toFixed(1)}</b> Hz<br>双音分辨 <b>${d.resolved?"是":"否"}</b>`;
  chart(document.getElementById("c_err"), p, d=>d.block, d=>d.band_err_db,
    {color:"var(--orange)", xlog:true, xlabel:"块长 B（样本）", ylabel:"误差 (dB)", aria:"误差下降律", tip:tipE});
  chart(document.getElementById("c_lat"), p, d=>d.block, d=>d.latency_ms,
    {color:"var(--blue)", xlog:true, xlabel:"块长 B（样本）", ylabel:"延迟 (ms)", aria:"延迟上升律", tip:tipE, ydec:0});
  chart(document.getElementById("c_trade"), p, d=>d.latency_ms, d=>d.band_err_db,
    {color:"var(--violet)", xlog:true, xlabel:"延迟 (ms)", ylabel:"误差 (dB)", aria:"权衡曲线", tip:tipE});

  const cand = DATA.candidates.filter(r=>r.rate_hz===rate);
  document.getElementById("cand").innerHTML =
    `<tr><th>块长 B</th><th>延迟 (ms)</th><th>误差 (dB)</th><th>峰值误差 (Hz)</th><th>双音分辨</th><th>单次 FFT (µs)</th><th>整段处理 (ms)</th></tr>` +
    cand.map(r=>`<tr><td>${r.block}</td><td>${r.latency_ms.toFixed(2)}</td><td>${r.band_err_db.toFixed(2)}</td>
      <td>${r.peak_err_hz.toFixed(2)}</td><td class="${r.pair_resolved==="yes"?"yes":"no"}">${r.pair_resolved==="yes"?"是":"否"}</td>
      <td>${r.fft_us.toFixed(2)}</td><td>${r.total_ms.toFixed(3)}</td></tr>`).join("");

  const tm = DATA.timing;
  document.getElementById("time").innerHTML =
    `<tr><th>采样率</th><th>块长 B</th><th>延迟 (ms)</th><th>单次 FFT (µs)</th><th>每音频秒 FFT (µs)</th><th>CPU 占比</th></tr>` +
    tm.map(r=>`<tr><td>${r.rate_hz} Hz</td><td>${r.block}</td><td>${r.latency_ms.toFixed(2)}</td>
      <td>${r.fft_us_per_transform.toFixed(2)}</td><td>${r.fft_us_per_audio_second.toFixed(0)}</td>
      <td>${(r.cpu_fraction*100).toFixed(3)}%</td></tr>`).join("");

  document.getElementById("fit").textContent =
    DATA.fits.split("=== ").filter(s=>s.startsWith(String(rate))).map(s=>"=== "+s).join("").trim();
}

document.querySelectorAll("[data-rate]").forEach(b=>b.addEventListener("click",()=>{
  rate=+b.dataset.rate;
  document.querySelectorAll("[data-rate]").forEach(x=>x.setAttribute("aria-pressed",x===b));
  render();
}));
document.getElementById("theme").addEventListener("click",()=>{
  const root=document.documentElement;
  const dark = root.getAttribute("data-theme")!=="dark";
  root.setAttribute("data-theme", dark?"dark":"light");
  document.getElementById("theme").setAttribute("aria-pressed", dark);
  document.getElementById("theme").textContent = dark?"切换浅色":"切换深色";
});
render();
</script>
</body>
</html>
"""
(here / "tradeoff.html").write_text(html.replace("__DATA__", DATA_JSON))
print("wrote bench/tradeoff.html")
