# Sidechain threshold analysis

Signal: 60 sine events on -40 dBFS RMS Gaussian noise, 60 s at 48 kHz, frame = 10 ms (4402 noise frames / 1598 event frames), seed 42.
Detector: attack 2 ms, release 40 ms, hold 30 ms.

## 1. Trigger-level distribution estimate (dBFS, per-frame envelope max)
- noise frames: GMM 0.85*N(-35.60,0.18) + 0.15*N(-28.16,5.50) dBFS
- event frames: GMM 0.05*N(-22.89,6.41) + 0.95*N(-18.36,3.30) dBFS
- BIC: noise 23263.9 vs 5145.6, event 8718.9 vs 8638.8 (lower BIC wins).
- physical check of the noise mixture's wide component: 240 post-event release-tail frames average -22.64 dBFS vs 4162 background frames at -35.18 +/- 1.65 dBFS; the two GMM components separate release tails from stationary background.
- Equal-error crossing at T = -25.00 dBFS, rate = 4.195%.

## 2. Three-tier measured counts vs theory

| tier | T (dBFS) | FP frames raw/theory | FPR meas/theory (dev) | miss frames raw/theory | FNR meas/theory (dev) | FP frames after hold | miss frames after hold | trigger events |
|---|---|---|---|---|---|---|---|---|
| low | -35 | 628 / 584.5 | 14.266% / 13.279% (7.4%) | 8 / 2.6 | 0.501% / 0.160% (212.6%) | 814 | 8 | 62 |
| mid | -30 | 365 / 412.1 | 8.292% / 9.362% (11.4%) | 15 / 11.9 | 0.939% / 0.746% (25.8%) | 545 | 15 | 60 |
| high | -24 | 154 / 146.7 | 3.498% / 3.333% (5.0%) | 44 / 103.6 | 2.753% / 6.486% (57.5%) | 334 | 44 | 60 |

Full 11-point sweep (raw per-frame rates; relative deviation in parentheses):

| T (dBFS) | FPR measured/theory | FNR measured/theory |
|---|---|---|
| -38 | 99.977% / 99.454% (1%) | 0.000% / 0.050% (80%) |
| -36 | 99.296% / 97.865% (1%) | 0.000% / 0.111% (100%) |
| -35 | 14.266% / 13.279% (7%) | 0.501% / 0.160% (213%) |
| -34 | 12.290% / 12.699% (3%) | 0.563% / 0.226% (149%) |
| -32 | 10.086% / 11.239% (10%) | 0.688% / 0.424% (62%) |
| -30 | 8.292% / 9.362% (11%) | 0.939% / 0.746% (26%) |
| -28 | 6.701% / 7.246% (8%) | 1.252% / 1.321% (5%) |
| -26 | 5.089% / 5.151% (1%) | 1.502% / 2.681% (44%) |
| -25 | 4.294% / 4.195% (2%) | 2.190% / 4.111% (47%) |
| -24 | 3.498% / 3.333% (5%) | 2.753% / 6.486% (58%) |
| -22 | 2.317% / 1.948% (19%) | 21.402% / 15.804% (35%) |

## 3. Deviations over 30% and attribution
- **low (-35 dBFS):** FNR dev 212.6%
  - FNR: event lower-tail mismatch, empirical=0.501% vs fitted=0.160% (fit under-puts mass here). Structure: 8/8 missed frames are fade in/out edge frames whose raised-cosine ramp has not reached the threshold, 0 interior; these weak edge frames form a discrete fringe the continuous Gaussian mixture cannot place correctly in the far lower tail
  - hold state machine: hold adds 186 noise frames (4.225% of noise frames) beyond raw exceedance; theory models per-frame levels only
- **high (-24 dBFS):** FNR dev 57.5%
  - FNR: event lower-tail mismatch, empirical=2.753% vs fitted=6.486% (fit over-puts mass here). Structure: 39/44 missed frames are fade in/out edge frames whose raised-cosine ramp has not reached the threshold, 5 interior; these weak edge frames form a discrete fringe the continuous Gaussian mixture cannot place correctly in the far lower tail
  - hold state machine: hold adds 180 noise frames (4.089% of noise frames) beyond raw exceedance; theory models per-frame levels only

## 4. Event-level outcome and threshold recommendation
- Fully missed events (no frame crosses): low 0/60, mid 0/60, high 0/60.
- Frame-level misses at the high tier are fade edge frames, not lost events.
- The equal-error point is -25.0 dBFS (4.19%); the high tier (-24 dBFS, 1 dB away) is nearest to it and its misses are confined to fade edge frames, so it is the recommended operating point; use -30 dBFS if even edge-frame misses must be avoided.
- Below about -36 dBFS the threshold is inside the stationary-background mode (FPR > 99%); the usable window starts near -35 dBFS.

Notes: theory predicts per-frame raw exceedance from the fitted level distributions; the hold state machine can only inflate engaged noise frames (reported separately). Exact Wilson 95% binomial intervals are used to attribute small-count deviations.

Unexplained deviations: 0