# paper7 — numerics provenance / fixtures

Reproducibility map for the quoted numbers in `paper7/main.tex`, so an audit does
not re-stall. Values below were reproduced by fresh audit agents (first gate audit
2026-07-14 and the full-wave numerics trace) running the committed scripts; the
large-M folds, the 0.8620 elimination, and the MPS floors already have committed
output files (`results_val_M*.txt`, `results_partial_folds_M8_cap6.txt`,
`mps_*_result*.txt`, `dt01_*.txt`, `diag_*_c64_dt02.txt`) and are not duplicated
here. **Reproduce with:** `python scripts/<name>.py`.

## Model + census — `t5_channel_ladder.py`, `t5_dimension_ladder.py`
| claim (main.tex) | paper value | reproduced |
|---|---|---|
| resonant channel census, M=2/3/4/5 | 10/21/36/55 | **10/21/36/55** (exact combinatorics) |
| analytic channel table vs exact vacuum column | 5.7e-16 | **5.68e-16** (max rel err) |

## Schur fold, small M — `t5_stage3.py`
| config | full/naive | full/Schur |
|---|---|---|
| M=2/3/4, cap 4, λ=0.02, t=40 | --- | **1.0000 / 0.9997 / 0.9993** |
| M=2/3/4, cap 6 | --- | **0.9998 / 0.9984 / 0.9963** |
| M=4, cap 4, λ=0.04 | **0.9434** | **0.9972** (orchestrator-verified 2026-07-14; the paper's 0.8964 was a wrong-cell copy, fixed) |
| λ-deficit series (M=3, cap 6), λ=0.005→0.04 | 0.0057/0.0227/0.0883/0.3144 | **0.0057/0.0227/0.0883/0.3144** (×3.98/3.89/3.56 per doubling) |
| cap convergence M=2 (caps 4/6/8/10): Schur flat / naive→ / dim | 0.9998–1.0000 / ~0.956 / 2568 | reproduced |

## Attribution + rate window — `t5_stage4.py` (M=8)
| claim | paper value | reproduced |
|---|---|---|
| rate window Γ(t)=a+bt: a / b / 2S_exact | 2.03e-4 / 4.28e-5 / 4.53e-5 | **2.032e-4 / 4.275e-5 / 4.528e-5** (b vs 2S = 5.6% ≈ "6%") |
| golden-rule plateau bracket | 1.56–1.84e-4 | reproduced |
| born-vacuum sustained-then-collapse, t=5→120 | 471→7 | **471→7** |
| d=2 spot check (48-channel shell, cap 4), full/Schur | 1.0000 | **1.0000** |

## Dimension ladder — `t5_dimension_ladder.py`
| prediction | paper value | reproduced |
|---|---|---|
| P1 d=1 resonant-shell saturation (last doubling, Λ=100) | +1.5% | **+1.51%** |
| P2 d=1 golden-rule sum, power-fit exponent (log-class) | 0.156 | **0.156** |
| P3 d=2 sum, fitted exponent (one shell deep) | 2.26 / 2.263 | **2.263** |
| channel-table extensivity under volume doubling | ~7.4% | **7.4%** |

## MPS L=32 diagnostic (§mps) — RESOLVED, confirms the paper
`dt01_L32.txt` (L=32, dt=0.01, χ=40, sign=−1, λ=0.06; walltime ~5.7 h) gives
**dE_end_per_site = −4.12e−3** and rate_per_site = −9.33e−4 — the per-site transfer
is genuinely NEGATIVE, confirming the §mps sentence "an L=32 point whose per-site
transfer came out negative, i.e. non-extensive under truncation." (The L=16/dt=0.04
floor diagnostic, −4.21e−2, is a distinct case.) completion_candidates-style
production runs are not part of paper7.
