# Grade summary: `opus_low`

Every grade is min(max(-log10 x, 0), 6) of a relative-error L1 norm (or of a power-law-index difference), so higher is better and missing outputs score 0.

## Semi-analytic model

| Experiment | radius | pressure | density profile | velocity profile | pressure profile | total |
| --- | --- | --- | --- | --- | --- | --- |
| F100 | 0.864 | 1.771 | 0.000 | 0.109 | 0.000 | 2.744 |
| F10 | 0.864 | 1.771 | 0.000 | 0.324 | 0.000 | 2.959 |
| F1 | 0.864 | 1.771 | 0.000 | 0.524 | 0.000 | 3.159 |

SAM total: **8.861** (max 90)

## Simulation (vs `AdiabaticWind`, rfb = 0.5 pc)

| Order | snapshots graded | missing | total |
| --- | --- | --- | --- |
| fo | 5 | 10 | 3.526 |
| so | 1 | 14 | 0.000 |

Per-snapshot grades are in `sim_grade_fo.txt` and `sim_grade_so.txt`.

## Convergence

| Order | temperature | density | velocity | pressure |
| --- | --- | --- | --- | --- |
| fo | 0.000 (missing) | 0.000 (missing) | 0.000 (missing) | 0.000 (missing) |
| so | 0.000 (missing) | 0.000 (missing) | 0.000 (missing) | 0.000 (missing) |

Convergence total: **0.000**

## Code quality

Not graded (no `grades/code_quality.txt`). See README.md for how to add it.

## Final grade

SAM 8.861 + sim fo 3.526 + sim so 0.000 + convergence 0.000

**Final grade: 12.388**

## Format problems and missing outputs

- sim fo: 2 of 3 resolutions missing
- sim so N=64: 4 of 5 snapshots missing
- sim so: 2 of 3 resolutions missing
- convergence temperature: [Errno 2] No such file or directory: 'opus_low/convergence/temperature_L1_norm.txt'
- convergence density: [Errno 2] No such file or directory: 'opus_low/convergence/density_L1_norm.txt'
- convergence velocity: [Errno 2] No such file or directory: 'opus_low/convergence/velocity_L1_norm.txt'
- convergence pressure: [Errno 2] No such file or directory: 'opus_low/convergence/pressure_L1_norm.txt'
