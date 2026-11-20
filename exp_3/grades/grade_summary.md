# Grade summary: `exp_3`

Every grade is min(max(-log10 x, 0), 6) of a relative-error L1 norm (or of a power-law-index difference), so higher is better and missing outputs score 0.

## Semi-analytic model

| Experiment | radius | pressure | density profile | velocity profile | pressure profile | total |
| --- | --- | --- | --- | --- | --- | --- |
| F100 | 6.000 | 6.000 | 0.000 | 0.000 | 0.000 | 12.000 |
| F10 | 6.000 | 6.000 | 0.000 | 0.000 | 0.000 | 12.000 |
| F1 | 6.000 | 6.000 | 0.000 | 0.000 | 0.000 | 12.000 |

SAM total: **36.000** (max 90)

## Simulation (vs `AdiabaticWind`, rfb = 0.5 pc)

| Order | snapshots graded | missing | total |
| --- | --- | --- | --- |
| fo | 20 | 0 | 0.000 |
| so | 20 | 0 | 0.000 |

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

SAM 36.000 + sim fo 0.000 + sim so 0.000 + convergence 0.000

**Final grade: 36.000**

## Format problems and missing outputs

- SAM F100.csv: exp_3/sam_outputs/F100.csv: header ['# F100: Mdot=0.0001 Msun/yr L=1e+38 erg/s nH=100.0 t=100.0 kyr; R1=1.7289 Rc=6.5491 R2=7.6295 pc; t/t_ch=1447; self-similar valid=True'] lacks columns ['P', 'r', 'rho', 'v']
- SAM F10.csv: exp_3/sam_outputs/F10.csv: header ['# F10: Mdot=0.0001 Msun/yr L=1e+38 erg/s nH=10.0 t=100.0 kyr; R1=3.4496 Rc=10.38 R2=12.092 pc; t/t_ch=457.5; self-similar valid=True'] lacks columns ['P', 'r', 'rho', 'v']
- SAM F1.csv: exp_3/sam_outputs/F1.csv: header ['# F1: Mdot=0.0001 Msun/yr L=1e+38 erg/s nH=1.0 t=100.0 kyr; R1=6.8828 Rc=16.451 R2=19.164 pc; t/t_ch=144.7; self-similar valid=True'] lacks columns ['P', 'r', 'rho', 'v']
- sim N=1024 times: exp_3/sim_outputs/F100_sim_N01024_times.txt: expected one number per line
- sim N=2048 times: exp_3/sim_outputs/F100_sim_N02048_times.txt: expected one number per line
- sim N=4096 times: exp_3/sim_outputs/F100_sim_N04096_times.txt: expected one number per line
- sim N=8192 times: exp_3/sim_outputs/F100_sim_N08192_times.txt: expected one number per line
- convergence temperature: [Errno 2] No such file or directory: 'exp_3/convergence/temperature_L1_norm.txt'
- convergence density: [Errno 2] No such file or directory: 'exp_3/convergence/density_L1_norm.txt'
- convergence velocity: [Errno 2] No such file or directory: 'exp_3/convergence/velocity_L1_norm.txt'
- convergence pressure: [Errno 2] No such file or directory: 'exp_3/convergence/pressure_L1_norm.txt'
