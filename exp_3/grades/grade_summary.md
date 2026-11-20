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

- SAM F100.csv: exp_3/sam_outputs/F100.csv: first line is '# F100: Mdot=0.0001 Msun/yr L=1e+38 erg/s nH=100.0 t=100.0 kyr; R1=1.7289 Rc=6.5', expected 'r,rho,v,P,T'
- SAM F10.csv: exp_3/sam_outputs/F10.csv: first line is '# F10: Mdot=0.0001 Msun/yr L=1e+38 erg/s nH=10.0 t=100.0 kyr; R1=3.4496 Rc=10.38', expected 'r,rho,v,P,T'
- SAM F1.csv: exp_3/sam_outputs/F1.csv: first line is '# F1: Mdot=0.0001 Msun/yr L=1e+38 erg/s nH=1.0 t=100.0 kyr; R1=6.8828 Rc=16.451 ', expected 'r,rho,v,P,T'
- sim N=1024 times: exp_3/sim_outputs/F100_sim_N01024_times.txt: line 1 is not numeric: '# snapshot  t_s  t_kyr'
- sim N=2048 times: exp_3/sim_outputs/F100_sim_N02048_times.txt: line 1 is not numeric: '# snapshot  t_s  t_kyr'
- sim N=4096 times: exp_3/sim_outputs/F100_sim_N04096_times.txt: line 1 is not numeric: '# snapshot  t_s  t_kyr'
- sim N=8192 times: exp_3/sim_outputs/F100_sim_N08192_times.txt: line 1 is not numeric: '# snapshot  t_s  t_kyr'
- convergence temperature: [Errno 2] No such file or directory: 'exp_3/convergence/temperature_L1_norm.txt'
- convergence density: [Errno 2] No such file or directory: 'exp_3/convergence/density_L1_norm.txt'
- convergence velocity: [Errno 2] No such file or directory: 'exp_3/convergence/velocity_L1_norm.txt'
- convergence pressure: [Errno 2] No such file or directory: 'exp_3/convergence/pressure_L1_norm.txt'
