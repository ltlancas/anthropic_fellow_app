# Code documentation

All code lives in `exp_3/` and runs in a standalone virtual environment `exp_3/.venv`
(Python 3.9, numpy, scipy, numba, pytest). Run all commands from `exp_3/`.

## Layout

```
sam/                 Task 1 semi-analytic model (frozen after Task 1)
  units.py           physical constants + unit conversions (cgs internally)
  model.py           WindBubble class: R2(t), Rc(t), R1(t), P_bubble(t), profile(r, t)
hydro/               Task 2 hydrodynamics code
  solver.py          Grid, Simulation, numba kernels (HLLC, reconstruction, time stepping)
scripts/
  run_sam.py         Task 1 driver -> sam_outputs/
  run_F100.py        Task 2 driver -> sim_outputs/
  convergence.py     Task 3 driver -> convergence/
tests/test_hydro.py  verification tests (<10 s; ~10 s extra the first time for numba compilation)
sam_outputs/ sim_outputs/ convergence/   results
```

## Hydro code (`hydro/solver.py`)

**Equations.** 1D Euler equations in spherical symmetry, conservative finite-volume form for
U = (ρ, ρv, E), with E = P/(γ−1) + ρv²/2 and γ = 5/3 by default:

d(U_i)/dt = −(A_{i+½} F_{i+½} − A_{i−½} F_{i−½})/V_i + S_i,

with A = 4πr², V_i = (4π/3)(r_{i+½}³ − r_{i−½}³). The geometric momentum source is
P_i (A_{i+½} − A_{i−½})/V_i. This form keeps a uniform static gas exactly static (it is
well-balanced). `geometry="cartesian"` (A = 1) exists for testing.

**Fluxes.** HLLC approximate Riemann solver (Toro), with Davis wave-speed estimates.

**Orders.**
- `order=1`: piecewise-constant states + forward Euler (first order in space and time).
- `order=2`: piecewise-linear reconstruction of primitive variables (ρ, v, P) with the
  monotonised-central (MC) limiter, and 2-stage SSP Runge–Kutta (Heun) in time. Second order in
  space and time for smooth flow. If a reconstructed ρ or P would be ≤ 0, that cell falls back
  to first order.

**Boundaries.** At r = 0, reflecting ghost cells (v → −v); the face area there is zero anyway.
At r_max, zero-gradient outflow.

**Sources.** Mass Ṁ and energy L are injected at rest, uniformly per unit volume, inside R_inj.
Each cell gets a share equal to its volume overlap with the sphere r < R_inj, so the totals are
exact even when R_inj is smaller than one cell.

**Time step.** dt = CFL · Δr / max(|v| + c) (CFL = 0.4 by default), further limited so that the
mass or energy injected in one step is at most 30% of a cell's content (`SRC_FRAC`). Steps are
cut so that output times are hit exactly.

**Floors.** ρ, P ≥ 1e−40 (cgs); these never activated in the tests.

### API
```python
from hydro.solver import Grid, Simulation
grid = Grid(N, r_max, geometry="spherical")       # uniform cells on [0, r_max]
sim  = Simulation(grid, W0, gamma=5/3, order=2, cfl=0.4, mdot=..., L=..., R_inj=...)
sim.run_to(t)                                     # advance to time t (s)
W = sim.primitives()                              # (3, N): rho, v, P ; radii grid.rc
```
`grid.rc` is the volume-weighted cell centre.

## Running

```bash
.venv/bin/python -m pytest tests -q                    # hydro tests
.venv/bin/python scripts/run_sam.py                    # Task 1
.venv/bin/python scripts/run_F100.py --N 1024 --order 1  # Task 2 (one run)
for N in 1024 2048 4096 8192; do for o in 1 2; do
  .venv/bin/python scripts/run_F100.py --N $N --order $o & done; done; wait
.venv/bin/python scripts/convergence.py                # Task 3
```

## F100 set-up (`scripts/run_F100.py`)
- Domain [0, 1.5 R_2,SAM(100 kyr)] = [0, 11.44 pc]. With this, R_2 sits at 2/3 of the domain at
  the final time and the outer boundary is never reached.
- Initial state: uniform ambient gas, n_H = 100 cm⁻³, ρ = 1.4 m_p n_H, T = 10 K, v = 0.
  μ = 1.4 is used in T = μ m_p P/(ρ k_B), the same as in the SAM.
- Ṁ = 1e−4 M☉/yr and L = 1e38 erg/s are injected within 100 AU. At these resolutions
  (Δr ≥ 290 AU) the injection region is smaller than one cell, so everything is deposited
  in cell 0.
- 5 snapshots at t = 20, 40, 60, 80, 100 kyr.
- N = 1024, 2048, 4096, 8192 cells; first (`fo`) and second (`so`) order.
- Output columns: r_cm, r_pc, rho_gcm3, nH_cm3, v_cms, P_dyncm2, T_K.

## Tests (`tests/test_hydro.py`)
1. Sod shock tube (Cartesian, γ = 1.4) against the exact Riemann solution. The error falls
   with resolution, and second order beats first order.
2. A uniform static sphere stays static to machine precision, at both orders.
3. Sedov–Taylor blast: shock radius within 4% of 1.15167 (E t²/ρ)^{1/5}, and total energy
   conserved to 1e−6.
4. Wind injection: the change in total mass and energy matches Ṁ t and L t to 1e−8.
