# Code documentation

Environment: `venv/` (numpy, scipy). Run everything with `./venv/bin/python` from this directory.

- `sam/__init__.py` — semi-analytic model (`Bubble(Mdot[Msun/yr], L[erg/s], nH[cm^-3])`; `R2(t)`, `P_hot(t)`, `profile(r, t)` in cgs). See `sam_outputs/assumptions.md`.
- `run_sam.py` — writes `sam_outputs/F*.csv|txt`.
- `hydro/__init__.py` — `Hydro(N, rmax, rho0, P0, order=1|2, Mdot, L, r_inj)`: 1D spherical finite-volume solver on a uniform grid from r=0 to rmax. Conserved variables (rho, rho v, E); HLLE Riemann solver; face-area-weighted flux divergence with geometric pressure source P dA/V. Order 1: piecewise-constant + forward Euler. Order 2: minmod-limited linear reconstruction of primitives + SSP-RK2 (Heun). Inner boundary reflecting, outer outflow. Mass and energy (thermal) injected uniformly in r < r_inj with volume-overlap weights. CFL 0.4. Methods: `evolve_to(t)`, `prim(U)`.
- `test_hydro.py` — fast tests (uniform state preservation, exact mass/energy injection budget).
- `run_sims.py` — F100 runs (r_inj = 0.5 pc, ambient T = 100 K, domain 1.5 x SAM R2(100 kyr)), N = 64/128/256, both orders; writes `sim_outputs/`.
- `convergence.py` — L1 norms vs SAM at the last snapshot, written to `convergence/`.
