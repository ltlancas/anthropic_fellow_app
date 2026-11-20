# Code documentation

All commands are run from `opus_high/`, using the standalone environment `venv/` (Python 3.9, numpy 2.0.2, scipy 1.13.1).

```
python3 -m venv venv && venv/bin/pip install numpy scipy   # one-time setup
venv/bin/python run_sam.py              # Task 1: writes sam_outputs/F*.csv, F*.txt
venv/bin/python tests/test_hydro.py     # hydro unit tests (~0.2 s)
venv/bin/python run_sims.py             # Task 2: writes sim_outputs/ (~5 s total)
venv/bin/python compute_convergence.py  # Task 3: writes convergence/*_L1_norm.txt
```

## Layout

| Path | Purpose |
|---|---|
| `sam/units.py` | cgs constants, `CodeUnits` (pc, kyr, Msun), $n_H\leftrightarrow\rho$, $T(P,\rho)$ with $\mu=1.4$ |
| `sam/cc85.py` | Chevalier & Clegg (1985) free-wind solution |
| `sam/shell.py` | self-similar adiabatic shell ODE integration |
| `sam/model.py` | `WindBubble` class: `state(t)` (R1, Rc, R2, P_b) and `profile(r, t)` |
| `sam_outputs/assumptions.md` | SAM derivations and assumptions |
| `hydro/solver.py` | `Hydro` class: the 1D finite-volume solver |
| `tests/test_hydro.py` | unit tests |
| `run_sims.py` | F100 runs (R_inj = 0.5 pc), fo/so × N = 64, 128, 256 |
| `compute_convergence.py` | L1 comparison with the SAM |

## Hydro solver (`hydro/solver.py`)

- **Equations.** The 1D Euler equations with ideal-gas $\gamma$, in `geometry="spherical"` (default) or `"cartesian"`.
  The finite-volume update in spherical coordinates is
  $\frac{dU_i}{dt} = -\frac{A_{i+1/2}F_{i+1/2} - A_{i-1/2}F_{i-1/2}}{V_i} + S_i$.
  Face areas are $A=4\pi r^2$ and exact shell volumes are $V_i = \frac{4\pi}{3}(r_{i+1/2}^3-r_{i-1/2}^3)$. The geometric
  momentum source is $P_i (A_{i+1/2}-A_{i-1/2})/V_i$, which preserves a uniform static state to machine precision.
- **Riemann solver.** HLLC with Davis wave-speed estimates.
- **First order** (`order=1`): piecewise-constant states and forward-Euler time stepping.
- **Second order** (`order=2`): piecewise-linear reconstruction of the primitive variables $(\rho, v, P)$ with minmod-limited slopes, and
  SSP-RK2 (Heun) time stepping.
- **Sources.** Uniform mass ($\dot M/V_{\rm inj}$) and energy ($L/V_{\rm inj}$) injection inside `R_inj`. A cell that straddles
  `R_inj` receives the exact fraction of its volume that lies inside. The injected mass carries no momentum.
- **Boundaries.** Two ghost cells per side. Inner: `reflect` (default, for r = 0) or `outflow`. Outer: `outflow`
  (zero-gradient) or `reflect`.
- **Timestep.** `cfl * min(dr / (|v| + c))` (default cfl = 0.4). It is also capped so that the injected energy per step is at most
  30% of a cell's energy. `evolve()` shortens the last step before each output so that outputs fall exactly on the requested times.
- **Floors.** `rho_floor` and `P_floor` are applied when converting conserved to primitive variables.
- **Vectorisation.** There are no loops over cells. All reconstruction, Riemann-solver and update operations are numpy array
  operations on `(3, N)` arrays.

### Usage
```python
from hydro.solver import Hydro
h = Hydro(N=128, r_max=10.0, order=2, geometry="spherical", Mdot=..., L=..., R_inj=0.5)
h.set_primitive(rho0, 0.0, P0)
h.evolve(t_end, output_times=[...], callback=lambda h, k: ...)
rho, v, P = h.W          # primitives at cell centres h.r
```

## Simulation setup (`run_sims.py`)

- Code units: [L] = pc, [t] = kyr, [M] = M$_\odot$ (`sam.units.CODE`).
- Domain: $0 \le r \le 10$ pc. The SAM gives $R_2(100\,{\rm kyr}) = 7.62$ pc, so the shock stays inside the domain (about $1.3R_2$) at all
  output times.
- Initial conditions: uniform $n_H=100$ cm$^{-3}$ ($\rho = 1.4 m_p n_H$), $T = 100$ K, $v = 0$.
- Wind: $\dot M = 10^{-4}$ M$_\odot$/yr and $L = 10^{38}$ erg/s, injected uniformly within 0.5 pc.
- Snapshots at 20, 40, 60, 80 and 100 kyr (`snap001`–`snap005`). The CSV columns are cell-centre $r$ [pc],
  $\rho$ [g cm$^{-3}$], $v$ [km/s], $P/k_B$ [K cm$^{-3}$], and $T$ [K] (computed with $\mu = 1.4$).

## Tests (`tests/test_hydro.py`, < 1 s)

1. **Sod shock tube** (cartesian, first and second order). The star-region values match the exact solution (within 3% at first order and 1% at second order),
   and the shock position is within 0.01.
2. **Accuracy.** The second-order error in the post-shock plateau is smaller than the first-order error.
3. **Spherical well-balancing.** A uniform static state stays static to $10^{-12}$.
4. **Injection conservation.** Total mass and energy grow exactly as $\dot M t$ and $L t$, to a relative precision of $10^{-10}$.

## Task 3 (`compute_convergence.py`)

The SAM is evaluated with `R_inj = 0.5 pc` at the cell centres of snapshot 005. The relative errors are averaged over cells with
$r < R_2^{\rm SAM}(100\,{\rm kyr}) = 7.62$ pc.
