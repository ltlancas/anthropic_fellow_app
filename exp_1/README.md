# wind_bubble

This is a semi-analytic, spherically symmetric model of a wind-blown bubble:
constant injection of mass and energy at the origin into a uniform medium. It
returns density, pressure, temperature and velocity as functions of radius and
time. It extends Weaver et al. (1977) to finite reverse-shock radius and
mass-loaded winds. See [assumptions.md](assumptions.md).

## Setup

The environment is standalone and used only for this project:

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[plot,test]"
```

## Usage

```python
from wind_bubble import WindBubble, units as u

t = 1e5 * u.yr
wb = WindBubble(L_w=1e38 * u.erg_s, Mdot_w=1e-2 * u.Msun / u.yr,
                n_H=100 * u.per_cm3, mu_H=1.4).solve(t_end=t)
prof = wb.profiles(t, r_max_factor=1.5)   # dict of CGS arrays: r, rho, n_H, P, T, v, region
print(prof["state"].R2 / u.pc)            # forward-shock radius in pc
hist = wb.history([1e3 * u.yr, 1e4 * u.yr, t])  # global quantities vs time
```

## Fiducial run (t = 1e5 yr)

```bash
.venv/bin/python scripts/run_fiducial.py
.venv/bin/pytest -q
```

The run writes these outputs:

- `output/profiles_t1e5yr.csv`: r, ρ, n_H, P, P/k, T and v, from 0 to 1.5 R2.
- `output/profiles_t1e5yr.png`: plots of the same profiles.
- `output/summary_t1e5yr.json`: radii, pressure, masses, the ratio to the
  Weaver et al. (1977) solution, and diagnostics.

## Layout

- `wind_bubble/units.py`: CGS constants and unit conversions.
- `wind_bubble/model.py`: the `WindBubble` model, covering the global ODEs and
  the Lagrangian reconstruction of the interior.
- `wind_bubble/weaver.py`: the Weaver et al. (1977) self-similar reference solution.
- `tests/test_model.py`: tests for the Weaver limit, conservation laws, jump
  conditions and insensitivity to initial conditions.
