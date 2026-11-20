# Wind Blown Bubble Grader

This prompt should have been provided with a reference directory (which will be referred to below as `ref_dir`), create a `grades/` subdirectory of that directory.

## Units and Comparison Conventions (apply to every section below)

**Units.** The `.csv` files in `ref_dir/` were requested with a header row and the following columns and units:

| Column name | Quantity | Units |
| --- | --- | --- |
| `r` | radius | ${\rm pc}$ |
| `rho` | mass density | ${\rm g}\,{\rm cm}^{-3}$ |
| `v` | radial velocity | ${\rm km}\,{\rm s}^{-1}$ |
| `P` | pressure divided by Boltzmann's constant ($P/k_B$) | ${\rm K}\,{\rm cm}^{-3}$ |
| `T` | temperature | ${\rm K}$ |

The `.txt` files in `sam_outputs/` give the outer bubble radius in ${\rm pc}$ and the hot-gas pressure $P/k_B$ in ${\rm K}\,{\rm cm}^{-3}$. The `sim_outputs/F100_sim_NXXXXX_times.txt` files give snapshot times in ${\rm kyr}$, one per line in snapshot order.

`AdiabaticWind()` works with `astropy` quantities. Construct its inputs with units attached: `Lwind` in erg/s, `Mdotw` in Msun/yr, `rfb` as a length, and `rho0` $= 1.4\,m_p\,n_{\rm H}$ as a mass density. Pass `r` in pc and `t` in kyr. Before comparing, convert every reference output to the units in the table above:
- `radius(t)` to `pc`
- `pressure(t)` to `K/cm3`
- `density_profile(r, t)` to `g/cm3` (it returns `solMass/pc3` by default)
- `velocity_profile(r, t)` to `km/s`
- `pressure_profile(r, t)` to `K/cm3`

**Relative errors.** Every comparison against `AdiabaticWind()` uses relative errors, $|q - q_{\rm ref}|/|q_{\rm ref}|$, where $q_{\rm ref}$ is the reference value.
- For the single numbers in the `.txt` files this is the score's L1 norm directly.
- For profiles, the "average L1 norm" is the mean of the relative error over the radii in the file that are smaller than the reference outer bubble radius, `AdiabaticWind.radius(t)`, at the relevant time. Radii at or beyond it are excluded, as are any points where $q_{\rm ref} = 0$.

**Zero radii and times.** The reference model is singular at $r = 0$ and $t = 0$ (its methods raise errors or return `nan`), so exclude these before calling it:
- Drop any rows with $r \le 0$ (e.g. a grid starting at the origin, or ghost cells) from every profile comparison.
- Skip any snapshot with $t \le 0$ entirely (e.g. an initial-condition snapshot); it gets no grade and does not count as missing.

**Logarithms.** Every grade written to a file in `grades/` as minus the (base-10) logarithm of a quantity $x$ (an L1 norm, or the absolute difference between a fitted and expected power-law index) is $\min\left(\max\left(-\log_{10} x,\ 0\right),\ 6\right)$. Clamping this way means every grade is non-negative, a higher grade is always better, and an exact match ($x = 0$) gives the maximum of 6 instead of infinity. Missing outputs score 0, the same as the worst possible grade.

## Semi-Analytic Model Task

### Specification for Grader

A previous instance of Claude has attempted to create a semi-analytic model for the evolution of wind-blown bubbles as a function of time and radius under the assumption of constant injection of mass ($\dot{M}_w$) and energy ($L_w$) within a given radius ($r_{\rm fb}$), into a background with a given constant number density of hydrogen nuclei ($n_{\rm H}$) and mean molecular weight of the gas (per hydrogen nucleus, $\mu_{\rm H}$). Below I list a set of experiments that the model should have carried out with corresponding names, parameters (wind mass loss, energy input, and injection radius), and the time at which it was asked to output the profiles.

| Experiment Name | $\dot{M}_w$ [$M_{\odot}/{\rm yr}$] | $L_w$ [${\rm erg}\,{\rm s}^{-1}$]  | $n_{\rm H}$ [${\rm cm}^{-3}$] | $t_{\rm out}$ [${\rm kyr}$]
| --- | --- | --- | --- | --- |
| F100 | $10^{-4}$ | $10^{38}$ | $100$ | $100$ |
| F10 | $10^{-4}$ | $10^{38}$ | $10$ | $100$ |
| F1 | $10^{-4}$ | $10^{38}$ | $1$ | $100$ |

For each of the experiments the injection radius is assumed to be $100\, {\rm AU}$ and the mean molecular weight of the gas is assumed to be $1.4 m_p$ where $m_p$ is the proton mass.

### How to Grade

For each of the models above the previous version of Claude should have output the model over a grid of radii at the times specified above in `ref_dir/` and under a sub-directory `sam_outputs/`. If this sub-directory does not exist, assign a score of zero for all of the grades in this section. For each of the models there should be two output files, the first (`sam_outputs/<experiment_name>.csv`) should specify the fluid variables over the radius range specified above. The second (`sam_outputs/<experiment_name>.txt`) should specify two numbers, the radius of the bubble in parsecs and the pressure of the hot gas within the bubble (divided by Boltzmann's constant) in units of Kelvin per cubic centimeter.

In the `grader_tools/` directory there is a python file (`ad_wind_spec.py`) that specifies a class `AdiabaticWind()` that takes as arguments

| Argument | Meaning |
| --- | --- |
| `Lwind` | $L_w$ or enery input rate |
| `Mdotw` | $\dot{M}w$ or mass loss rate|
| `rfb` | $r_{\rm fb}$ or injection radius |
| `rho0` | $\rho$ or background mass density |

It has methods `radius`, `pressure`, `density_profile`, `velocity_profile`, and `pressure_profile`. For each of the models above use the first two methods to compare against the numbers contained in the `sam_outputs/<experiment_name>.txt` (a relative-error L1 norm comparison) and the last three methods to compare against the density, velocity, and pressure in the `sam_outputs/<experiment_name>.csv` at each of the radii smaller than the reference outer bubble radius (an average relative-error L1 norm over those radii). Use the units and conventions above. Save minus the base-10 logarithm of each of these in a text file in `ref_dir/grades/<model_name>_sam_grade.txt` with a header explaining which test each number refers to.

## Hydrodynamics Code Task

### Specification for Grader

There should also be a directory `ref_dir/sim_outputs/`, if this directory does not exist, assign a score of zero for each of the grades in this section. This directory should contain the outputs in a `sim_outputs/F100_sim_VV_NXXXXX_snapYYY.csv` where `VV` is either `fo` or `so` for first or second order, `XXXXX` specifies the number of resolution elements in the simulation, `YYY` specifies the snapshot number, and the file contains all primitive fluid variables and the radial position of each grid point. There should also be a short `sim_outputs/F100_sim_NXXXXX_times.txt` text file that contains the times of each output.

### How to Grade

For each snapshot, at each order, and at each resolution use the `AdiabaticWind()` module in `ad_wind_spec.py` (with `rfb` as 0.5 pc) to compute the average relative-error L1 norm of the density, pressure, and velocity using the `density_profile`, `pressure_profile`, and `velocity_profile` methods in comparison to those values provided in the `sim_outputs/F100_sim_VV_NXXXXX_snapYYY.csv` files, averaging only over grid points with radii smaller than the reference outer bubble radius at that snapshot's time. Use the units and conventions above. Use the `sim_outputs/F100_sim_NXXXXX_times.txt` (times in kyr) to reference the relevant time for each comparison. Save minus the base-10 logarithm of the average L1 norms for each snapshot and resolution in a text file `ref_dir/grades/sim_grade_VV.txt`. where `VV` stands for `fo` or `so` for the first and second order results.


## Agreement Task

The final task was to do a convergence study to check for agreement between the semi-analytic model and the code base. In `ref_dir/convergence/` (if this directory doesn't exist, assign a score of zero for all of the grades in this section) there should be average relative-error L1 norms (averaged over radii smaller than the outer bubble radius) for temperature, density, velocity, and pressure as a function of resolution. Check that the simulations converge at the proper order for each fluid variable by fitting a power-law to the L1 norms for the first and second order simulations. Compare the resulting power law against -1 and -2 for first and second order resepectively. Place minus the log base 10 of the differences in `ref_dir/grades/convergence_grades.txt` which has 2 rows (one per order) × 4 columns (one per variable).

## Code Quality Task

Read the `coding_rubric.md` file and evaluate a score for the code created in `ref_dir/` for both the Semi-Analytic Model code and the simulation code. Place this score in `ref_dir/grades/code_quality.txt`.

## Summary

Write a summary `ref_dir/grades/grade_summary.md` file that summarizes these grades and provides a single number which is the sum of these grades as a final grade.