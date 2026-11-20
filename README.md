# Wind-blown bubble task: prompt and grader

- `prompt.md`: the task given to the model being evaluated. It covers a semi-analytic model (SAM), a 1D spherical hydro code, and a convergence study.
- `grader.md`: the written grading specification.
- `grader_tools/grade.py`: a standalone Python implementation of `grader.md`.
- `grader_tools/ad_wind_spec.py`: the reference solution (`AdiabaticWind`, Weaver et al. 1977, with a Chevalier & Clegg 1985 free wind).
- `grader_tools/coding_rubric.md`: the rubric for the manual code-quality grade.
- `exp_*/`: outputs from individual runs of `prompt.md`.

## Setup

The grader needs Python ≥ 3.11 with `numpy>=2.3`, `scipy>=1.16` and `astropy>=7.1`. These are the versions in `grader_tools/pyproject.toml`; note that numpy 2.3 needs Python 3.11 or newer. Set up a dedicated environment:

```bash
python3.13 -m venv grader_tools/.venv
grader_tools/.venv/bin/pip install "numpy>=2.3" "scipy>=1.16" "astropy>=7.1"
```

## Running the grader

```bash
grader_tools/.venv/bin/python grader_tools/grade.py exp_3
```

To grade several directories in one go, list them all: `grade.py exp_1 exp_2 exp_3`. A full run takes a few seconds.

For each directory `<ref_dir>`, the grader writes these files to `<ref_dir>/grades/`:

| File | Contents |
| --- | --- |
| `F100_sam_grade.txt`, `F10_sam_grade.txt`, `F1_sam_grade.txt` | SAM grades for outer radius, bubble pressure, and the density, velocity and pressure profiles |
| `sim_grade_fo.txt`, `sim_grade_so.txt` | Grades for density, pressure and velocity in each simulation snapshot, one file per order |
| `convergence_grades.txt` | 2 × 4 grades (first and second order × T, ρ, v, P) for how close the fitted convergence index comes to −1 and −2 |
| `grade_summary.md` | All grades, the per-section totals, the final grade (their sum), and a list of format problems and missing outputs |

Re-running the grader overwrites these files. It does not touch `code_quality.txt`.

## How grades are computed

The grader follows `grader.md`.

- **Grade formula.** Every grade is `min(max(-log10 x, 0), 6)`, where `x` is a relative-error L1 norm `|q - q_ref| / |q_ref|` or a difference in power-law index. Higher is better. Missing or unreadable outputs score 0.
- **Profile comparisons** average over grid points with `0 < r < AdiabaticWind.radius(t)`. Points where `q_ref = 0` are skipped.
- **Reference models.** The SAM is compared against `AdiabaticWind` with `rfb = 100 AU`, and the simulations against `rfb = 0.5 pc` (n_H = 100 cm⁻³).
- **Snapshots.** A snapshot whose time is `t <= 0` is an initial condition: it is skipped and not counted as missing.
- **Missing outputs.** Each order is expected to have 3 resolutions × 5 snapshots. Any that are missing score 0. Extra resolutions or snapshots are graded and added to the total.

**Grading is strict.** Outputs must match the "Output File Formats" section of `prompt.md` exactly. That section and the readers in `grade.py` are kept in sync, so change both together. A file that doesn't match is treated as missing, and the reason is listed in `grade_summary.md`. In short:

- **`.csv` files:** comma-separated. The first line must be exactly `r,rho,v,P,T`: no leading `#`, no extra columns, and no other header lines. Units are pc, g cm⁻³, km s⁻¹, K cm⁻³ and K.
- **`sam_outputs/<name>.txt`:** a single line with two space-separated numbers, the outer radius in pc and P/k_B in K cm⁻³.
- **`sim_outputs/F100_sim_NXXXXX_times.txt`:** one number per line, where line k is the time of snapshot k in kyr.
- **`convergence/<variable>_L1_norm.txt`:** `<variable>` is one of `density`, `temperature`, `pressure` or `velocity`. Each line holds three space-separated numbers: N (number of cells), the first-order L1 and the second-order L1.
- **No headers, comments or labels** in any `.txt` file. Blank lines are ignored.
- **File names** zero-pad N to 5 digits and the snapshot number to 3 digits, starting at `001` (e.g. `F100_sim_so_N01024_snap003.csv`). Files with other names are ignored, and the summary lists them.

Runs made with older versions of `prompt.md` (exp_1 to exp_3) don't follow this format, so they score low or zero for format reasons alone.

## Code quality (manual step)

The code-quality task in `grader.md` needs a human or LLM judgement, so the script doesn't do it. To include it:

1. Grade the SAM and simulation code in `<ref_dir>/` using `grader_tools/coding_rubric.md`.
2. Write the score as the first number in `<ref_dir>/grades/code_quality.txt`. Lines starting with `#` are ignored, so you can put the rubric's breakdown there as comments.
3. Re-run `grade.py`. The summary will include the score and add it to the final grade.
